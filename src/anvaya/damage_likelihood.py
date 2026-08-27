"""Candidate-conditioned terminal-damage likelihood fitting."""

from collections.abc import Sequence
from dataclasses import dataclass
from math import exp, lgamma, log


@dataclass(slots=True, frozen=True)
class ParameterEstimate:
    """One fitted parameter with a conditional likelihood interval."""

    value: float
    lower_95_support: float
    upper_95_support: float
    interval_method: str = "conditional_likelihood_drop_1.92"


@dataclass(slots=True, frozen=True)
class CandidateDamageFit:
    """Geometric damage fit compared with a constant-background null."""

    status: str
    reasons: tuple[str, ...]
    observation_scope: str
    formula: str
    loci: int
    locus_distance_observations: int
    observed_distances: int
    total_alternative: int
    total_reference: int
    amplitude: ParameterEstimate | None
    decay: ParameterEstimate | None
    background: ParameterEstimate | None
    concentration: ParameterEstimate | None
    fitted_probabilities: tuple[float, ...]
    log_likelihood: float | None
    null_background: float | None
    null_concentration: float | None
    null_log_likelihood: float | None
    log_likelihood_ratio: float | None
    likelihood_ratio_statistic: float | None


_EPSILON = 1e-9
_LOG_LIKELIHOOD_95_DROP = 1.920729410347062


def _beta_binomial_log_probability(
    successes: int,
    trials: int,
    probability: float,
    concentration: float,
) -> float:
    probability = min(1.0 - _EPSILON, max(_EPSILON, probability))
    alpha = probability * concentration
    beta = (1.0 - probability) * concentration
    return (
        lgamma(trials + 1)
        - lgamma(successes + 1)
        - lgamma(trials - successes + 1)
        + lgamma(successes + alpha)
        + lgamma(trials - successes + beta)
        - lgamma(trials + concentration)
        + lgamma(concentration)
        - lgamma(alpha)
        - lgamma(beta)
    )


def _coordinate_maximize(objective, starts, bounds):
    best_values = None
    best_score = float("-inf")
    initial_steps = tuple((upper - lower) / 5 for lower, upper in bounds)
    for start in starts:
        values = [
            min(upper, max(lower, value))
            for value, (lower, upper) in zip(start, bounds, strict=True)
        ]
        score = objective(values)
        steps = list(initial_steps)
        while max(steps) > 1e-5:
            improved = False
            for index, step in enumerate(steps):
                for direction in (-1.0, 1.0):
                    candidate = values.copy()
                    lower, upper = bounds[index]
                    candidate[index] = min(
                        upper,
                        max(lower, candidate[index] + direction * step),
                    )
                    candidate_score = objective(candidate)
                    if candidate_score > score + 1e-10:
                        values = candidate
                        score = candidate_score
                        improved = True
            if not improved:
                steps = [step / 2 for step in steps]
        if score > best_score:
            best_values = tuple(values)
            best_score = score
    assert best_values is not None
    return best_values, best_score


def _support_interval(
    score,
    optimum: float,
    lower: float,
    upper: float,
    best_score: float,
    *,
    logarithmic: bool = False,
) -> tuple[float, float]:
    """Return a bounded conditional likelihood-support interval."""
    if logarithmic:
        low_scale = log(lower)
        high_scale = log(upper)
        values = tuple(
            exp(low_scale + index * (high_scale - low_scale) / 200)
            for index in range(201)
        )
    else:
        values = tuple(
            lower + index * (upper - lower) / 200
            for index in range(201)
        )
    accepted = [
        value
        for value in values + (optimum,)
        if score(value) >= best_score - _LOG_LIKELIHOOD_95_DROP
    ]
    return min(accepted), max(accepted)


def fit_candidate_damage_model(
    locus_counts: Sequence[
        tuple[tuple[int, ...], tuple[int, ...]]
    ],
    end_window: int,
) -> CandidateDamageFit:
    """Fit a beta-binomial geometric curve to candidate-locus counts."""
    if not isinstance(end_window, int) or isinstance(end_window, bool):
        raise TypeError("end_window must be an integer")
    if end_window < 1:
        raise ValueError("end_window must be at least 1")
    records = []
    observed_loci = 0
    distances = set()
    total_alternative = 0
    total_reference = 0
    for alternative, reference in locus_counts:
        if len(alternative) != end_window or len(reference) != end_window:
            raise ValueError("locus count vectors must match end_window")
        if any(
            not isinstance(value, int) or isinstance(value, bool)
            for value in alternative + reference
        ):
            raise TypeError("locus counts must be integers")
        if any(value < 0 for value in alternative + reference):
            raise ValueError("locus counts must not be negative")
        locus_observed = False
        for distance, (damaged, unchanged) in enumerate(
            zip(alternative, reference, strict=True)
        ):
            trials = damaged + unchanged
            if not trials:
                continue
            records.append((distance, damaged, trials))
            distances.add(distance)
            total_alternative += damaged
            total_reference += unchanged
            locus_observed = True
        observed_loci += locus_observed

    reasons = []
    if observed_loci < 10:
        reasons.append("fewer_than_10_observed_loci")
    if len(distances) < 3:
        reasons.append("fewer_than_3_observed_distances")
    if total_alternative < 3:
        reasons.append("fewer_than_3_alternative_observations")
    if total_alternative + total_reference < 20:
        reasons.append("fewer_than_20_total_observations")
    common = dict(
        reasons=tuple(reasons),
        observation_scope="matched_candidate_loci",
        formula="p[d] = background + amplitude * decay**d",
        loci=observed_loci,
        locus_distance_observations=len(records),
        observed_distances=len(distances),
        total_alternative=total_alternative,
        total_reference=total_reference,
    )
    if reasons:
        return CandidateDamageFit(
            status="insufficient_evidence",
            amplitude=None,
            decay=None,
            background=None,
            concentration=None,
            fitted_probabilities=(),
            log_likelihood=None,
            null_background=None,
            null_concentration=None,
            null_log_likelihood=None,
            log_likelihood_ratio=None,
            likelihood_ratio_statistic=None,
            **common,
        )

    pooled = total_alternative / (total_alternative + total_reference)
    by_distance = {
        distance: sum(k for d, k, _ in records if d == distance)
        / sum(n for d, _, n in records if d == distance)
        for distance in distances
    }
    terminal = by_distance[min(distances)]
    interior = by_distance[max(distances)]
    background_guess = min(interior, pooled)
    amplitude_guess = max(0.001, terminal - background_guess)
    amplitude_fraction = min(
        0.99, amplitude_guess / max(_EPSILON, 1 - background_guess)
    )

    def alternative_objective(values):
        background, fraction, decay, log_concentration = values
        amplitude = (1 - background) * fraction
        concentration = exp(log_concentration)
        return sum(
            _beta_binomial_log_probability(
                successes,
                trials,
                background + amplitude * decay**distance,
                concentration,
            )
            for distance, successes, trials in records
        )

    alternative_bounds = (
        (_EPSILON, 0.5),
        (_EPSILON, 0.999),
        (0.01, 0.999),
        (log(1.0), log(100000.0)),
    )
    starts = (
        (background_guess, amplitude_fraction, 0.5, log(20.0)),
        (background_guess / 2, amplitude_fraction, 0.8, log(100.0)),
        (min(0.49, pooled), 0.05, 0.5, log(10.0)),
        (0.001, min(0.99, max(0.05, pooled)), 0.9, log(1000.0)),
    )
    optimum, log_likelihood = _coordinate_maximize(
        alternative_objective,
        starts,
        alternative_bounds,
    )
    background, fraction, decay, log_concentration = optimum
    amplitude = (1 - background) * fraction
    concentration = exp(log_concentration)

    def null_objective(values):
        probability, null_log_concentration = values
        null_concentration = exp(null_log_concentration)
        return sum(
            _beta_binomial_log_probability(
                successes, trials, probability, null_concentration
            )
            for _, successes, trials in records
        )

    null_bounds = ((_EPSILON, 0.999), (log(1.0), log(100000.0)))
    null_optimum, null_log_likelihood = _coordinate_maximize(
        null_objective,
        ((pooled, log(20.0)), (pooled, log(1000.0))),
        null_bounds,
    )
    null_background, null_log_concentration = null_optimum
    amplitude_interval = _support_interval(
        lambda value: alternative_objective(
            (background, value / (1 - background), decay, log_concentration)
        ),
        amplitude,
        _EPSILON,
        (1 - background) * 0.999,
        log_likelihood,
    )
    decay_interval = _support_interval(
        lambda value: alternative_objective(
            (background, fraction, value, log_concentration)
        ),
        decay,
        0.01,
        0.999,
        log_likelihood,
    )
    maximum_background = min(0.5, 1 - amplitude / 0.999)
    background_interval = _support_interval(
        lambda value: alternative_objective(
            (value, amplitude / (1 - value), decay, log_concentration)
        ),
        background,
        _EPSILON,
        maximum_background,
        log_likelihood,
    )
    concentration_interval = _support_interval(
        lambda value: alternative_objective(
            (background, fraction, decay, log(value))
        ),
        concentration,
        1.0,
        100000.0,
        log_likelihood,
        logarithmic=True,
    )
    log_ratio = max(0.0, log_likelihood - null_log_likelihood)
    return CandidateDamageFit(
        status="fitted",
        amplitude=ParameterEstimate(amplitude, *amplitude_interval),
        decay=ParameterEstimate(decay, *decay_interval),
        background=ParameterEstimate(background, *background_interval),
        concentration=ParameterEstimate(
            concentration, *concentration_interval
        ),
        fitted_probabilities=tuple(
            background + amplitude * decay**distance
            for distance in range(end_window)
        ),
        log_likelihood=log_likelihood,
        null_background=null_background,
        null_concentration=exp(null_log_concentration),
        null_log_likelihood=null_log_likelihood,
        log_likelihood_ratio=log_ratio,
        likelihood_ratio_statistic=2 * log_ratio,
        **common,
    )
