"""Held-out likelihood comparison for matched graph alternatives."""

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from math import exp, log, log1p

from anvaya.bidirected import BidirectedDeBruijnGraph, _successor
from anvaya.damage_likelihood import fit_candidate_damage_model
from anvaya.damage_profile import DamageProfile, DamageProfileLocus
from anvaya.tip_matching import TipBackboneMatch


_EPSILON = 1e-12
_FOLD_COUNT = 5
_MAX_SINGLETON_ERROR_QUALITY = 20
LocusKey = tuple[int, int, str]


@dataclass(slots=True, frozen=True)
class NucleotideObservation:
    """One reference or alternative call with cycle and Phred quality."""

    alternative: bool
    distance: int
    quality: int
    damage_probability: float | None = None
    damage_probability_ensemble: tuple[float, ...] = ()
    molecule_id: int | None = None


@dataclass(slots=True, frozen=True)
class EventLikelihood:
    """Comparable predictive scores for three explanations of an event."""

    status: str
    reasons: tuple[str, ...]
    profile_scope: str
    conditioning_scope: str
    observations: int
    alternative_observations: int
    reference_observations: int
    missing_quality_observations: int
    damage_log_likelihood: float | None
    damage_log_likelihood_min: float | None
    damage_log_likelihood_max: float | None
    damage_log_likelihood_spread: float | None
    damage_conditioning_log_probability: float | None
    error_log_likelihood: float | None
    error_conditioning_log_probability: float | None
    variation_log_likelihood: float | None
    variation_conditioning_log_probability: float | None
    variation_frequency: float | None
    variation_complexity_penalty: float | None
    variation_penalized_log_likelihood: float | None
    best_explanation: str | None
    log_likelihood_margin: float | None
    damage_error_log_contrast: float | None
    damage_variation_log_contrast: float | None
    error_variation_log_contrast: float | None


@dataclass(slots=True, frozen=True)
class CrossFittedDamageModels:
    """Damage curves fitted without the candidate loci they score."""

    status: str
    reasons: tuple[str, ...]
    folds: int
    fitted_folds: int
    locus_folds: dict[LocusKey, int]
    probabilities_by_fold: dict[int, tuple[float, ...]]
    independent_probabilities: tuple[float, ...]

    def probabilities_for(
        self,
        key: LocusKey,
    ) -> tuple[tuple[float, ...] | None, str]:
        """Return held-out probabilities and their observation scope."""
        fold = self.locus_folds.get(key)
        if fold is None:
            return (
                self.independent_probabilities or None,
                "independent_tip_candidate_profile",
            )
        return (
            self.probabilities_by_fold.get(fold),
            f"held_out_candidate_fold_{fold + 1}_of_{self.folds}",
        )


def fit_cross_fitted_damage_models(
    profile: DamageProfile,
    folds: int = _FOLD_COUNT,
) -> CrossFittedDamageModels:
    """Fit deterministic folds so a candidate never trains its own curve."""
    if not isinstance(folds, int) or isinstance(folds, bool):
        raise TypeError("folds must be an integer")
    if folds < 2:
        raise ValueError("folds must be at least 2")

    ordered = sorted(
        profile.loci,
        key=lambda locus: (
            locus.alternative_edge_id,
            locus.reference_edge_id,
            locus.channel,
        ),
    )
    locus_folds = {
        (
            locus.alternative_edge_id,
            locus.reference_edge_id,
            locus.channel,
        ): index % folds
        for index, locus in enumerate(ordered)
    }
    probabilities_by_fold = {}
    reasons = []
    for fold in range(folds):
        training = [
            (locus.alternative, locus.reference)
            for locus in ordered
            if locus_folds[
                (
                    locus.alternative_edge_id,
                    locus.reference_edge_id,
                    locus.channel,
                )
            ]
            != fold
        ]
        fit = fit_candidate_damage_model(training, profile.end_window)
        if fit.status == "fitted":
            probabilities_by_fold[fold] = fit.fitted_probabilities
        else:
            reasons.append(f"fold_{fold + 1}_insufficient_evidence")

    full_fit = profile.candidate_damage_fit
    independent = (
        full_fit.fitted_probabilities
        if full_fit.status == "fitted"
        else ()
    )
    if not independent:
        reasons.append("full_candidate_profile_insufficient_evidence")
    fitted_folds = len(probabilities_by_fold)
    return CrossFittedDamageModels(
        status=(
            "fitted"
            if fitted_folds == folds and independent
            else "insufficient_evidence"
        ),
        reasons=tuple(reasons),
        folds=folds,
        fitted_folds=fitted_folds,
        locus_folds=locus_folds,
        probabilities_by_fold=probabilities_by_fold,
        independent_probabilities=independent,
    )


def _emission_probability(
    alternative: bool,
    latent_alternative_probability: float,
    quality: int,
) -> float:
    error = 10.0 ** (-quality / 10.0)
    if alternative:
        probability = (
            latent_alternative_probability * (1.0 - error)
            + (1.0 - latent_alternative_probability) * error / 3.0
        )
    else:
        probability = (
            (1.0 - latent_alternative_probability) * (1.0 - error)
            + latent_alternative_probability * error / 3.0
        )
    return min(1.0 - _EPSILON, max(_EPSILON, probability))


def _logaddexp(left: float, right: float) -> float:
    maximum = max(left, right)
    return maximum + log(exp(left - maximum) + exp(right - maximum))


def _observation_groups(
    observations: Sequence[NucleotideObservation],
) -> tuple[tuple[NucleotideObservation, ...], ...]:
    groups: dict[tuple[str, object], list[NucleotideObservation]] = {}
    for index, observation in enumerate(observations):
        key = (
            ("molecule", observation.molecule_id)
            if observation.molecule_id is not None
            else ("observation", index)
        )
        groups.setdefault(key, []).append(observation)
    return tuple(tuple(group) for group in groups.values())


def _conditioned_log_likelihood(
    observations: Sequence[NucleotideObservation],
    probabilities: Sequence[float],
    *,
    use_observation_damage: bool = False,
    ensemble_index: int | None = None,
    groups: Sequence[Sequence[NucleotideObservation]] | None = None,
) -> tuple[float, float]:
    """Condition a two-path likelihood on consistent, observed alleles."""
    observed_log_likelihood = 0.0
    consistency_log_probability = 0.0
    log_all_alternative = 0.0
    log_all_reference = 0.0
    if groups is None:
        groups = _observation_groups(observations)
    for group in groups:
        alternative_log_probability = 0.0
        reference_log_probability = 0.0
        for observation in group:
            if ensemble_index is not None:
                latent_probability = observation.damage_probability_ensemble[
                    ensemble_index
                ]
            else:
                latent_probability = (
                    observation.damage_probability
                    if use_observation_damage
                    and observation.damage_probability is not None
                    else probabilities[observation.distance]
                )
            alternative_log_probability += log(
                _emission_probability(
                    True,
                    latent_probability,
                    observation.quality,
                )
            )
            reference_log_probability += log(
                _emission_probability(
                    False,
                    latent_probability,
                    observation.quality,
                )
            )
        normalization = _logaddexp(
            alternative_log_probability,
            reference_log_probability,
        )
        consistency_log_probability += normalization
        log_all_alternative += alternative_log_probability - normalization
        log_all_reference += reference_log_probability - normalization
        observed_log_likelihood += (
            alternative_log_probability
            if group[0].alternative
            else reference_log_probability
        )

    excluded_probability = exp(log_all_alternative) + exp(log_all_reference)
    selection_log_probability = log1p(
        -min(1.0 - _EPSILON, excluded_probability)
    )
    conditioning_log_probability = (
        consistency_log_probability + selection_log_probability
    )
    return (
        observed_log_likelihood - conditioning_log_probability,
        conditioning_log_probability,
    )


def _conditioned_constant_log_likelihood(
    group_counts: Counter[tuple[bool, tuple[int, ...]]],
    frequency: float,
) -> tuple[float, float]:
    """Evaluate a constant-frequency model over compressed fragments."""
    observed_log_likelihood = 0.0
    consistency_log_probability = 0.0
    log_all_alternative = 0.0
    log_all_reference = 0.0
    for (observed_alternative, qualities), count in group_counts.items():
        alternative_log_probability = sum(
            log(_emission_probability(True, frequency, quality))
            for quality in qualities
        )
        reference_log_probability = sum(
            log(_emission_probability(False, frequency, quality))
            for quality in qualities
        )
        normalization = _logaddexp(
            alternative_log_probability,
            reference_log_probability,
        )
        consistency_log_probability += count * normalization
        log_all_alternative += count * (
            alternative_log_probability - normalization
        )
        log_all_reference += count * (
            reference_log_probability - normalization
        )
        observed_log_likelihood += count * (
            alternative_log_probability
            if observed_alternative
            else reference_log_probability
        )
    excluded_probability = exp(log_all_alternative) + exp(log_all_reference)
    selection_log_probability = log1p(
        -min(1.0 - _EPSILON, excluded_probability)
    )
    conditioning_log_probability = (
        consistency_log_probability + selection_log_probability
    )
    return (
        observed_log_likelihood - conditioning_log_probability,
        conditioning_log_probability,
    )


def _ensemble_damage_log_likelihoods(
    observations: Sequence[NucleotideObservation],
    groups: Sequence[Sequence[NucleotideObservation]],
) -> tuple[tuple[float, float], ...]:
    """Score aligned cross-fit curves without treating folds as molecules."""
    ensemble_size = min(
        (len(observation.damage_probability_ensemble) for observation in observations),
        default=0,
    )
    return tuple(
        _conditioned_log_likelihood(
            observations,
            (),
            ensemble_index=index,
            groups=groups,
        )
        for index in range(ensemble_size)
    )


def _maximize_variation_frequency(
    observations: Sequence[NucleotideObservation],
    groups: Sequence[Sequence[NucleotideObservation]],
) -> tuple[float, float, float]:
    """Fit one constant alternative frequency by bounded golden search."""
    group_counts = Counter(
        (
            group[0].alternative,
            tuple(sorted(observation.quality for observation in group)),
        )
        for group in groups
    )
    lower = _EPSILON
    upper = 0.5
    ratio = (5.0**0.5 - 1.0) / 2.0
    left = upper - ratio * (upper - lower)
    right = lower + ratio * (upper - lower)

    def objective(frequency: float) -> float:
        likelihood, _ = _conditioned_constant_log_likelihood(
            group_counts,
            frequency,
        )
        return likelihood

    left_score = objective(left)
    right_score = objective(right)
    # Forty iterations resolve a frequency interval below 1e-8, far beyond
    # the six decimal places written to reports.
    for _ in range(40):
        if left_score > right_score:
            upper = right
            right = left
            right_score = left_score
            left = upper - ratio * (upper - lower)
            left_score = objective(left)
        else:
            lower = left
            left = right
            left_score = right_score
            right = lower + ratio * (upper - lower)
            right_score = objective(right)
    frequency = (lower + upper) / 2.0
    likelihood, conditioning = _conditioned_constant_log_likelihood(
        group_counts,
        frequency,
    )
    return frequency, likelihood, conditioning


def compare_event_likelihoods(
    observations: Sequence[NucleotideObservation],
    damage_probabilities: Sequence[float],
    *,
    profile_scope: str = "held_out_candidate_profile",
    missing_quality_observations: int = 0,
    damage_explanation_applicable: bool = True,
) -> EventLikelihood:
    """Compare damage, error, and constant-frequency variation models."""
    reasons = []
    keyed_observations = [
        (
            (
                ("molecule", observation.molecule_id)
                if observation.molecule_id is not None
                else ("observation", index)
            ),
            observation,
        )
        for index, observation in enumerate(observations)
    ]
    alleles_by_molecule: dict[tuple[str, object], set[bool]] = {}
    for key, observation in keyed_observations:
        alleles_by_molecule.setdefault(key, set()).add(observation.alternative)
    conflicting_molecules = {
        key for key, alleles in alleles_by_molecule.items() if len(alleles) > 1
    }
    keyed_observations = [
        item for item in keyed_observations if item[0] not in conflicting_molecules
    ]
    molecule_keys = [key for key, _ in keyed_observations]
    observations = [observation for _, observation in keyed_observations]
    alternatives = len(
        {
            key
            for key, observation in zip(molecule_keys, observations, strict=True)
            if observation.alternative
        }
    )
    references = len(
        {
            key
            for key, observation in zip(molecule_keys, observations, strict=True)
            if not observation.alternative
        }
    )
    molecule_count = len(set(molecule_keys))
    if not observations:
        reasons.append("no_eligible_terminal_observations")
    if alternatives == 0:
        reasons.append("no_alternative_observations")
    if references == 0:
        reasons.append("no_reference_observations")
    if missing_quality_observations:
        reasons.append("missing_base_quality")
    has_embedded_damage = bool(observations) and all(
        observation.damage_probability is not None
        for observation in observations
    )
    if not damage_probabilities and not has_embedded_damage:
        reasons.append("damage_profile_unavailable")
    if observations and damage_probabilities and not has_embedded_damage and any(
        observation.distance >= len(damage_probabilities)
        for observation in observations
    ):
        reasons.append("damage_profile_window_mismatch")
    if reasons:
        return EventLikelihood(
            status="insufficient_evidence",
            reasons=tuple(reasons),
            profile_scope=profile_scope,
            conditioning_scope="two_path_fragment_ascertainment",
            observations=molecule_count,
            alternative_observations=alternatives,
            reference_observations=references,
            missing_quality_observations=missing_quality_observations,
            damage_log_likelihood=None,
            damage_log_likelihood_min=None,
            damage_log_likelihood_max=None,
            damage_log_likelihood_spread=None,
            damage_conditioning_log_probability=None,
            error_log_likelihood=None,
            error_conditioning_log_probability=None,
            variation_log_likelihood=None,
            variation_conditioning_log_probability=None,
            variation_frequency=None,
            variation_complexity_penalty=None,
            variation_penalized_log_likelihood=None,
            best_explanation=None,
            log_likelihood_margin=None,
            damage_error_log_contrast=None,
            damage_variation_log_contrast=None,
            error_variation_log_contrast=None,
        )

    observation_groups = _observation_groups(observations)
    damage, damage_conditioning = _conditioned_log_likelihood(
        observations,
        damage_probabilities,
        use_observation_damage=has_embedded_damage,
        groups=observation_groups,
    )
    ensemble_damage = _ensemble_damage_log_likelihoods(
        observations,
        observation_groups,
    )
    damage_candidates = (damage, *(score for score, _ in ensemble_damage))
    damage_min = min(damage_candidates)
    damage_max = max(damage_candidates)
    error, error_conditioning = _conditioned_log_likelihood(
        observations,
        (0.0,) * (max(observation.distance for observation in observations) + 1),
        groups=observation_groups,
    )
    variation_frequency, variation, variation_conditioning = (
        _maximize_variation_frequency(observations, observation_groups)
    )
    variation_penalty = 0.5 * log(len(observations))
    variation_penalized = variation - variation_penalty
    scores = {
        "sequencing_error": error,
        "variation": variation_penalized,
    }
    if damage_explanation_applicable:
        scores["damage"] = damage
    ranked = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    best_explanation = ranked[0][0]
    ranking_reasons = ["conditioned_on_two_path_fragment_ascertainment"]
    if (
        best_explanation == "sequencing_error"
        and alternatives == 1
        and min(
            observation.quality
            for observation in observations
            if observation.alternative
        ) > _MAX_SINGLETON_ERROR_QUALITY
    ):
        best_explanation = "ambiguous"
        ranking_reasons.append("high_quality_singleton_cannot_exclude_variation")
    ranking_reasons.append("ranking_not_calibrated_for_classification")
    return EventLikelihood(
        status="scored",
        reasons=tuple(ranking_reasons),
        profile_scope=profile_scope,
        conditioning_scope="two_path_fragment_ascertainment",
        observations=molecule_count,
        alternative_observations=alternatives,
        reference_observations=references,
        missing_quality_observations=0,
        damage_log_likelihood=damage,
        damage_log_likelihood_min=damage_min,
        damage_log_likelihood_max=damage_max,
        damage_log_likelihood_spread=damage_max - damage_min,
        damage_conditioning_log_probability=damage_conditioning,
        error_log_likelihood=error,
        error_conditioning_log_probability=error_conditioning,
        variation_log_likelihood=variation,
        variation_conditioning_log_probability=variation_conditioning,
        variation_frequency=variation_frequency,
        variation_complexity_penalty=variation_penalty,
        variation_penalized_log_likelihood=variation_penalized,
        best_explanation=best_explanation,
        log_likelihood_margin=ranked[0][1] - ranked[1][1],
        damage_error_log_contrast=damage - error,
        damage_variation_log_contrast=damage - variation_penalized,
        error_variation_log_contrast=error - variation_penalized,
    )


def _edge_observations_at(
    graph: BidirectedDeBruijnGraph,
    start: int,
    edge_ids: tuple[int, ...],
    edge_index: int,
    expected_end: str,
    *,
    alternative: bool,
    damage_probabilities: Sequence[float],
    damage_probability_ensemble: Sequence[Sequence[float]] = (),
) -> tuple[list[NucleotideObservation], int]:
    current = start
    for index, edge_id in enumerate(edge_ids):
        if index == edge_index:
            side = expected_end
            if current != graph.edge_sources[edge_id]:
                side = "right" if side == "left" else "left"
            offset = graph.node_length if expected_end == "left" else 0
            observations_by_molecule: dict[int, NucleotideObservation] = {}
            missing_molecules: set[int] = set()
            for link in graph.molecule_end_links(edge_id):
                distance = link.distance + offset
                if link.end != side or distance >= graph.end_window:
                    continue
                if link.base_quality is None:
                    missing_molecules.add(link.molecule_id)
                else:
                    observation = NucleotideObservation(
                        alternative=alternative,
                        distance=distance,
                        quality=link.base_quality,
                        damage_probability=damage_probabilities[distance],
                        damage_probability_ensemble=tuple(
                            probabilities[distance]
                            for probabilities in damage_probability_ensemble
                            if distance < len(probabilities)
                        ),
                        molecule_id=link.molecule_id,
                    )
                    previous = observations_by_molecule.get(link.molecule_id)
                    if previous is None or (
                        observation.distance,
                        -observation.quality,
                    ) < (previous.distance, -previous.quality):
                        observations_by_molecule[link.molecule_id] = observation
            missing_molecules.difference_update(observations_by_molecule)
            return list(observations_by_molecule.values()), len(missing_molecules)
        current = _successor(graph, current, edge_id)
    return [], 0


def _ordinary_edge_observations_at(
    graph: BidirectedDeBruijnGraph,
    match: TipBackboneMatch,
    edge_index: int,
    ensemble_size: int,
) -> tuple[list[NucleotideObservation], int]:
    """Choose the terminal side with the strongest paired path evidence."""
    zero_probabilities = (0.0,) * graph.end_window
    zero_ensemble = (zero_probabilities,) * ensemble_size
    candidates = []
    for end in ("left", "right"):
        alternative, alternative_missing = _edge_observations_at(
            graph,
            match.branch_handle,
            match.tip_edge_ids,
            edge_index,
            end,
            alternative=True,
            damage_probabilities=zero_probabilities,
            damage_probability_ensemble=zero_ensemble,
        )
        reference, reference_missing = _edge_observations_at(
            graph,
            match.branch_handle,
            match.backbone_edge_ids,
            edge_index,
            end,
            alternative=False,
            damage_probabilities=zero_probabilities,
            damage_probability_ensemble=zero_ensemble,
        )
        candidates.append(
            (
                (
                    bool(alternative) and bool(reference),
                    min(len(alternative), len(reference)),
                    len(alternative) + len(reference),
                    -(alternative_missing + reference_missing),
                    end == "right",
                ),
                alternative + reference,
                alternative_missing + reference_missing,
            )
        )
    _, observations, missing = max(candidates, key=lambda candidate: candidate[0])
    return observations, missing


def score_matched_event(
    graph: BidirectedDeBruijnGraph,
    match: TipBackboneMatch,
    models: CrossFittedDamageModels | None,
) -> EventLikelihood:
    """Score one matched graph event using exact-base molecule evidence."""
    if models is None:
        return compare_event_likelihoods(
            (),
            (),
            profile_scope="damage_profile_not_requested",
        )

    observations: list[NucleotideObservation] = []
    missing = 0
    scopes = []
    ensemble = tuple(
        models.probabilities_by_fold[fold]
        for fold in sorted(models.probabilities_by_fold)
    )
    for change in match.substitutions:
        edge_index = change.position - graph.node_length - 1
        if not (
            0 <= edge_index < len(match.tip_edge_ids)
            and edge_index < len(match.backbone_edge_ids)
        ):
            continue
        if change.reference == "C" and change.alternative == "T":
            expected_end = "left"
            channel = "five_prime_ct"
        elif change.reference == "G" and change.alternative == "A":
            expected_end = "right"
            channel = "three_prime_ga"
        else:
            ordinary, ordinary_missing = _ordinary_edge_observations_at(
                graph,
                match,
                edge_index,
                len(ensemble),
            )
            observations.extend(ordinary)
            missing += ordinary_missing
            scopes.append("ordinary_substitution_no_damage_channel")
            continue
        key = (
            match.tip_edge_ids[edge_index],
            match.backbone_edge_ids[edge_index],
            channel,
        )
        probabilities, scope = models.probabilities_for(key)
        if probabilities is None:
            return compare_event_likelihoods(
                (),
                (),
                profile_scope=scope,
            )
        scopes.append(scope)
        alternative, alternative_missing = _edge_observations_at(
            graph,
            match.branch_handle,
            match.tip_edge_ids,
            edge_index,
            expected_end,
            alternative=True,
            damage_probabilities=probabilities,
            damage_probability_ensemble=ensemble,
        )
        reference, reference_missing = _edge_observations_at(
            graph,
            match.branch_handle,
            match.backbone_edge_ids,
            edge_index,
            expected_end,
            alternative=False,
            damage_probabilities=probabilities,
            damage_probability_ensemble=ensemble,
        )
        observations.extend(alternative)
        observations.extend(reference)
        missing += alternative_missing + reference_missing

    unique_scopes = tuple(dict.fromkeys(scopes))
    scope = (
        unique_scopes[0]
        if len(unique_scopes) == 1
        else "multiple_held_out_candidate_folds"
    ) if unique_scopes else "no_damage_compatible_locus"
    return compare_event_likelihoods(
        observations,
        (),
        profile_scope=scope,
        missing_quality_observations=missing,
        damage_explanation_applicable=any(
            item != "ordinary_substitution_no_damage_channel"
            for item in unique_scopes
        ),
    )
