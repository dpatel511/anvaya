"""Empirical, report-only confidence calibration for graph events."""

import csv
import json
from bisect import bisect_left
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from math import isfinite, sqrt
from pathlib import Path

from anvaya.event_likelihood import EventLikelihood


EXPLANATIONS = ("damage", "sequencing_error", "variation")
CALIBRATION_FIELDS = (
    "calibration_status",
    "calibration_reasons",
    "calibration_prediction_set",
    "calibration_decision",
    "calibration_damage_pvalue",
    "calibration_error_pvalue",
    "calibration_variation_pvalue",
    "calibration_damage_qvalue",
    "calibration_variation_qvalue",
    "calibration_damage_instability",
)


@dataclass(slots=True, frozen=True)
class EventCalibrationScores:
    """The conservative score representation consumed by calibration."""

    observations: int
    damage_min: float
    damage_max: float
    error: float
    variation: float
    alternative_observations: int = 0
    reference_observations: int = 0


@dataclass(slots=True, frozen=True)
class CalibrationExample:
    """One labelled event from an independently simulated or validated sample."""

    truth: str
    scores: EventCalibrationScores
    sample_id: str


@dataclass(slots=True, frozen=True)
class ConformalCalibrationModel:
    """Class-conditional nonconformity distributions from held-out samples."""

    schema_version: int
    method: str
    normalization: str
    calibration_samples: tuple[str, ...]
    nonconformity_by_class: dict[str, tuple[float, ...]]


@dataclass(slots=True, frozen=True)
class EventConfidence:
    """Multiple-testing-aware, conservative event decision."""

    status: str
    reasons: tuple[str, ...]
    prediction_set: tuple[str, ...]
    decision: str
    damage_pvalue: float | None
    error_pvalue: float | None
    variation_pvalue: float | None
    damage_qvalue: float | None
    variation_qvalue: float | None
    damage_instability: float | None


@dataclass(slots=True, frozen=True)
class EventCalibrationSummary:
    rows: int
    scored: int
    eligible_error: int
    protect_damage: int
    protect_variation: int
    insufficient: int


def scores_from_likelihood(
    likelihood: EventLikelihood,
) -> EventCalibrationScores | None:
    """Use lower damage support for protection and upper support for removal."""
    values = (
        likelihood.damage_log_likelihood_min,
        likelihood.damage_log_likelihood_max,
        likelihood.error_log_likelihood,
        likelihood.variation_penalized_log_likelihood,
    )
    if likelihood.status != "scored" or any(value is None for value in values):
        return None
    assert all(value is not None for value in values)
    return EventCalibrationScores(
        observations=likelihood.observations,
        damage_min=values[0],
        damage_max=values[1],
        error=values[2],
        variation=values[3],
        alternative_observations=likelihood.alternative_observations,
        reference_observations=likelihood.reference_observations,
    )


def _nonconformity(scores: EventCalibrationScores, label: str) -> float:
    scale = sqrt(max(1, scores.observations))
    if label == "damage":
        return (max(scores.error, scores.variation) - scores.damage_min) / scale
    if label == "sequencing_error":
        return (max(scores.damage_max, scores.variation) - scores.error) / scale
    if label == "variation":
        return (max(scores.damage_max, scores.error) - scores.variation) / scale
    raise ValueError(f"unknown calibration class: {label}")


def fit_conformal_calibration(
    examples: Sequence[CalibrationExample],
    *,
    minimum_per_class: int = 100,
) -> ConformalCalibrationModel:
    """Fit class-conditional calibration without mixing events across samples."""
    if minimum_per_class < 1:
        raise ValueError("minimum_per_class must be positive")
    grouped = {label: [] for label in EXPLANATIONS}
    sample_ids = set()
    for example in examples:
        if example.truth not in grouped:
            raise ValueError(f"unknown calibration truth: {example.truth}")
        if not example.sample_id:
            raise ValueError("calibration examples require a sample_id")
        grouped[example.truth].append(_nonconformity(example.scores, example.truth))
        sample_ids.add(example.sample_id)
    insufficient = [
        label for label, values in grouped.items() if len(values) < minimum_per_class
    ]
    if insufficient:
        raise ValueError(
            "insufficient calibration examples for " + ", ".join(insufficient)
        )
    return ConformalCalibrationModel(
        schema_version=1,
        method="mondrian_class_conditional_conformal",
        normalization="winning_score_gap_divided_by_sqrt_observations",
        calibration_samples=tuple(sorted(sample_ids)),
        nonconformity_by_class={
            label: tuple(sorted(values)) for label, values in grouped.items()
        },
    )


def write_calibration_model(
    model: ConformalCalibrationModel,
    path: str | Path,
) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(asdict(model), indent=2) + "\n", encoding="utf-8")


def load_calibration_model(path: str | Path) -> ConformalCalibrationModel:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported event calibration schema")
    distributions = payload.get("nonconformity_by_class", {})
    if set(distributions) != set(EXPLANATIONS):
        raise ValueError("calibration model must contain all three classes")
    converted = {
        label: tuple(float(value) for value in distributions[label])
        for label in EXPLANATIONS
    }
    if any(not values for values in converted.values()):
        raise ValueError("calibration class distributions must not be empty")
    return ConformalCalibrationModel(
        schema_version=1,
        method=str(payload.get("method", "")),
        normalization=str(payload.get("normalization", "")),
        calibration_samples=tuple(payload.get("calibration_samples", ())),
        nonconformity_by_class=converted,
    )


def _conformal_pvalue(values: Sequence[float], observed: float) -> float:
    return (1 + len(values) - bisect_left(values, observed)) / (len(values) + 1)


def _benjamini_hochberg(values: Sequence[float]) -> list[float]:
    """Return monotone BH-adjusted values in original order."""
    count = len(values)
    if count == 0:
        return []
    order = sorted(range(count), key=values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank_index in range(count - 1, -1, -1):
        original = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[original] * count / rank)
        adjusted[original] = min(1.0, running)
    return adjusted


def calibrate_score_batch(
    scores: Sequence[EventCalibrationScores | None],
    model: ConformalCalibrationModel,
    *,
    alpha: float = 0.01,
    minimum_observations: int = 5,
    minimum_alternative_observations: int = 2,
    minimum_reference_observations: int = 5,
    maximum_damage_instability: float = 1.0,
) -> list[EventConfidence]:
    """Calibrate a report together so protective nulls receive BH correction."""
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be between zero and one")
    if minimum_observations < 1:
        raise ValueError("minimum_observations must be positive")
    if minimum_alternative_observations < 1:
        raise ValueError("minimum_alternative_observations must be positive")
    if minimum_reference_observations < 1:
        raise ValueError("minimum_reference_observations must be positive")
    unresolved = [
        label
        for label in EXPLANATIONS
        if 1.0 / (len(model.nonconformity_by_class[label]) + 1) > alpha
    ]
    if unresolved:
        raise ValueError(
            "calibration model lacks p-value resolution for alpha in "
            + ", ".join(unresolved)
        )
    raw: list[dict[str, float] | None] = []
    valid_indices = []
    damage_pvalues = []
    variation_pvalues = []
    for index, item in enumerate(scores):
        if item is None:
            raw.append(None)
            continue
        pvalues = {
            label: _conformal_pvalue(
                model.nonconformity_by_class[label],
                _nonconformity(item, label),
            )
            for label in EXPLANATIONS
        }
        raw.append(pvalues)
        valid_indices.append(index)
        damage_pvalues.append(pvalues["damage"])
        variation_pvalues.append(pvalues["variation"])

    damage_qvalues = _benjamini_hochberg(damage_pvalues)
    variation_qvalues = _benjamini_hochberg(variation_pvalues)
    q_by_index = {
        index: (damage_qvalues[offset], variation_qvalues[offset])
        for offset, index in enumerate(valid_indices)
    }
    results = []
    for index, item in enumerate(scores):
        pvalues = raw[index]
        if item is None or pvalues is None:
            results.append(
                EventConfidence(
                    "insufficient_evidence",
                    ("likelihood_not_scored",),
                    (),
                    "insufficient",
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                )
            )
            continue
        damage_qvalue, variation_qvalue = q_by_index[index]
        prediction_set = tuple(
            label for label in EXPLANATIONS if pvalues[label] > alpha
        )
        instability = (item.damage_max - item.damage_min) / sqrt(
            max(1, item.observations)
        )
        reasons = []
        decision = "insufficient"
        if item.observations < minimum_observations:
            reasons.append("too_few_observations")
        if item.alternative_observations < minimum_alternative_observations:
            reasons.append("too_few_alternative_observations")
        if item.reference_observations < minimum_reference_observations:
            reasons.append("too_few_reference_observations")
        if instability > maximum_damage_instability:
            reasons.append("damage_profile_unstable")
        if not reasons:
            if (
                pvalues["sequencing_error"] > alpha
                and damage_qvalue <= alpha
                and variation_qvalue <= alpha
            ):
                decision = "eligible_error"
            elif "damage" in prediction_set:
                decision = "protect_damage"
                if len(prediction_set) > 1:
                    reasons.append("multiple_plausible_explanations")
            elif "variation" in prediction_set:
                decision = "protect_variation"
                if len(prediction_set) > 1:
                    reasons.append("multiple_plausible_explanations")
            else:
                reasons.append("ambiguous_conformal_prediction_set")
        results.append(
            EventConfidence(
                status="calibrated" if decision != "insufficient" else "insufficient_evidence",
                reasons=tuple(reasons),
                prediction_set=prediction_set,
                decision=decision,
                damage_pvalue=pvalues["damage"],
                error_pvalue=pvalues["sequencing_error"],
                variation_pvalue=pvalues["variation"],
                damage_qvalue=damage_qvalue,
                variation_qvalue=variation_qvalue,
                damage_instability=instability,
            )
        )
    return results


def _parse_report_scores(row: dict[str, str]) -> EventCalibrationScores | None:
    if row.get("likelihood_status") != "scored":
        return None
    try:
        damage = float(row["damage_model_log_likelihood"])
        damage_min = float(row.get("damage_model_log_likelihood_min") or damage)
        damage_max = float(row.get("damage_model_log_likelihood_max") or damage)
        item = EventCalibrationScores(
            observations=int(row["likelihood_observations"]),
            damage_min=damage_min,
            damage_max=damage_max,
            error=float(row["error_model_log_likelihood"]),
            variation=float(row["variation_model_penalized_log_likelihood"]),
            alternative_observations=int(
                row["likelihood_alternative_observations"]
            ),
            reference_observations=int(row["likelihood_reference_observations"]),
        )
    except (KeyError, ValueError) as error:
        raise ValueError("event report has invalid likelihood fields") from error
    if not all(
        isfinite(value)
        for value in (item.damage_min, item.damage_max, item.error, item.variation)
    ):
        raise ValueError("event report contains non-finite likelihoods")
    return item


def calibrate_event_report(
    input_path: str | Path,
    model_path: str | Path,
    output_path: str | Path,
    *,
    alpha: float = 0.01,
    minimum_observations: int = 5,
    minimum_alternative_observations: int = 2,
    minimum_reference_observations: int = 5,
) -> EventCalibrationSummary:
    """Append calibrated confidence fields without rebuilding the graph."""
    source = Path(input_path)
    output = Path(output_path)
    if source.resolve() == output.resolve():
        raise ValueError("calibrated output must differ from the input report")
    model = load_calibration_model(model_path)
    with source.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("event report is missing a header")
        fieldnames = tuple(reader.fieldnames)
        scores = [_parse_report_scores(row) for row in reader]
    confidence = calibrate_score_batch(
        scores,
        model,
        alpha=alpha,
        minimum_observations=minimum_observations,
        minimum_alternative_observations=minimum_alternative_observations,
        minimum_reference_observations=minimum_reference_observations,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with source.open("r", encoding="utf-8", newline="") as source_handle, output.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(source_handle, delimiter="\t")
        writer = csv.DictWriter(
            handle,
            (*fieldnames, *CALIBRATION_FIELDS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row, result in zip(reader, confidence, strict=True):
            row.update(
                {
                    "calibration_status": result.status,
                    "calibration_reasons": ";".join(result.reasons),
                    "calibration_prediction_set": ";".join(result.prediction_set),
                    "calibration_decision": result.decision,
                    "calibration_damage_pvalue": _format_optional(result.damage_pvalue),
                    "calibration_error_pvalue": _format_optional(result.error_pvalue),
                    "calibration_variation_pvalue": _format_optional(result.variation_pvalue),
                    "calibration_damage_qvalue": _format_optional(result.damage_qvalue),
                    "calibration_variation_qvalue": _format_optional(result.variation_qvalue),
                    "calibration_damage_instability": _format_optional(result.damage_instability),
                }
            )
            writer.writerow(row)
    decisions = [result.decision for result in confidence]
    return EventCalibrationSummary(
        rows=len(scores),
        scored=sum(item is not None for item in scores),
        eligible_error=decisions.count("eligible_error"),
        protect_damage=decisions.count("protect_damage"),
        protect_variation=decisions.count("protect_variation"),
        insufficient=decisions.count("insufficient"),
    )


def _format_optional(value: float | None) -> str:
    return "" if value is None else f"{value:.8g}"
