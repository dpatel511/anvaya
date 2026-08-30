"""Report-only gate for projecting matched damage paths onto a backbone."""

from dataclasses import dataclass

from anvaya.event_likelihood import EventLikelihood
from anvaya.tip_classification import TipClassification
from anvaya.tip_matching import TipBackboneMatch


@dataclass(slots=True, frozen=True)
class ProjectionDecision:
    """One non-destructive decision about a matched path."""

    decision: str
    reasons: tuple[str, ...]
    target_sequence: str


def decide_damage_projection(
    match: TipBackboneMatch,
    classification: TipClassification,
    likelihood: EventLikelihood,
    *,
    minimum_log_likelihood_margin: float = 0.5,
) -> ProjectionDecision:
    """Return whether a damage-like path would be projected to its backbone.

    This is intentionally report-only. A projection requires independent
    likelihood support for damage in addition to the existing conservative
    topology and terminal-evidence classification.
    """
    reasons: list[str] = []
    if classification.label != "damage-like":
        reasons.append("not_damage_like")
    if likelihood.status != "scored":
        reasons.append("likelihood_not_scored")
    elif likelihood.best_explanation != "damage":
        reasons.append("damage_not_best_explanation")
    elif (
        likelihood.log_likelihood_margin is None
        or likelihood.log_likelihood_margin < minimum_log_likelihood_margin
    ):
        reasons.append("damage_likelihood_margin_too_small")
    if (
        likelihood.damage_variation_log_contrast is None
        or likelihood.damage_variation_log_contrast <= 0.0
    ):
        reasons.append("damage_does_not_beat_variation")
    if likelihood.missing_quality_observations:
        reasons.append("missing_base_quality")

    approved = not reasons
    return ProjectionDecision(
        decision="would_project" if approved else "leave_unresolved",
        reasons=tuple(reasons) if reasons else ("all_projection_gates_passed",),
        target_sequence=match.backbone_sequence if approved else "",
    )
