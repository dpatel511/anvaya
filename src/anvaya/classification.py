"""Explainable, non-destructive classification of alternative unitig paths."""

from dataclasses import dataclass
from typing import Protocol


class PathEvidence(Protocol):
    """Evidence fields required by the path classifier."""

    left_terminal_observations: int
    right_terminal_observations: int
    relative_coverage: float
    strand_balance: float | None
    terminal_enrichment: float
    similarity_to_strongest: float


@dataclass(slots=True, frozen=True)
class ClassificationThresholds:
    """Conservative defaults derived from the controlled validation matrix."""

    maximum_error_coverage: float = 0.30
    minimum_error_similarity: float = 0.95
    maximum_error_terminal_enrichment: float = -0.15
    maximum_damage_coverage: float = 0.30
    maximum_damage_strand_balance: float = 0.50
    minimum_damage_terminal_enrichment: float = 0.0
    minimum_variation_coverage: float = 0.30
    minimum_variation_strand_balance: float = 0.50
    minimum_variation_terminal_enrichment: float = -0.15


@dataclass(slots=True, frozen=True)
class Substitution:
    """One sequence difference relative to the strongest local path."""

    position: int
    reference: str
    alternative: str


@dataclass(slots=True, frozen=True)
class PathClassification:
    """One label with auditable reasons and substitution evidence."""

    label: str
    reasons: tuple[str, ...]
    substitutions: tuple[Substitution, ...]
    damage_compatible: bool


def substitutions(
    reference: str, alternative: str
) -> tuple[Substitution, ...]:
    """Return substitutions for equal-length paths, or no calls otherwise."""
    if len(reference) != len(alternative):
        return ()
    return tuple(
        Substitution(index, reference_base, alternative_base)
        for index, (reference_base, alternative_base) in enumerate(
            zip(reference, alternative, strict=True),
            start=1,
        )
        if reference_base != alternative_base
    )


def _damage_compatible(
    changes: tuple[Substitution, ...], evidence: PathEvidence
) -> bool:
    if not changes:
        return False
    for change in changes:
        if (
            change.reference == "C"
            and change.alternative == "T"
            and evidence.left_terminal_observations > 0
        ):
            continue
        if (
            change.reference == "G"
            and change.alternative == "A"
            and evidence.right_terminal_observations > 0
        ):
            continue
        return False
    return True


def classify_unitig_path(
    reference_sequence: str,
    alternative_sequence: str,
    evidence: PathEvidence,
    *,
    dominant: bool = False,
    thresholds: ClassificationThresholds = ClassificationThresholds(),
) -> PathClassification:
    """Classify one path without changing its graph or support values."""
    changes = substitutions(reference_sequence, alternative_sequence)
    damage_compatible = _damage_compatible(changes, evidence)

    if dominant:
        return PathClassification(
            "dominant",
            ("highest_local_coverage",),
            changes,
            damage_compatible,
        )

    if (
        damage_compatible
        and evidence.relative_coverage
        <= thresholds.maximum_damage_coverage
        and evidence.strand_balance is not None
        and evidence.strand_balance
        <= thresholds.maximum_damage_strand_balance
        and evidence.terminal_enrichment
        >= thresholds.minimum_damage_terminal_enrichment
    ):
        return PathClassification(
            "damage-like",
            (
                "damage_compatible_substitutions",
                "low_local_coverage",
                "strand_skewed",
                "terminal_enriched",
            ),
            changes,
            True,
        )

    if (
        evidence.relative_coverage <= thresholds.maximum_error_coverage
        and evidence.similarity_to_strongest
        >= thresholds.minimum_error_similarity
        and evidence.terminal_enrichment
        <= thresholds.maximum_error_terminal_enrichment
    ):
        return PathClassification(
            "error-like",
            ("low_local_coverage", "high_sequence_similarity", "terminal_depleted"),
            changes,
            damage_compatible,
        )

    if (
        evidence.relative_coverage >= thresholds.minimum_variation_coverage
        and evidence.strand_balance is not None
        and evidence.strand_balance
        >= thresholds.minimum_variation_strand_balance
        and evidence.terminal_enrichment
        >= thresholds.minimum_variation_terminal_enrichment
    ):
        return PathClassification(
            "variation-like",
            ("supported_local_coverage", "strand_balanced", "not_terminal_depleted"),
            changes,
            damage_compatible,
        )

    return PathClassification(
        "ambiguous",
        ("conflicting_or_insufficient_evidence",),
        changes,
        damage_compatible,
    )


def format_substitutions(changes: tuple[Substitution, ...]) -> str:
    """Format substitutions for compact TSV output."""
    return ";".join(
        f"{change.position}:{change.reference}>{change.alternative}"
        for change in changes
    )
