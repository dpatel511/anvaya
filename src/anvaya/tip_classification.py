"""Explainable, non-destructive classification of matched weak tips."""

import math
from dataclasses import dataclass

from anvaya.bidirected import (
    BidirectedDeBruijnGraph,
    _successor,
)
from anvaya.tip_matching import TipBackboneMatch


@dataclass(slots=True, frozen=True)
class TipClassificationThresholds:
    """Conservative defaults derived from the controlled validation matrix."""

    minimum_damage_score: float = 0.80
    minimum_error_score: float = 0.75
    minimum_variation_score: float = 0.65
    minimum_score_margin: float = 0.15
    low_coverage_reference: float = 0.30
    maximum_error_terminal_enrichment: float = -0.15
    minimum_variation_coverage: float = 0.30


@dataclass(slots=True, frozen=True)
class TipSubstitutionEvidence:
    """End and internal support associated with one tip substitution."""

    position: int
    reference: str
    alternative: str
    expected_end: str | None
    expected_terminal_observations: int
    other_terminal_observations: int
    internal_observations: int
    mean_terminal_distance: float | None
    expected_molecule_count: int | None
    quality_observations: int
    missing_quality_observations: int
    mean_base_quality: float | None
    sequencing_error_log_likelihood: float | None


@dataclass(slots=True, frozen=True)
class TipEvidence:
    """Combined path and substitution evidence for one matched weak tip."""

    tip_terminal_observations: int
    tip_internal_observations: int
    tip_terminal_fraction: float
    tip_strand_balance: float | None
    backbone_terminal_fraction: float
    terminal_enrichment: float
    compatible_substitutions: int
    substitution_terminal_observations: int
    substitution_other_terminal_observations: int
    substitution_internal_observations: int
    substitution_terminal_fraction: float
    mean_damage_distance: float | None
    molecule_links_collected: bool
    joint_molecule_observations: int | None
    joint_molecule_fraction: float | None
    quality_observations: int
    missing_quality_observations: int
    mean_base_quality: float | None
    sequencing_error_log_likelihood: float | None


@dataclass(slots=True, frozen=True)
class TipClassification:
    """One conservative label with normalized evidence scores."""

    label: str
    reasons: tuple[str, ...]
    damage_score: float
    error_score: float
    variation_score: float
    evidence: TipEvidence
    substitutions: tuple[TipSubstitutionEvidence, ...]


@dataclass(slots=True, frozen=True)
class _OrientedEdgeEvidence:
    left: tuple[int, ...]
    right: tuple[int, ...]
    internal: int
    ambiguous_terminal: int
    forward: int
    reverse: int
    left_molecules: tuple[tuple[int, int, int | None], ...] | None
    right_molecules: tuple[tuple[int, int, int | None], ...] | None


def _clamp(value: float) -> float:
    return max(0.0, min(1.0, value))


def _sequencing_error_log_likelihood(qualities: list[int]) -> float | None:
    """Return the log likelihood of specific wrong calls under Phred error."""
    if not qualities:
        return None
    return sum(
        math.log((10.0 ** (-quality / 10.0)) / 3.0)
        for quality in qualities
    )


def _oriented_path_evidence(
    graph: BidirectedDeBruijnGraph,
    start: int,
    edge_ids: tuple[int, ...],
) -> tuple[_OrientedEdgeEvidence, ...]:
    current = start
    values: list[_OrientedEdgeEvidence] = []
    for edge_id in edge_ids:
        end_support = graph.end_support(edge_id)
        left = end_support.left
        right = end_support.right
        forward = graph.forward_support[edge_id]
        reverse = graph.reverse_support[edge_id]
        left_molecules: tuple[tuple[int, int, int | None], ...] | None = None
        right_molecules: tuple[tuple[int, int, int | None], ...] | None = None
        if graph.molecule_links_collected:
            links = graph.molecule_end_links(edge_id)
            left_molecules = tuple(
                (link.molecule_id, link.distance, link.base_quality)
                for link in links
                if link.end == "left"
            )
            right_molecules = tuple(
                (link.molecule_id, link.distance, link.base_quality)
                for link in links
                if link.end == "right"
            )
        if current != graph.edge_sources[edge_id]:
            left, right = right, left
            forward, reverse = reverse, forward
            left_molecules, right_molecules = (
                right_molecules,
                left_molecules,
            )
        values.append(
            _OrientedEdgeEvidence(
                left=left,
                right=right,
                internal=end_support.internal,
                ambiguous_terminal=end_support.ambiguous_terminal,
                forward=forward,
                reverse=reverse,
                left_molecules=left_molecules,
                right_molecules=right_molecules,
            )
        )
        current = _successor(graph, current, edge_id)
    return tuple(values)


def _path_summary(
    values: tuple[_OrientedEdgeEvidence, ...],
) -> tuple[int, int, float, float | None]:
    terminal = sum(
        sum(value.left)
        + sum(value.right)
        + value.ambiguous_terminal
        for value in values
    )
    internal = sum(value.internal for value in values)
    total = terminal + internal
    forward = sum(value.forward for value in values)
    reverse = sum(value.reverse for value in values)
    directional = forward + reverse
    return (
        terminal,
        internal,
        terminal / total if total else 0.0,
        2 * min(forward, reverse) / directional if directional else None,
    )


def _exact_end_counts(
    graph: BidirectedDeBruijnGraph,
    counts: tuple[int, ...],
    expected_end: str,
) -> tuple[int, ...]:
    """Move k-mer-boundary evidence to the appended nucleotide's cycle."""
    offset = graph.node_length if expected_end == "left" else 0
    exact = [0] * graph.end_window
    for distance, count in enumerate(counts):
        base_distance = distance + offset
        if base_distance < graph.end_window:
            exact[base_distance] += count
    return tuple(exact)


def collect_tip_evidence(
    graph: BidirectedDeBruijnGraph,
    match: TipBackboneMatch,
) -> tuple[TipEvidence, tuple[TipSubstitutionEvidence, ...]]:
    """Associate each substitution with evidence on its exact tip edge."""
    tip_edges = _oriented_path_evidence(
        graph,
        match.branch_handle,
        match.tip_edge_ids,
    )
    backbone_edges = _oriented_path_evidence(
        graph,
        match.branch_handle,
        match.backbone_edge_ids,
    )
    tip_terminal, tip_internal, tip_fraction, strand_balance = _path_summary(
        tip_edges
    )
    _, _, backbone_fraction, _ = _path_summary(backbone_edges)

    substitutions: list[TipSubstitutionEvidence] = []
    compatible = 0
    expected_total = 0
    other_total = 0
    internal_total = 0
    weighted_distance = 0.0
    expected_molecule_sets: list[set[int]] = []
    observed_qualities: list[int] = []
    missing_quality_total = 0

    for change in match.substitutions:
        edge_index = change.position - graph.node_length - 1
        expected_end: str | None = None
        expected_counts: tuple[int, ...] = ()
        other_terminal = 0
        internal = 0
        expected_molecules: set[int] | None = (
            set() if graph.molecule_links_collected else None
        )
        substitution_qualities: list[int] = []
        missing_qualities = 0
        if 0 <= edge_index < len(tip_edges):
            edge = tip_edges[edge_index]
            internal = edge.internal
            if change.reference == "C" and change.alternative == "T":
                expected_end = "left"
                expected_counts = _exact_end_counts(graph, edge.left, "left")
                if edge.left_molecules is not None:
                    expected_molecules = {
                        read_index
                        for read_index, distance, _ in edge.left_molecules
                        if distance + graph.node_length < graph.end_window
                    }
                    substitution_qualities = [
                        quality
                        for _, distance, quality in edge.left_molecules
                        if distance + graph.node_length < graph.end_window
                        and quality is not None
                    ]
                    missing_qualities = sum(
                        quality is None
                        for _, distance, quality in edge.left_molecules
                        if distance + graph.node_length < graph.end_window
                    )
                other_terminal = sum(edge.right) + edge.ambiguous_terminal
            elif change.reference == "G" and change.alternative == "A":
                expected_end = "right"
                expected_counts = _exact_end_counts(graph, edge.right, "right")
                if edge.right_molecules is not None:
                    expected_molecules = {
                        read_index
                        for read_index, distance, _ in edge.right_molecules
                        if distance < graph.end_window
                    }
                    substitution_qualities = [
                        quality
                        for _, distance, quality in edge.right_molecules
                        if distance < graph.end_window and quality is not None
                    ]
                    missing_qualities = sum(
                        quality is None
                        for _, distance, quality in edge.right_molecules
                        if distance < graph.end_window
                    )
                other_terminal = sum(edge.left) + edge.ambiguous_terminal
            else:
                other_terminal = (
                    sum(edge.left)
                    + sum(edge.right)
                    + edge.ambiguous_terminal
                )

        expected = sum(expected_counts)
        if expected_end is not None and expected:
            compatible += 1
        distance_sum = sum(
            distance * count
            for distance, count in enumerate(expected_counts)
        )
        expected_total += expected
        other_total += other_terminal
        internal_total += internal
        weighted_distance += distance_sum
        if expected_molecules is not None:
            expected_molecule_sets.append(expected_molecules)
        observed_qualities.extend(substitution_qualities)
        missing_quality_total += missing_qualities
        error_log_likelihood = _sequencing_error_log_likelihood(
            substitution_qualities
        )
        substitutions.append(
            TipSubstitutionEvidence(
                position=change.position,
                reference=change.reference,
                alternative=change.alternative,
                expected_end=expected_end,
                expected_terminal_observations=expected,
                other_terminal_observations=other_terminal,
                internal_observations=internal,
                mean_terminal_distance=(
                    distance_sum / expected if expected else None
                ),
                expected_molecule_count=(
                    None
                    if expected_molecules is None
                    else len(expected_molecules)
                ),
                quality_observations=len(substitution_qualities),
                missing_quality_observations=missing_qualities,
                mean_base_quality=(
                    sum(substitution_qualities) / len(substitution_qualities)
                    if substitution_qualities
                    else None
                ),
                sequencing_error_log_likelihood=error_log_likelihood,
            )
        )

    substitution_total = expected_total + other_total + internal_total
    joint_molecules: set[int] | None = None
    joint_fraction: float | None = None
    if graph.molecule_links_collected:
        joint_molecules = (
            set.intersection(*expected_molecule_sets)
            if expected_molecule_sets
            else set()
        )
        minimum_molecules = min(
            (len(values) for values in expected_molecule_sets),
            default=0,
        )
        joint_fraction = (
            len(joint_molecules) / minimum_molecules
            if minimum_molecules
            else 0.0
        )
    evidence = TipEvidence(
        tip_terminal_observations=tip_terminal,
        tip_internal_observations=tip_internal,
        tip_terminal_fraction=tip_fraction,
        tip_strand_balance=strand_balance,
        backbone_terminal_fraction=backbone_fraction,
        terminal_enrichment=tip_fraction - backbone_fraction,
        compatible_substitutions=compatible,
        substitution_terminal_observations=expected_total,
        substitution_other_terminal_observations=other_total,
        substitution_internal_observations=internal_total,
        substitution_terminal_fraction=(
            expected_total / substitution_total if substitution_total else 0.0
        ),
        mean_damage_distance=(
            weighted_distance / expected_total if expected_total else None
        ),
        molecule_links_collected=graph.molecule_links_collected,
        joint_molecule_observations=(
            None if joint_molecules is None else len(joint_molecules)
        ),
        joint_molecule_fraction=joint_fraction,
        quality_observations=len(observed_qualities),
        missing_quality_observations=missing_quality_total,
        mean_base_quality=(
            sum(observed_qualities) / len(observed_qualities)
            if observed_qualities
            else None
        ),
        sequencing_error_log_likelihood=_sequencing_error_log_likelihood(
            observed_qualities
        ),
    )
    return evidence, tuple(substitutions)


def classify_tip_match(
    graph: BidirectedDeBruijnGraph,
    match: TipBackboneMatch,
    thresholds: TipClassificationThresholds = TipClassificationThresholds(),
    *,
    include_molecule_linkage: bool = True,
) -> TipClassification:
    """Score damage, error, and variation explanations without graph changes."""
    if not isinstance(include_molecule_linkage, bool):
        raise TypeError("include_molecule_linkage must be a boolean")
    evidence, substitution_evidence = collect_tip_evidence(graph, match)
    change_count = len(match.substitutions)
    compatibility = (
        evidence.compatible_substitutions / change_count
        if change_count
        else 0.0
    )
    low_coverage = _clamp(
        1.0
        - match.relative_coverage / thresholds.low_coverage_reference
    )
    coverage_support = _clamp(
        match.relative_coverage / thresholds.low_coverage_reference
    )
    internal_fraction = (
        evidence.tip_internal_observations
        / (
            evidence.tip_terminal_observations
            + evidence.tip_internal_observations
        )
        if evidence.tip_terminal_observations
        + evidence.tip_internal_observations
        else 0.0
    )
    proximity = (
        1.0
        - evidence.mean_damage_distance / graph.end_window
        if evidence.mean_damage_distance is not None and graph.end_window
        else 0.0
    )
    strand_balance = evidence.tip_strand_balance or 0.0

    damage_components = [
        compatibility,
        match.ry_identity,
        evidence.substitution_terminal_fraction,
        _clamp(proximity),
        low_coverage,
    ]
    if (
        include_molecule_linkage
        and evidence.joint_molecule_fraction is not None
        and len(substitution_evidence) > 1
    ):
        damage_components.append(evidence.joint_molecule_fraction)
    damage_score = sum(damage_components) / len(damage_components)
    error_score = sum(
        (
            low_coverage,
            match.sequence_identity,
            1.0 - compatibility,
            1.0 - evidence.substitution_terminal_fraction,
        )
    ) / 4
    variation_score = sum(
        (
            coverage_support,
            internal_fraction,
            strand_balance,
            1.0 - low_coverage,
        )
    ) / 4

    scores = {
        "damage-like": damage_score,
        "error-like": error_score,
        "variation-like": variation_score,
    }
    label = "ambiguous"
    reasons: tuple[str, ...] = ("score_threshold_or_margin_not_met",)
    ranked = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    highest_label, highest_score = ranked[0]
    second_score = ranked[1][1]
    required = {
        "damage-like": thresholds.minimum_damage_score,
        "error-like": thresholds.minimum_error_score,
        "variation-like": thresholds.minimum_variation_score,
    }[highest_label]
    evidence_gate = True
    if highest_label == "error-like":
        evidence_gate = (
            evidence.terminal_enrichment
            <= thresholds.maximum_error_terminal_enrichment
        )
    elif highest_label == "variation-like":
        evidence_gate = (
            match.relative_coverage
            >= thresholds.minimum_variation_coverage
        )
    if not evidence_gate:
        reasons = (
            "terminal_depletion_not_observed"
            if highest_label == "error-like"
            else "supported_coverage_not_observed",
        )
    if (
        highest_score >= required
        and highest_score - second_score >= thresholds.minimum_score_margin
        and evidence_gate
    ):
        label = highest_label
        if label == "damage-like":
            reasons = (
                "damage_compatible_substitutions",
                "terminal_substitution_support",
                "end_proximal_support",
                "low_local_coverage",
                *(
                    ("molecule_linked_substitutions",)
                    if (
                        include_molecule_linkage
                        and len(substitution_evidence) > 1
                        and evidence.joint_molecule_fraction == 1.0
                    )
                    else ()
                ),
            )
        elif label == "error-like":
            reasons = (
                "low_local_coverage",
                "high_sequence_identity",
                "damage_evidence_depleted",
                "terminally_depleted",
            )
        else:
            reasons = (
                "supported_local_coverage",
                "internal_or_balanced_support",
            )

    return TipClassification(
        label=label,
        reasons=reasons,
        damage_score=damage_score,
        error_score=error_score,
        variation_score=variation_score,
        evidence=evidence,
        substitutions=substitution_evidence,
    )
