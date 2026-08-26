"""Sequence and read-end evidence reports for graph alternatives."""

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from anvaya.bidirected import (
    BidirectedDeBruijnGraph,
    _edge_support_total,
    _successor,
    spell_edge_path,
)
from anvaya.bubbles import Bubble, BubblePath
from anvaya.cleaning import TipCandidate
from anvaya.classification import format_substitutions
from anvaya.tip_matching import match_tip_to_backbone
from anvaya.tip_classification import classify_tip_match


@dataclass(slots=True, frozen=True)
class EventReportSummary:
    """Counts written to one graph-event report."""

    tips: int = 0
    matched_tips: int = 0
    bubbles: int = 0
    paths: int = 0


_TIP_MATCH_FIELDS = (
    "tip_match_sequence",
    "backbone_matched",
    "backbone_sequence",
    "backbone_edge_count",
    "backbone_observations",
    "backbone_minimum_edge_support",
    "relative_coverage",
    "sequence_identity",
    "ry_identity",
)

_TIP_CLASSIFICATION_FIELDS = (
    "tip_classification",
    "tip_classification_reasons",
    "tip_damage_score",
    "tip_error_score",
    "tip_variation_score",
    "tip_terminal_fraction",
    "tip_strand_balance",
    "tip_terminal_enrichment",
    "substitution_terminal_observations",
    "substitution_other_terminal_observations",
    "substitution_internal_observations",
    "substitution_terminal_fraction",
    "mean_damage_distance",
)


def _oriented_edge_evidence(
    graph: BidirectedDeBruijnGraph,
    start: int,
    edge_ids: Sequence[int],
) -> list[tuple[int, int, int, int]]:
    values: list[tuple[int, int, int, int]] = []
    current = start

    for edge_id in edge_ids:
        evidence = graph.end_support(edge_id)
        edge_left = sum(evidence.left)
        edge_right = sum(evidence.right)
        if current != graph.edge_sources[edge_id]:
            edge_left, edge_right = edge_right, edge_left
        values.append(
            (
                edge_left,
                edge_right,
                evidence.internal,
                evidence.ambiguous_terminal,
            )
        )
        current = _successor(graph, current, edge_id)

    return values


def _path_evidence(
    graph: BidirectedDeBruijnGraph,
    start: int,
    edge_ids: Sequence[int],
) -> tuple[int, int, int, int]:
    values = _oriented_edge_evidence(graph, start, edge_ids)
    return tuple(sum(value[index] for value in values) for index in range(4))


def _substitutions(reference: str, alternative: str) -> list[tuple[int, str, str]]:
    if len(reference) != len(alternative):
        return []
    return [
        (position, reference_base, alternative_base)
        for position, (reference_base, alternative_base) in enumerate(
            zip(reference, alternative, strict=True),
            start=1,
        )
        if reference_base != alternative_base
    ]


def _format_substitutions(
    substitutions: Sequence[tuple[int, str, str]],
) -> str:
    return ";".join(
        f"{position}:{reference}>{alternative}"
        for position, reference, alternative in substitutions
    )


def _damage_compatibility(
    graph: BidirectedDeBruijnGraph,
    start: int,
    edge_ids: Sequence[int],
    substitutions: Sequence[tuple[int, str, str]],
) -> str:
    if not substitutions:
        return "unknown"
    edge_evidence = _oriented_edge_evidence(graph, start, edge_ids)

    def is_compatible(
        position: int,
        reference: str,
        alternative: str,
    ) -> bool:
        edge_index = position - graph.node_length - 1
        if not 0 <= edge_index < len(edge_evidence):
            return False
        left_terminal, right_terminal, _, _ = edge_evidence[edge_index]
        return (
            reference == "C"
            and alternative == "T"
            and left_terminal > 0
        ) or (
            reference == "G"
            and alternative == "A"
            and right_terminal > 0
        )

    compatible = all(
        is_compatible(position, reference, alternative)
        for position, reference, alternative in substitutions
    )
    return str(compatible).lower()


def _bubble_rows(
    graph: BidirectedDeBruijnGraph,
    event_id: str,
    bubble: Bubble,
) -> list[list[object]]:
    sequences = [
        spell_edge_path(graph, bubble.start, path.edge_ids)
        for path in bubble.paths
    ]
    reference_index = max(
        range(len(bubble.paths)),
        key=lambda index: (
            bubble.paths[index].minimum_edge_support,
            bubble.paths[index].observations,
            -index,
        ),
    )
    reference_sequence = sequences[reference_index]
    rows: list[list[object]] = []

    for index, (path, sequence) in enumerate(
        zip(bubble.paths, sequences, strict=True)
    ):
        left, right, internal, ambiguous = _path_evidence(
            graph,
            bubble.start,
            path.edge_ids,
        )
        substitutions = (
            []
            if index == reference_index
            else _substitutions(reference_sequence, sequence)
        )
        compatibility = (
            "reference"
            if index == reference_index
            else _damage_compatibility(
                graph,
                bubble.start,
                path.edge_ids,
                substitutions,
            )
        )
        rows.append(
            [
                event_id,
                "bubble",
                index,
                str(index == reference_index).lower(),
                bubble.start,
                bubble.end,
                len(path.edge_ids),
                path.observations,
                path.minimum_edge_support,
                left,
                right,
                internal,
                ambiguous,
                sequence,
                _format_substitutions(substitutions),
                compatibility,
            ]
            + [""]
            * (len(_TIP_MATCH_FIELDS) + len(_TIP_CLASSIFICATION_FIELDS))
        )
    return rows


def _tip_row(
    graph: BidirectedDeBruijnGraph,
    event_id: str,
    tip: TipCandidate,
) -> tuple[list[object], bool]:
    left, right, internal, ambiguous = _path_evidence(
        graph,
        tip.start,
        tip.edge_ids,
    )
    observations = sum(
        _edge_support_total(graph, edge_id) for edge_id in tip.edge_ids
    )
    minimum_support = min(
        _edge_support_total(graph, edge_id) for edge_id in tip.edge_ids
    )
    match = match_tip_to_backbone(graph, tip)
    decision = None if match is None else classify_tip_match(graph, match)
    substitutions = "" if match is None else format_substitutions(match.substitutions)
    compatibility = (
        "unknown" if match is None else str(match.damage_compatible).lower()
    )
    row = [
        event_id,
        "tip",
        0,
        "false",
        tip.start,
        tip.end,
        len(tip.edge_ids),
        observations,
        minimum_support,
        left,
        right,
        internal,
        ambiguous,
        spell_edge_path(graph, tip.start, tip.edge_ids),
        substitutions,
        compatibility,
        "" if match is None else match.tip_sequence,
        str(match is not None).lower(),
        "" if match is None else match.backbone_sequence,
        "" if match is None else len(match.backbone_edge_ids),
        "" if match is None else match.backbone_observations,
        "" if match is None else match.backbone_minimum_edge_support,
        "" if match is None else f"{match.relative_coverage:.6f}",
        "" if match is None else f"{match.sequence_identity:.6f}",
        "" if match is None else f"{match.ry_identity:.6f}",
        "" if decision is None else decision.label,
        "" if decision is None else ";".join(decision.reasons),
        "" if decision is None else f"{decision.damage_score:.6f}",
        "" if decision is None else f"{decision.error_score:.6f}",
        "" if decision is None else f"{decision.variation_score:.6f}",
        (
            ""
            if decision is None
            else f"{decision.evidence.tip_terminal_fraction:.6f}"
        ),
        (
            ""
            if decision is None
            or decision.evidence.tip_strand_balance is None
            else f"{decision.evidence.tip_strand_balance:.6f}"
        ),
        (
            ""
            if decision is None
            else f"{decision.evidence.terminal_enrichment:.6f}"
        ),
        (
            ""
            if decision is None
            else decision.evidence.substitution_terminal_observations
        ),
        (
            ""
            if decision is None
            else decision.evidence.substitution_other_terminal_observations
        ),
        (
            ""
            if decision is None
            else decision.evidence.substitution_internal_observations
        ),
        (
            ""
            if decision is None
            else f"{decision.evidence.substitution_terminal_fraction:.6f}"
        ),
        (
            ""
            if decision is None
            or decision.evidence.mean_damage_distance is None
            else f"{decision.evidence.mean_damage_distance:.6f}"
        ),
    ]
    return row, match is not None


def write_event_report(
    graph: BidirectedDeBruijnGraph,
    tips: Sequence[TipCandidate],
    bubbles: Sequence[Bubble],
    path: str | Path,
) -> EventReportSummary:
    """Write non-destructive tip and bubble-path evidence as TSV."""
    if graph.end_window == 0:
        raise ValueError("event reporting requires read-end evidence")

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "event_id",
                "event_type",
                "path_index",
                "reference_path",
                "start_handle",
                "end_handle",
                "edge_count",
                "observations",
                "minimum_edge_support",
                "left_terminal",
                "right_terminal",
                "internal",
                "ambiguous_terminal",
                "sequence",
                "substitutions",
                "damage_compatible",
                *_TIP_MATCH_FIELDS,
                *_TIP_CLASSIFICATION_FIELDS,
            ]
        )
        path_count = 0
        matched_tip_count = 0
        for index, tip in enumerate(tips, start=1):
            row, matched = _tip_row(graph, f"tip-{index}", tip)
            writer.writerow(row)
            matched_tip_count += matched
            path_count += 1
        for index, bubble in enumerate(bubbles, start=1):
            rows = _bubble_rows(graph, f"bubble-{index}", bubble)
            writer.writerows(rows)
            path_count += len(rows)

    return EventReportSummary(
        tips=len(tips),
        matched_tips=matched_tip_count,
        bubbles=len(bubbles),
        paths=path_count,
    )
