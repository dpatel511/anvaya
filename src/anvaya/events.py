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
from anvaya.incomplete_branches import (
    IncompleteBranchCandidate,
    match_incomplete_branch_to_backbone,
)
from anvaya.tip_classification import classify_tip_match
from anvaya.tip_matching import TipBackboneMatch, match_tip_to_backbone


@dataclass(slots=True, frozen=True)
class EventReportSummary:
    """Counts written to one graph-event report."""

    tips: int = 0
    matched_tips: int = 0
    incomplete_branches: int = 0
    matched_incomplete_branches: int = 0
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


def _match_columns(
    graph: BidirectedDeBruijnGraph,
    match: TipBackboneMatch | None,
) -> tuple[str, str, list[object]]:
    if match is None:
        return (
            "",
            "unknown",
            [""]
            + ["false"]
            + [""] * (len(_TIP_MATCH_FIELDS) - 2)
            + [""] * len(_TIP_CLASSIFICATION_FIELDS),
        )

    decision = classify_tip_match(graph, match)
    columns: list[object] = [
        match.tip_sequence,
        "true",
        match.backbone_sequence,
        len(match.backbone_edge_ids),
        match.backbone_observations,
        match.backbone_minimum_edge_support,
        f"{match.relative_coverage:.6f}",
        f"{match.sequence_identity:.6f}",
        f"{match.ry_identity:.6f}",
        decision.label,
        ";".join(decision.reasons),
        f"{decision.damage_score:.6f}",
        f"{decision.error_score:.6f}",
        f"{decision.variation_score:.6f}",
        f"{decision.evidence.tip_terminal_fraction:.6f}",
        (
            ""
            if decision.evidence.tip_strand_balance is None
            else f"{decision.evidence.tip_strand_balance:.6f}"
        ),
        f"{decision.evidence.terminal_enrichment:.6f}",
        decision.evidence.substitution_terminal_observations,
        decision.evidence.substitution_other_terminal_observations,
        decision.evidence.substitution_internal_observations,
        f"{decision.evidence.substitution_terminal_fraction:.6f}",
        (
            ""
            if decision.evidence.mean_damage_distance is None
            else f"{decision.evidence.mean_damage_distance:.6f}"
        ),
    ]
    return (
        format_substitutions(match.substitutions),
        str(match.damage_compatible).lower(),
        columns,
    )


def _matched_path_row(
    graph: BidirectedDeBruijnGraph,
    event_id: str,
    event_type: str,
    start: int,
    end: int,
    edge_ids: tuple[int, ...],
    observations: int,
    minimum_support: int,
    sequence: str,
    match: TipBackboneMatch | None,
) -> list[object]:
    left, right, internal, ambiguous = _path_evidence(
        graph,
        start,
        edge_ids,
    )
    substitutions, compatibility, match_columns = _match_columns(graph, match)
    return [
        event_id,
        event_type,
        0,
        "false",
        start,
        end,
        len(edge_ids),
        observations,
        minimum_support,
        left,
        right,
        internal,
        ambiguous,
        sequence,
        substitutions,
        compatibility,
        *match_columns,
    ]


def _tip_row(
    graph: BidirectedDeBruijnGraph,
    event_id: str,
    tip: TipCandidate,
) -> tuple[list[object], bool]:
    observations = sum(
        _edge_support_total(graph, edge_id) for edge_id in tip.edge_ids
    )
    minimum_support = min(
        _edge_support_total(graph, edge_id) for edge_id in tip.edge_ids
    )
    match = match_tip_to_backbone(graph, tip)
    return (
        _matched_path_row(
            graph,
            event_id,
            "tip",
            tip.start,
            tip.end,
            tip.edge_ids,
            observations,
            minimum_support,
            spell_edge_path(graph, tip.start, tip.edge_ids),
            match,
        ),
        match is not None,
    )


def _incomplete_branch_row(
    graph: BidirectedDeBruijnGraph,
    event_id: str,
    candidate: IncompleteBranchCandidate,
) -> tuple[list[object], bool]:
    match = match_incomplete_branch_to_backbone(graph, candidate)
    return (
        _matched_path_row(
            graph,
            event_id,
            "incomplete_branch",
            candidate.start,
            candidate.end,
            candidate.edge_ids,
            candidate.observations,
            candidate.minimum_edge_support,
            spell_edge_path(graph, candidate.start, candidate.edge_ids),
            match,
        ),
        match is not None,
    )


def write_event_report(
    graph: BidirectedDeBruijnGraph,
    tips: Sequence[TipCandidate],
    bubbles: Sequence[Bubble],
    path: str | Path,
    incomplete_branches: Sequence[IncompleteBranchCandidate] = (),
) -> EventReportSummary:
    """Write non-destructive tip, branch, and bubble-path evidence as TSV."""
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
        matched_incomplete_count = 0
        for index, candidate in enumerate(incomplete_branches, start=1):
            row, matched = _incomplete_branch_row(
                graph,
                f"incomplete-branch-{index}",
                candidate,
            )
            writer.writerow(row)
            matched_incomplete_count += matched
            path_count += 1
        for index, bubble in enumerate(bubbles, start=1):
            rows = _bubble_rows(graph, f"bubble-{index}", bubble)
            writer.writerows(rows)
            path_count += len(rows)

    return EventReportSummary(
        tips=len(tips),
        matched_tips=matched_tip_count,
        incomplete_branches=len(incomplete_branches),
        matched_incomplete_branches=matched_incomplete_count,
        bubbles=len(bubbles),
        paths=path_count,
    )
