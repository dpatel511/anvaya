"""Non-destructive bubble scoring on the compacted unitig graph."""

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path

from anvaya.unitig_graph import CompactedUnitigGraph


@dataclass(slots=True, frozen=True)
class UnitigBubblePath:
    """Support and similarity scores for one alternative unitig path."""

    handles: tuple[int, ...]
    nucleotide_length: int
    edge_count: int
    observations: int
    minimum_edge_support: int
    mean_edge_support: float
    relative_coverage: float
    similarity_to_strongest: float


@dataclass(slots=True, frozen=True)
class UnitigBubble:
    """Bounded alternatives that leave and rejoin the same unitigs."""

    start: int
    end: int
    strongest_path_index: int
    paths: tuple[UnitigBubblePath, ...]


@dataclass(slots=True, frozen=True)
class UnitigBubbleReportSummary:
    """Counts written to a unitig-bubble report."""

    bubbles: int = 0
    paths: int = 0


@dataclass(slots=True, frozen=True)
class _RawPath:
    handles: tuple[int, ...]
    sequence: str
    edge_count: int
    observations: int
    minimum_edge_support: int
    mean_edge_support: float


def _linear_branch(
    graph: CompactedUnitigGraph,
    start: int,
    first: int,
    max_unitigs: int,
    max_nucleotides: int,
) -> tuple[int, tuple[int, ...]] | None:
    """Follow one branch to its first non-linear oriented unitig."""
    path: list[int] = []
    current = first
    visited = {start}
    nucleotide_length = 0

    while current not in visited:
        incoming_degree = len(graph.outgoing(current ^ 1))
        outgoing = graph.outgoing(current)
        if incoming_degree != 1 or len(outgoing) != 1:
            return current, tuple(path)

        if len(path) == max_unitigs:
            return None
        sequence_length = len(graph.oriented_sequence(current))
        nucleotide_length += (
            sequence_length if not path else sequence_length - graph.k + 1
        )
        if nucleotide_length > max_nucleotides:
            return None

        visited.add(current)
        path.append(current)
        current = outgoing[0]

    return None


def _spell_path(graph: CompactedUnitigGraph, handles: Sequence[int]) -> str:
    sequence = graph.oriented_sequence(handles[0])
    overlap = graph.k - 1
    return sequence + "".join(
        graph.oriented_sequence(handle)[overlap:] for handle in handles[1:]
    )


def _score_path(
    graph: CompactedUnitigGraph, handles: tuple[int, ...]
) -> _RawPath:
    supports = [graph.support(handle >> 1) for handle in handles]
    edge_count = sum(support.edge_count for support in supports)
    observations = sum(support.observations for support in supports)
    return _RawPath(
        handles=handles,
        sequence=_spell_path(graph, handles),
        edge_count=edge_count,
        observations=observations,
        minimum_edge_support=min(
            support.minimum_edge_support for support in supports
        ),
        mean_edge_support=observations / edge_count,
    )


def find_unitig_bubbles(
    graph: CompactedUnitigGraph,
    max_unitigs: int = 32,
    max_nucleotides: int = 2_000,
) -> list[UnitigBubble]:
    """Score bounded, unitig-disjoint alternatives without changing the graph."""
    for value, name in (
        (max_unitigs, "max_unitigs"),
        (max_nucleotides, "max_nucleotides"),
    ):
        if not isinstance(value, int) or isinstance(value, bool):
            raise TypeError(f"{name} must be an integer")
        if value < 1:
            raise ValueError(f"{name} must be at least 1")

    bubbles: list[UnitigBubble] = []
    seen: set[tuple[tuple[int, ...], ...]] = set()

    for start in range(graph.unitig_count * 2):
        outgoing = graph.outgoing(start)
        if len(outgoing) < 2:
            continue

        branches = [
            _linear_branch(
                graph,
                start,
                first,
                max_unitigs,
                max_nucleotides,
            )
            for first in sorted(outgoing)
        ]
        if any(branch is None for branch in branches):
            continue

        completed = [branch for branch in branches if branch is not None]
        ends = {end for end, _ in completed}
        paths = [path for _, path in completed]
        if len(ends) != 1 or any(not path for path in paths):
            continue
        end = completed[0][0]
        if end == start or len(graph.outgoing(end ^ 1)) < 2:
            continue

        used_unitigs: set[int] = set()
        disjoint = True
        for path in paths:
            physical_unitigs = {handle >> 1 for handle in path}
            if len(physical_unitigs) != len(path) or used_unitigs & physical_unitigs:
                disjoint = False
                break
            used_unitigs.update(physical_unitigs)
        if not disjoint:
            continue

        signature = tuple(
            sorted(tuple(sorted(handle >> 1 for handle in path)) for path in paths)
        )
        if signature in seen:
            continue
        seen.add(signature)

        raw_paths = [_score_path(graph, path) for path in paths]
        strongest_index = max(
            range(len(raw_paths)),
            key=lambda index: (
                raw_paths[index].mean_edge_support,
                raw_paths[index].minimum_edge_support,
                raw_paths[index].observations,
                -index,
            ),
        )
        strongest = raw_paths[strongest_index]
        scored_paths = tuple(
            UnitigBubblePath(
                handles=path.handles,
                nucleotide_length=len(path.sequence),
                edge_count=path.edge_count,
                observations=path.observations,
                minimum_edge_support=path.minimum_edge_support,
                mean_edge_support=path.mean_edge_support,
                relative_coverage=(
                    path.mean_edge_support / strongest.mean_edge_support
                ),
                similarity_to_strongest=(
                    1.0
                    if index == strongest_index
                    else SequenceMatcher(
                        None,
                        strongest.sequence,
                        path.sequence,
                        autojunk=False,
                    ).ratio()
                ),
            )
            for index, path in enumerate(raw_paths)
        )
        bubbles.append(
            UnitigBubble(
                start=start,
                end=end,
                strongest_path_index=strongest_index,
                paths=scored_paths,
            )
        )

    return bubbles


def write_unitig_bubble_report(
    bubbles: Sequence[UnitigBubble], path: str | Path
) -> UnitigBubbleReportSummary:
    """Write raw path scores for threshold validation as TSV."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "bubble_id",
                "path_index",
                "strongest_path",
                "start_handle",
                "end_handle",
                "unitig_handles",
                "nucleotide_length",
                "edge_count",
                "observations",
                "minimum_edge_support",
                "mean_edge_support",
                "relative_coverage",
                "similarity_to_strongest",
            ]
        )
        path_count = 0
        for bubble_index, bubble in enumerate(bubbles, start=1):
            for path_index, bubble_path in enumerate(bubble.paths):
                writer.writerow(
                    [
                        f"unitig-bubble-{bubble_index}",
                        path_index,
                        str(path_index == bubble.strongest_path_index).lower(),
                        bubble.start,
                        bubble.end,
                        ",".join(map(str, bubble_path.handles)),
                        bubble_path.nucleotide_length,
                        bubble_path.edge_count,
                        bubble_path.observations,
                        bubble_path.minimum_edge_support,
                        f"{bubble_path.mean_edge_support:.6f}",
                        f"{bubble_path.relative_coverage:.6f}",
                        f"{bubble_path.similarity_to_strongest:.6f}",
                    ]
                )
                path_count += 1

    return UnitigBubbleReportSummary(
        bubbles=len(bubbles),
        paths=path_count,
    )
