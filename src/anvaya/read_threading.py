"""Non-destructive read-threading audit for compacted graph junctions."""

import csv
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, replace
from pathlib import Path

from anvaya.bidirected import (
    BidirectedDeBruijnGraph,
    Handle,
    _positioned_encoded_kmers,
    _reverse_complement_code,
)
from anvaya.reads import Read
from anvaya.unitig_graph import CompactedUnitigGraph


@dataclass(slots=True, frozen=True)
class ThreadTransition:
    """Molecule-level evidence for one physical unitig transition."""

    source: Handle
    target: Handle
    molecules: int
    reads: int
    forward_molecules: int
    reverse_molecules: int
    terminal_molecules: int
    internal_molecules: int
    terminal_fraction: float
    quality_observations: int
    mean_quality: float | None
    source_mean_support: float
    target_mean_support: float
    target_support_ratio: float | None
    source_out_degree: int
    rank: int
    runner_up_molecules: int
    dominance_ratio: float | None
    resolvable: bool
    reciprocal: bool


@dataclass(slots=True, frozen=True)
class ReadThreadingSummary:
    """Counts describing junction coverage and conservative resolvability."""

    molecules: int = 0
    reads: int = 0
    junctions: int = 0
    candidate_transitions: int = 0
    supported_transitions: int = 0
    threaded_molecules: int = 0
    threaded_reads: int = 0
    resolvable_junctions: int = 0
    reciprocal_resolvable_links: int = 0


@dataclass(slots=True)
class _MutableEvidence:
    molecules: int = 0
    reads: int = 0
    forward_molecules: int = 0
    reverse_molecules: int = 0
    terminal_molecules: int = 0
    quality_sum: int = 0
    quality_observations: int = 0


def _canonical_transition(
    source: Handle,
    target: Handle,
) -> tuple[tuple[Handle, Handle], bool]:
    forward = (source, target)
    reverse = (target ^ 1, source ^ 1)
    return (forward, True) if forward <= reverse else (reverse, False)


def _edge_code(graph: BidirectedDeBruijnGraph, edge_id: int) -> int:
    """Reconstruct a canonical physical-edge code without decoding DNA."""
    source = graph.edge_sources[edge_id]
    target = graph.edge_targets[edge_id]
    source_code = graph.node_codes[source >> 1]
    target_code = graph.node_codes[target >> 1]
    if source & 1:
        source_code = _reverse_complement_code(source_code, graph.node_length)
    if target & 1:
        target_code = _reverse_complement_code(target_code, graph.node_length)
    return (source_code << 2) | (target_code & 3)


def _junction_edge_index(
    graph: BidirectedDeBruijnGraph,
    unitigs: CompactedUnitigGraph,
) -> dict[int, int]:
    """Index only physical edges immediately adjacent to branching links."""
    edge_ids: set[int] = set()
    for source in range(unitigs.unitig_count * 2):
        targets = unitigs.outgoing(source)
        if len(targets) < 2:
            continue
        unitig_id = source >> 1
        edge_ids.add(
            unitigs.last_edge_ids[unitig_id]
            if source & 1 == 0
            else unitigs.first_edge_ids[unitig_id]
        )
        for target in targets:
            target_id = target >> 1
            edge_ids.add(
                unitigs.first_edge_ids[target_id]
                if target & 1 == 0
                else unitigs.last_edge_ids[target_id]
            )
    return {_edge_code(graph, edge_id): edge_id for edge_id in edge_ids}


def audit_read_threads(
    graph: BidirectedDeBruijnGraph,
    unitigs: CompactedUnitigGraph,
    molecules: Iterable[Sequence[Read]],
    *,
    end_window: int = 0,
    minimum_support: int = 5,
    dominance_ratio: float = 3.0,
) -> tuple[tuple[ThreadTransition, ...], ReadThreadingSummary]:
    """Count direct read crossings of compacted-graph branch transitions."""
    if len(unitigs.edge_handles) != graph.edge_count:
        raise ValueError("read threading requires edge-to-unitig mappings")
    if end_window < 0:
        raise ValueError("end_window must not be negative")
    if minimum_support < 1:
        raise ValueError("minimum_support must be at least 1")
    if dominance_ratio < 1.0:
        raise ValueError("dominance_ratio must be at least 1")

    edge_index = _junction_edge_index(graph, unitigs)
    evidence: dict[tuple[Handle, Handle], _MutableEvidence] = {}
    molecule_count = 0
    read_count = 0
    threaded_reads = 0
    threaded_molecules = 0
    k = graph.node_length + 1

    for molecule in molecules:
        molecule_count += 1
        molecule_hits: dict[tuple[Handle, Handle], list[int | bool | None]] = {}
        molecule_threaded = False
        for read in molecule:
            read_count += 1
            read_hits: set[tuple[Handle, Handle]] = set()
            previous: tuple[int, int, int] | None = None
            last_start = len(read.sequence) - k
            for start, code, orientation in _positioned_encoded_kmers(
                read.sequence,
                k,
            ):
                edge_id = edge_index.get(code)
                if edge_id is None:
                    previous = None
                    continue
                handle = unitigs.edge_handles[edge_id]
                if handle < 0:
                    previous = None
                    continue
                handle ^= orientation & 1
                if previous is not None and start == previous[0] + 1:
                    source = previous[1]
                    if (source >> 1) != (handle >> 1) and handle in unitigs.outgoing(source):
                        key, is_forward = _canonical_transition(source, handle)
                        terminal = (
                            previous[0] < end_window
                            or last_start - start < end_window
                        )
                        quality = (
                            None
                            if read.qualities is None
                            else read.qualities[start + k - 1]
                        )
                        hit = molecule_hits.setdefault(
                            key,
                            [False, False, terminal, quality],
                        )
                        hit[0 if is_forward else 1] = True
                        hit[2] = bool(hit[2]) or terminal
                        if quality is not None and (hit[3] is None or quality > hit[3]):
                            hit[3] = quality
                        read_hits.add(key)
                previous = (start, handle, orientation)
            if read_hits:
                threaded_reads += 1
                molecule_threaded = True
                for key in read_hits:
                    evidence.setdefault(key, _MutableEvidence()).reads += 1
        if molecule_threaded:
            threaded_molecules += 1
        for key, (forward, reverse, terminal, quality) in molecule_hits.items():
            item = evidence.setdefault(key, _MutableEvidence())
            item.molecules += 1
            if forward:
                item.forward_molecules += 1
            if reverse:
                item.reverse_molecules += 1
            if terminal:
                item.terminal_molecules += 1
            if quality is not None:
                item.quality_sum += int(quality)
                item.quality_observations += 1

    junctions = 0
    candidates = 0
    resolvable = 0
    rows: list[ThreadTransition] = []
    chosen: dict[Handle, Handle] = {}
    for source in range(unitigs.unitig_count * 2):
        targets = unitigs.outgoing(source)
        if len(targets) == 1:
            chosen[source] = targets[0]
            continue
        if len(targets) < 2:
            continue
        junctions += 1
        candidates += len(targets)
        scored = sorted(
            (
                (
                    evidence.get(
                        _canonical_transition(source, target)[0],
                        _MutableEvidence(),
                    ).molecules,
                    target,
                )
                for target in targets
            ),
            key=lambda item: (-item[0], item[1]),
        )
        winner_support = scored[0][0]
        runner_up_support = scored[1][0]
        is_resolvable = winner_support >= minimum_support and (
            runner_up_support == 0
            or winner_support >= dominance_ratio * runner_up_support
        )
        if is_resolvable:
            resolvable += 1
            chosen[source] = scored[0][1]
        for rank, (molecule_support, target) in enumerate(scored, start=1):
            key, forward = _canonical_transition(source, target)
            item = evidence.get(key, _MutableEvidence())
            source_id = source >> 1
            target_id = target >> 1
            source_support = (
                unitigs.observations[source_id] / unitigs.edge_counts[source_id]
            )
            target_support = (
                unitigs.observations[target_id] / unitigs.edge_counts[target_id]
            )
            rows.append(
                ThreadTransition(
                    source=source,
                    target=target,
                    molecules=molecule_support,
                    reads=item.reads,
                    forward_molecules=(
                        item.forward_molecules
                        if forward
                        else item.reverse_molecules
                    ),
                    reverse_molecules=(
                        item.reverse_molecules
                        if forward
                        else item.forward_molecules
                    ),
                    terminal_molecules=item.terminal_molecules,
                    internal_molecules=(
                        item.molecules - item.terminal_molecules
                    ),
                    terminal_fraction=(
                        item.terminal_molecules / item.molecules
                        if item.molecules
                        else 0.0
                    ),
                    quality_observations=item.quality_observations,
                    mean_quality=(
                        item.quality_sum / item.quality_observations
                        if item.quality_observations
                        else None
                    ),
                    source_mean_support=source_support,
                    target_mean_support=target_support,
                    target_support_ratio=(
                        None
                        if source_support == 0
                        else target_support / source_support
                    ),
                    source_out_degree=len(targets),
                    rank=rank,
                    runner_up_molecules=runner_up_support,
                    dominance_ratio=(
                        None
                        if runner_up_support == 0
                        else winner_support / runner_up_support
                    ),
                    resolvable=is_resolvable and rank == 1,
                    reciprocal=False,
                )
            )

    reciprocal_links = {
        _canonical_transition(source, target)[0]
        for source, target in chosen.items()
        if len(unitigs.outgoing(source)) >= 2
        and chosen.get(target ^ 1) == (source ^ 1)
    }
    rows = [
        replace(
            row,
            reciprocal=(
                row.resolvable
                and _canonical_transition(row.source, row.target)[0]
                in reciprocal_links
            ),
        )
        for row in rows
    ]

    return tuple(rows), ReadThreadingSummary(
        molecules=molecule_count,
        reads=read_count,
        junctions=junctions,
        candidate_transitions=candidates,
        supported_transitions=len(evidence),
        threaded_molecules=threaded_molecules,
        threaded_reads=threaded_reads,
        resolvable_junctions=resolvable,
        reciprocal_resolvable_links=len(reciprocal_links),
    )


def write_read_thread_report(
    transitions: Sequence[ThreadTransition],
    path: str | Path,
) -> None:
    """Write deterministic molecule-level transition evidence as TSV."""
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "source_unitig",
                "source_orientation",
                "target_unitig",
                "target_orientation",
                "molecules",
                "reads",
                "forward_molecules",
                "reverse_molecules",
                "terminal_molecules",
                "internal_molecules",
                "terminal_fraction",
                "quality_observations",
                "mean_quality",
                "source_mean_support",
                "target_mean_support",
                "target_support_ratio",
                "source_out_degree",
                "rank",
                "runner_up_molecules",
                "dominance_ratio",
                "resolvable",
                "reciprocal",
            )
        )
        for item in transitions:
            writer.writerow(
                (
                    item.source // 2 + 1,
                    "reverse" if item.source & 1 else "forward",
                    item.target // 2 + 1,
                    "reverse" if item.target & 1 else "forward",
                    item.molecules,
                    item.reads,
                    item.forward_molecules,
                    item.reverse_molecules,
                    item.terminal_molecules,
                    item.internal_molecules,
                    f"{item.terminal_fraction:.6f}",
                    item.quality_observations,
                    "" if item.mean_quality is None else f"{item.mean_quality:.3f}",
                    f"{item.source_mean_support:.3f}",
                    f"{item.target_mean_support:.3f}",
                    (
                        ""
                        if item.target_support_ratio is None
                        else f"{item.target_support_ratio:.3f}"
                    ),
                    item.source_out_degree,
                    item.rank,
                    item.runner_up_molecules,
                    (
                        ""
                        if item.dominance_ratio is None
                        else f"{item.dominance_ratio:.3f}"
                    ),
                    str(item.resolvable).lower(),
                    str(item.reciprocal).lower(),
                )
            )


def resolve_read_thread_extensions(
    transitions: Sequence[ThreadTransition],
    *,
    repeat_coverage_ratio: float = 2.0,
) -> dict[Handle, Handle]:
    """Return reciprocal winners, breaking chains through coverage repeats."""
    if repeat_coverage_ratio < 1.0:
        raise ValueError("repeat_coverage_ratio must be at least 1")

    reciprocal = {
        _canonical_transition(item.source, item.target)[0]: item
        for item in transitions
        if item.reciprocal
    }
    accepted: dict[Handle, Handle] = {}
    coverage: dict[Handle, float] = {}
    for item in reciprocal.values():
        accepted[item.source] = item.target
        accepted[item.target ^ 1] = item.source ^ 1
        coverage[item.source] = item.source_mean_support
        coverage[item.source ^ 1] = item.source_mean_support
        coverage[item.target] = item.target_mean_support
        coverage[item.target ^ 1] = item.target_mean_support

    predecessor = {target: source for source, target in accepted.items()}
    rejected: set[tuple[Handle, Handle]] = set()
    for middle, target in accepted.items():
        source = predecessor.get(middle)
        if source is None:
            continue
        first = _canonical_transition(source, middle)[0]
        second = _canonical_transition(middle, target)[0]
        if first == second:
            continue
        middle_coverage = coverage[middle]
        flank_coverage = max(coverage[source], coverage[target])
        if middle_coverage < repeat_coverage_ratio * flank_coverage:
            continue
        first_support = reciprocal[first].molecules
        second_support = reciprocal[second].molecules
        rejected.add(
            first
            if (first_support, first) < (second_support, second)
            else second
        )

    successors: dict[Handle, Handle] = {}
    for transition in transitions:
        physical = _canonical_transition(transition.source, transition.target)[0]
        if not transition.reciprocal or physical in rejected:
            continue
        source = transition.source
        target = transition.target
        reverse_source = target ^ 1
        reverse_target = source ^ 1
        if source in successors and successors[source] != target:
            raise ValueError("conflicting threaded successor assignments")
        if (
            reverse_source in successors
            and successors[reverse_source] != reverse_target
        ):
            raise ValueError("conflicting reverse threaded assignments")
        successors[source] = target
        successors[reverse_source] = reverse_target
    return successors
