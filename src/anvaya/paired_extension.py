"""Conservative paired-read extension of compacted unitig paths."""

from collections import deque
from dataclasses import dataclass

from anvaya.bidirected import BidirectedDeBruijnGraph, Handle
from anvaya.sequences import canonical_sequence
from anvaya.unitig_graph import CompactedUnitigGraph


@dataclass(slots=True, frozen=True)
class PairedExtensionSummary:
    """Counts describing paired-read path extension."""

    paired_molecules: int = 0
    mapped_pairs: int = 0
    evidence_links: int = 0
    resolved_oriented_junctions: int = 0
    evaluated_oriented_junctions: int = 0
    work_limited_oriented_junctions: int = 0
    searched_states: int = 0
    joined_unitig_links: int = 0
    output_contigs: int = 0


@dataclass(slots=True, frozen=True)
class PairedExtensionResult:
    """Extended sequences and their audit summary."""

    sequences: tuple[str, ...]
    summary: PairedExtensionSummary


def collect_paired_unitig_links(
    source_graph: BidirectedDeBruijnGraph,
    unitig_graph: CompactedUnitigGraph,
) -> tuple[dict[tuple[Handle, Handle], int], int, int]:
    """Map paired read ends to directional oriented-unitig links."""
    if not source_graph.molecule_links_collected:
        raise ValueError("paired extension requires molecule links")
    if len(unitig_graph.edge_handles) != source_graph.edge_count:
        raise ValueError("paired extension requires edge-to-unitig mappings")

    # One nearest retained k-mer anchors each read's three-prime end. A tie
    # between different unitigs is ambiguous and therefore discarded.
    endpoints: dict[int, tuple[int, Handle] | None] = {}
    for edge_id, base_handle in enumerate(unitig_graph.edge_handles):
        if base_handle < 0:
            continue
        for link in source_graph.molecule_end_links(edge_id):
            normalized_side = 0 if link.end == "left" else 1
            actual_side = normalized_side ^ link.orientation
            if actual_side != 1:
                continue
            handle = base_handle ^ link.orientation
            previous = endpoints.get(link.read_index)
            current = (link.distance, handle)
            if previous is None and link.read_index in endpoints:
                continue
            if previous is None or current[0] < previous[0]:
                endpoints[link.read_index] = current
            elif current[0] == previous[0] and current[1] != previous[1]:
                endpoints[link.read_index] = None

    by_molecule: dict[int, list[tuple[int, Handle]]] = {}
    for read_index, endpoint in endpoints.items():
        if endpoint is None:
            continue
        molecule_id = source_graph.read_molecule_ids[read_index]
        by_molecule.setdefault(molecule_id, []).append(
            (read_index, endpoint[1])
        )

    molecule_reads: dict[int, int] = {}
    for molecule_id in source_graph.read_molecule_ids:
        molecule_reads[molecule_id] = molecule_reads.get(molecule_id, 0) + 1

    evidence: dict[tuple[Handle, Handle], int] = {}
    paired_molecules = sum(count == 2 for count in molecule_reads.values())
    mapped_pairs = 0
    for observations in by_molecule.values():
        read_indices = {read_index for read_index, _ in observations}
        if len(observations) != 2 or len(read_indices) != 2:
            continue
        observations.sort()
        first = observations[0][1]
        second = observations[1][1]
        evidence[(first, second ^ 1)] = (
            evidence.get((first, second ^ 1), 0) + 1
        )
        evidence[(second, first ^ 1)] = (
            evidence.get((second, first ^ 1), 0) + 1
        )
        mapped_pairs += 1
    return evidence, paired_molecules, mapped_pairs


def _candidate_ownership(
    graph: CompactedUnitigGraph,
    candidates: tuple[Handle, ...],
    targets: set[Handle],
    maximum_distance: int,
    maximum_states: int,
) -> tuple[dict[Handle, Handle | None], int, bool]:
    """Assign observed targets to candidates within a bounded work budget."""
    owners: dict[Handle, Handle | None] = {}
    reached_by: dict[Handle, Handle] = {}
    searched_states = 0
    for candidate in candidates:
        queue = deque([(candidate, 0)])
        best_distance: dict[Handle, int] = {}
        remaining = set(targets)
        while queue:
            handle, distance = queue.popleft()
            previous_distance = best_distance.get(handle)
            if previous_distance is not None and previous_distance <= distance:
                continue
            best_distance[handle] = distance
            searched_states += 1
            if searched_states > maximum_states:
                return owners, searched_states, True

            if handle in remaining:
                remaining.remove(handle)
                previous_owner = reached_by.get(handle)
                if previous_owner is None:
                    reached_by[handle] = candidate
                    owners[handle] = candidate
                elif previous_owner != candidate:
                    owners[handle] = None
                if not remaining:
                    break

            added = max(
                1,
                len(graph.oriented_sequence(handle)) - graph.k + 1,
            )
            next_distance = distance + added
            if next_distance > maximum_distance:
                continue
            for successor in graph.outgoing(handle):
                queue.append((successor, next_distance))
    return owners, searched_states, False


@dataclass(slots=True, frozen=True)
class _ResolutionAudit:
    evaluated: int = 0
    work_limited: int = 0
    searched_states: int = 0


def _resolve_paired_extensions(
    graph: CompactedUnitigGraph,
    evidence: dict[tuple[Handle, Handle], int],
    *,
    minimum_support: int,
    dominance_ratio: float,
    maximum_distance: int,
    maximum_states: int,
) -> tuple[dict[Handle, Handle], int, _ResolutionAudit]:
    """Resolve links and retain search-work diagnostics."""
    evidence_by_source: dict[Handle, list[tuple[Handle, int]]] = {}
    for (source, target), count in evidence.items():
        evidence_by_source.setdefault(source, []).append((target, count))

    chosen: dict[Handle, Handle] = {}
    resolved = 0
    evaluated = 0
    work_limited = 0
    searched_states = 0
    for source in range(graph.unitig_count * 2):
        outgoing = graph.outgoing(source)
        if len(outgoing) == 1:
            chosen[source] = outgoing[0]

    for source, observations in evidence_by_source.items():
        outgoing = graph.outgoing(source)
        if len(outgoing) == 1:
            continue
        if len(outgoing) < 2 or sum(count for _, count in observations) < minimum_support:
            continue

        evaluated += 1
        scores = {candidate: 0 for candidate in outgoing}
        owners, states, limited = _candidate_ownership(
            graph,
            outgoing,
            {target for target, _ in observations},
            maximum_distance,
            maximum_states,
        )
        searched_states += states
        if limited:
            work_limited += 1
            continue
        for observed_target, count in observations:
            candidate = owners.get(observed_target)
            if candidate is not None:
                scores[candidate] += count

        ranked = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
        winner, winner_score = ranked[0]
        runner_up_score = ranked[1][1]
        if winner_score < minimum_support:
            continue
        if runner_up_score and winner_score < dominance_ratio * runner_up_score:
            continue
        chosen[source] = winner
        resolved += 1

    accepted = {
        source: target
        for source, target in chosen.items()
        if (source >> 1) != (target >> 1)
        and chosen.get(target ^ 1) == (source ^ 1)
    }
    return accepted, resolved, _ResolutionAudit(
        evaluated=evaluated,
        work_limited=work_limited,
        searched_states=searched_states,
    )


def resolve_paired_extensions(
    graph: CompactedUnitigGraph,
    evidence: dict[tuple[Handle, Handle], int],
    *,
    minimum_support: int = 5,
    dominance_ratio: float = 3.0,
    maximum_distance: int = 1000,
    maximum_states: int = 1_000,
) -> tuple[dict[Handle, Handle], int]:
    """Choose reciprocal unitig links using an exSPAnder-style strong winner."""
    if minimum_support < 1:
        raise ValueError("minimum_support must be at least 1")
    if dominance_ratio < 1.0:
        raise ValueError("dominance_ratio must be at least 1")
    if maximum_distance < 1:
        raise ValueError("maximum_distance must be at least 1")
    if maximum_states < 1:
        raise ValueError("maximum_states must be at least 1")
    accepted, resolved, _ = _resolve_paired_extensions(
        graph,
        evidence,
        minimum_support=minimum_support,
        dominance_ratio=dominance_ratio,
        maximum_distance=maximum_distance,
        maximum_states=maximum_states,
    )
    return accepted, resolved


def spell_extended_paths(
    graph: CompactedUnitigGraph,
    successors: dict[Handle, Handle],
) -> tuple[str, ...]:
    """Spell a physical-unitig path cover from reciprocal successor links."""
    predecessors = {target: source for source, target in successors.items()}
    seen_unitigs: set[int] = set()
    sequences: list[str] = []

    def append_path(start: Handle) -> None:
        path: list[Handle] = []
        current = start
        while current >> 1 not in seen_unitigs:
            path.append(current)
            seen_unitigs.add(current >> 1)
            successor = successors.get(current)
            if successor is None:
                break
            current = successor
        if not path:
            return
        sequence = graph.oriented_sequence(path[0])
        for handle in path[1:]:
            sequence += graph.oriented_sequence(handle)[graph.k - 1 :]
        sequences.append(canonical_sequence(sequence))

    for handle in range(graph.unitig_count * 2):
        if handle not in predecessors:
            append_path(handle)
    for handle in range(graph.unitig_count * 2):
        append_path(handle)
    return tuple(sequences)


def extend_paired_unitig_paths(
    source_graph: BidirectedDeBruijnGraph,
    unitig_graph: CompactedUnitigGraph,
    *,
    minimum_support: int = 5,
    dominance_ratio: float = 3.0,
    maximum_distance: int = 1000,
    maximum_states: int = 1_000,
) -> PairedExtensionResult:
    """Extend compacted unitigs only across reciprocal strong-winner links."""
    evidence, paired_molecules, mapped_pairs = collect_paired_unitig_links(
        source_graph,
        unitig_graph,
    )
    successors, resolved, audit = _resolve_paired_extensions(
        unitig_graph,
        evidence,
        minimum_support=minimum_support,
        dominance_ratio=dominance_ratio,
        maximum_distance=maximum_distance,
        maximum_states=maximum_states,
    )
    sequences = spell_extended_paths(unitig_graph, successors)
    return PairedExtensionResult(
        sequences=sequences,
        summary=PairedExtensionSummary(
            paired_molecules=paired_molecules,
            mapped_pairs=mapped_pairs,
            evidence_links=len(evidence),
            resolved_oriented_junctions=resolved,
            evaluated_oriented_junctions=audit.evaluated,
            work_limited_oriented_junctions=audit.work_limited,
            searched_states=audit.searched_states,
            joined_unitig_links=len(successors) // 2,
            output_contigs=len(sequences),
        ),
    )
