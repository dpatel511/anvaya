"""Research-only exact string graph for variable-length reads.

All offsets are 0-based. Topology is built from unique canonical strings while
raw-read placements remain immutable evidence for later consensus integration.
This module is intentionally unavailable from the production CLI: its Phase B
damage integration failed the frozen coordinate-safety gate.
"""

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass

from anvaya.damage_string_graph import _damage_compatible_edge
from anvaya.overlap_assembly import _CandidateDiagnostics, _candidate_alignments
from anvaya.overlap_index import _anchor_index
from anvaya.overlap_progressive import (
    ProgressiveSequencePool,
    SequenceRecord,
    SequenceState,
)
from anvaya.raw_consensus import DamageProfile, RawPlacement
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


OrientedNode = tuple[int, bool]


@dataclass(frozen=True, slots=True)
class ExactEdge:
    source: OrientedNode
    target: OrientedNode
    shift: int
    overlap: int
    corrections: tuple[tuple[int, str], ...] = ()


@dataclass(frozen=True, slots=True)
class EvidencePlacement:
    raw_index: int
    topology_node: int
    reverse: bool
    offset: int


@dataclass(frozen=True, slots=True)
class ExactUnitig:
    sequence: str
    layout: tuple[tuple[int, bool, int], ...]


@dataclass(frozen=True, slots=True)
class ExactStringGraph:
    topology_sequences: tuple[str, ...]
    candidate_edges: tuple[ExactEdge, ...]
    reduced_edges: tuple[ExactEdge, ...]
    unitigs: tuple[ExactUnitig, ...]
    evidence_placements: tuple[EvidencePlacement, ...]
    contained_groups: int
    ambiguous_evidence_reads: int


@dataclass(frozen=True, slots=True)
class UnitigEvidencePool:
    """Immutable raw records plus independently stored derived unitigs."""

    records: tuple[SequenceRecord, ...]
    derived: tuple[SequenceRecord, ...]

    @property
    def active_derived(self) -> tuple[SequenceRecord, ...]:
        return self.derived


def _canonical_groups(reads: list[Read]):
    groups: dict[str, list[tuple[int, bool]]] = {}
    for index, read in enumerate(reads):
        reverse = reverse_complement(read.sequence)
        canonical = min(read.sequence, reverse)
        groups.setdefault(canonical, []).append((index, read.sequence != canonical))
    return tuple(groups), tuple(tuple(groups[sequence]) for sequence in groups)


def _containments(sequences: tuple[str, ...], maximum_seed_length: int = 31):
    """Return all strict exact child placements in canonical parent strings."""
    seeds: dict[int, dict[str, list[tuple[int, bool, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for child_index, child in enumerate(sequences):
        oriented = [(False, child)]
        reverse = reverse_complement(child)
        if reverse != child:
            oriented.append((True, reverse))
        seed_length = min(maximum_seed_length, len(child))
        for is_reverse, sequence in oriented:
            seeds[seed_length][sequence[:seed_length]].append(
                (child_index, is_reverse, sequence)
            )

    placements: dict[int, set[tuple[int, bool, int]]] = defaultdict(set)
    for parent_index, parent in enumerate(sequences):
        for seed_length, index in seeds.items():
            if seed_length >= len(parent):
                continue
            for offset in range(len(parent) - seed_length + 1):
                for child_index, reverse, child in index.get(
                    parent[offset:offset + seed_length], ()
                ):
                    if len(child) < len(parent) and parent.startswith(child, offset):
                        placements[child_index].add((parent_index, reverse, offset))
    return {child: placements[child] for child in sorted(placements)}


def _oriented_sequence(sequences: tuple[str, ...], node: OrientedNode) -> str:
    sequence = sequences[node[0]]
    return reverse_complement(sequence) if node[1] else sequence


def _mirror(edge: ExactEdge, sequences: tuple[str, ...]) -> ExactEdge:
    target_length = len(sequences[edge.target[0]])
    return ExactEdge(
        (edge.target[0], not edge.target[1]),
        (edge.source[0], not edge.source[1]),
        target_length - edge.overlap,
        edge.overlap,
        tuple(sorted(
            (
                target_length - 1 - (position - edge.shift),
                reverse_complement(base),
            )
            for position, base in edge.corrections
        )),
    )


def _candidate_edges(
    sequences: tuple[str, ...], minimum_overlap: int
) -> dict[tuple[OrientedNode, OrientedNode], ExactEdge]:
    """Find all maximal proper exact dovetails using one suffix seed per source."""
    occurrences: dict[str, list[tuple[OrientedNode, int]]] = defaultdict(list)
    oriented = [
        ((index, reverse), _oriented_sequence(sequences, (index, reverse)))
        for index in range(len(sequences)) for reverse in (False, True)
    ]
    for node, sequence in oriented:
        for position in range(len(sequence) - minimum_overlap + 1):
            occurrences[sequence[position:position + minimum_overlap]].append(
                (node, position)
            )

    edges: dict[tuple[OrientedNode, OrientedNode], ExactEdge] = {}
    for source, source_sequence in oriented:
        if len(source_sequence) <= minimum_overlap:
            continue
        seed = source_sequence[-minimum_overlap:]
        for target, position in occurrences.get(seed, ()):
            if source[0] == target[0]:
                continue
            target_sequence = _oriented_sequence(sequences, target)
            overlap = position + minimum_overlap
            if overlap >= len(source_sequence) or overlap >= len(target_sequence):
                continue
            if source_sequence[-overlap:] != target_sequence[:overlap]:
                continue
            edge = ExactEdge(
                source, target, len(source_sequence) - overlap, overlap
            )
            key = (source, target)
            if key not in edges or overlap > edges[key].overlap:
                edges[key] = edge

    for edge in tuple(edges.values()):
        mirror = _mirror(edge, sequences)
        current = edges.get((mirror.source, mirror.target))
        if current is None or mirror.overlap > current.overlap:
            edges[(mirror.source, mirror.target)] = mirror
    return edges


def _transitive_reduction(
    edges: dict[tuple[OrientedNode, OrientedNode], ExactEdge],
    sequences: tuple[str, ...],
) -> dict[tuple[OrientedNode, OrientedNode], ExactEdge]:
    outgoing: dict[OrientedNode, list[ExactEdge]] = defaultdict(list)
    for edge in edges.values():
        outgoing[edge.source].append(edge)
    for values in outgoing.values():
        values.sort(key=lambda edge: (edge.shift, edge.target))

    def spell(path: tuple[ExactEdge, ...]) -> str | None:
        sequence = list(_oriented_sequence(sequences, path[0].source))
        current = path[0].source
        source_start = 0
        for edge in path:
            if edge.source != current:
                return None
            if source_start + edge.shift != len(sequence) - edge.overlap:
                return None
            for position, base in edge.corrections:
                assembled_position = source_start + position
                if not 0 <= assembled_position < len(sequence):
                    return None
                sequence[assembled_position] = base
            sequence.extend(
                _oriented_sequence(sequences, edge.target)[edge.overlap:]
            )
            source_start += edge.shift
            current = edge.target
        return "".join(sequence)

    def has_alternative(direct: ExactEdge) -> bool:
        expected = spell((direct,))
        stack = [
            (edge.target, edge.shift, (edge,), frozenset((direct.source, edge.target)))
            for edge in outgoing[direct.source]
            if edge != direct and edge.shift <= direct.shift
        ]
        while stack:
            node, shift, path, visited = stack.pop()
            if (
                node == direct.target
                and shift == direct.shift
                and spell(path) == expected
            ):
                return True
            for edge in outgoing.get(node, ()):
                total = shift + edge.shift
                if total <= direct.shift and edge.target not in visited:
                    stack.append((
                        edge.target,
                        total,
                        (*path, edge),
                        visited | {edge.target},
                    ))
        return False

    removed: set[tuple[OrientedNode, OrientedNode]] = set()
    for edge in sorted(edges.values(), key=lambda item: (item.shift, item.source, item.target)):
        key = (edge.source, edge.target)
        if key not in removed and has_alternative(edge):
            removed.add(key)
            mirror = _mirror(edge, sequences)
            removed.add((mirror.source, mirror.target))
    return {key: edge for key, edge in edges.items() if key not in removed}


def _linear_key(path: tuple[OrientedNode, ...]):
    mirror = tuple((index, not reverse) for index, reverse in reversed(path))
    return min(path, mirror)


def _cycle_key(path: tuple[OrientedNode, ...]):
    mirror = tuple((index, not reverse) for index, reverse in reversed(path))
    rotations = [path[i:] + path[:i] for i in range(len(path))]
    rotations.extend(mirror[i:] + mirror[:i] for i in range(len(mirror)))
    return min(rotations)


def _canonical_cycle(
    path: tuple[OrientedNode, ...], sequences: tuple[str, ...]
) -> tuple[OrientedNode, ...]:
    """Choose a cycle break from sequence content, independent of input IDs."""
    mirror = tuple((index, not reverse) for index, reverse in reversed(path))
    rotations = [path[i:] + path[:i] for i in range(len(path))]
    rotations.extend(mirror[i:] + mirror[:i] for i in range(len(mirror)))
    return min(
        rotations,
        key=lambda candidate: tuple(
            _oriented_sequence(sequences, node) for node in candidate
        ),
    )


def _spell(
    path: tuple[OrientedNode, ...],
    edges: dict[tuple[OrientedNode, OrientedNode], ExactEdge],
    sequences: tuple[str, ...],
) -> ExactUnitig:
    assembled = _oriented_sequence(sequences, path[0])
    layout = [(path[0][0], path[0][1], 0)]
    offset = 0
    for source, target in zip(path, path[1:]):
        edge = edges[(source, target)]
        for position, base in edge.corrections:
            assembled_position = offset + position
            if not 0 <= assembled_position < len(assembled):
                raise ValueError("unitig correction is outside its source")
            assembled = (
                assembled[:assembled_position]
                + base
                + assembled[assembled_position + 1:]
            )
        offset += edge.shift
        target_sequence = _oriented_sequence(sequences, target)
        if not edge.corrections and (
            assembled[offset:offset + edge.overlap] != target_sequence[:edge.overlap]
        ):
            raise ValueError("exact unitig edge does not spell its declared overlap")
        assembled += target_sequence[edge.overlap:]
        layout.append((target[0], target[1], offset))

    reverse = reverse_complement(assembled)
    if reverse < assembled:
        length = len(assembled)
        layout = [
            (index, not is_reverse, length - offset - len(sequences[index]))
            for index, is_reverse, offset in reversed(layout)
        ]
        assembled = reverse
    return ExactUnitig(assembled, tuple(layout))


def _unitigs(
    sequences: tuple[str, ...],
    edges: dict[tuple[OrientedNode, OrientedNode], ExactEdge],
) -> tuple[ExactUnitig, ...]:
    outgoing: dict[OrientedNode, list[ExactEdge]] = defaultdict(list)
    incoming: dict[OrientedNode, list[ExactEdge]] = defaultdict(list)
    for edge in edges.values():
        outgoing[edge.source].append(edge)
        incoming[edge.target].append(edge)
    for values in (*outgoing.values(), *incoming.values()):
        values.sort(key=lambda edge: (edge.target, edge.shift))

    paths: dict[tuple[OrientedNode, ...], tuple[OrientedNode, ...]] = {}
    visited_edges: set[tuple[OrientedNode, OrientedNode]] = set()
    nodes = sorted(set(outgoing) | set(incoming))
    for start in nodes:
        if len(incoming[start]) == 1 and len(outgoing[start]) == 1:
            continue
        for first in outgoing[start]:
            path = [start, first.target]
            visited_edges.add((first.source, first.target))
            while len(incoming[path[-1]]) == 1 and len(outgoing[path[-1]]) == 1:
                edge = outgoing[path[-1]][0]
                if edge.target in path:
                    break
                visited_edges.add((edge.source, edge.target))
                path.append(edge.target)
            candidate = tuple(path)
            paths.setdefault(_linear_key(candidate), candidate)

    for first in sorted(edges.values(), key=lambda edge: (edge.source, edge.target)):
        if (first.source, first.target) in visited_edges:
            continue
        path = [first.source]
        node = first.source
        while True:
            edge = outgoing[node][0]
            visited_edges.add((edge.source, edge.target))
            if edge.target == path[0]:
                break
            if edge.target in path:
                break
            path.append(edge.target)
            node = edge.target
        candidate = tuple(path)
        paths.setdefault(
            _cycle_key(candidate), _canonical_cycle(candidate, sequences)
        )

    connected = {node[0] for node in nodes}
    unitigs = [_spell(path, edges, sequences) for path in paths.values()]
    unitigs.extend(
        ExactUnitig(sequence, ((index, False, 0),))
        for index, sequence in enumerate(sequences) if index not in connected
    )
    return tuple(sorted(unitigs, key=lambda unitig: (unitig.sequence, unitig.layout)))


def assemble_exact_unitigs(
    reads: list[Read], *, minimum_overlap: int = 30
) -> ExactStringGraph:
    """Build an exact containment-free graph without changing raw evidence."""
    if not reads:
        raise ValueError("exact string graph requires at least one read")
    if minimum_overlap < 1:
        raise ValueError("minimum overlap must be at least 1")

    groups, raw_groups = _canonical_groups(reads)
    containments = _containments(groups)
    maximal_groups = [index for index in range(len(groups)) if index not in containments]
    topology = tuple(groups[index] for index in maximal_groups)
    topology_by_group = {group: node for node, group in enumerate(maximal_groups)}

    evidence: list[EvidencePlacement] = []
    for group_index, raw_members in enumerate(raw_groups):
        if group_index in topology_by_group:
            placements = ((topology_by_group[group_index], False, 0),)
        else:
            placements = tuple(
                (topology_by_group[parent], reverse, offset)
                for parent, reverse, offset in sorted(containments[group_index])
                if parent in topology_by_group
            )
            if not placements:
                raise ValueError("contained sequence has no maximal placement")
        for raw_index, raw_reverse in raw_members:
            evidence.extend(
                EvidencePlacement(
                    raw_index, topology_node, raw_reverse != child_reverse, offset
                )
                for topology_node, child_reverse, offset in placements
            )

    candidate = _candidate_edges(topology, minimum_overlap)
    reduced = _transitive_reduction(candidate, topology)
    placements_per_read: dict[int, int] = defaultdict(int)
    for placement in evidence:
        placements_per_read[placement.raw_index] += 1
    return ExactStringGraph(
        topology,
        tuple(sorted(candidate.values(), key=lambda edge: (edge.source, edge.target))),
        tuple(sorted(reduced.values(), key=lambda edge: (edge.source, edge.target))),
        _unitigs(topology, reduced),
        tuple(sorted(evidence, key=lambda item: (
            item.raw_index, item.topology_node, item.reverse, item.offset
        ))),
        len(containments),
        sum(count > 1 for count in placements_per_read.values()),
    )


def _add_damage_edge(
    edges: dict[tuple[OrientedNode, OrientedNode], ExactEdge],
    edge: ExactEdge,
    sequences: tuple[str, ...],
    ambiguous: set[tuple[OrientedNode, OrientedNode]],
) -> bool:
    """Insert one physical bidirected edge, abstaining on conflicting data."""
    physical = min(
        (edge.source, edge.target),
        ((edge.target[0], not edge.target[1]),
         (edge.source[0], not edge.source[1])),
    )
    if physical in ambiguous:
        return False
    mirror = _mirror(edge, sequences)
    current = edges.get((edge.source, edge.target))
    if current is not None:
        if current == edge:
            return True
        ambiguous.add(physical)
        edges.pop((edge.source, edge.target), None)
        edges.pop((mirror.source, mirror.target), None)
        return False
    edges[(edge.source, edge.target)] = edge
    edges[(mirror.source, mirror.target)] = mirror
    return True


def assemble_damage_aware_unitigs(
    reads: list[Read],
    profile: DamageProfile,
    molecule_ids: list[int] | None = None,
    *,
    minimum_overlap: int = 30,
) -> tuple[ExactStringGraph, dict]:
    """Add one evidence-gated damage-compatible layer to the exact topology."""
    molecules = list(range(len(reads))) if molecule_ids is None else molecule_ids
    if len(molecules) != len(reads):
        raise ValueError("molecule IDs must align one-to-one with reads")
    exact = assemble_exact_unitigs(reads, minimum_overlap=minimum_overlap)
    topology_reads = [
        Read(f"topology_{index}", sequence)
        for index, sequence in enumerate(exact.topology_sequences)
    ]
    maximum_length = max(len(read.sequence) for read in reads)
    position_bits = maximum_length.bit_length()
    topology_anchors = _anchor_index(topology_reads, 15, 0, 8, 100)
    raw_anchors = _anchor_index(reads, 15, 0, 8, 100)
    topology_ids = list(range(len(topology_reads)))
    raw_ids = list(range(len(reads)))
    candidate_diagnostics = _CandidateDiagnostics()
    raw_diagnostics = _CandidateDiagnostics()
    classifications = Counter()
    edges = {(edge.source, edge.target): edge for edge in exact.candidate_edges}
    ambiguous: set[tuple[OrientedNode, OrientedNode]] = set()
    accepted_damage_edges = 0

    for source_index, source in enumerate(topology_reads):
        raw_candidates = _candidate_alignments(
            source.sequence,
            reads,
            raw_ids,
            raw_anchors,
            set(),
            anchor_k=15,
            anchors_per_read=8,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=1,
            minimum_overlap=minimum_overlap,
            minimum_identity=.9,
            minimum_ry_identity=.99,
            position_bits=position_bits,
            target_window=maximum_length,
            diagnostics=raw_diagnostics,
        )
        candidates = _candidate_alignments(
            source.sequence,
            topology_reads,
            topology_ids,
            topology_anchors,
            {source_index},
            anchor_k=15,
            anchors_per_read=8,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=1,
            minimum_overlap=minimum_overlap,
            minimum_identity=.9,
            minimum_ry_identity=.99,
            position_bits=position_bits,
            target_window=maximum_length,
            diagnostics=candidate_diagnostics,
        )
        for candidate in candidates:
            corrections, reason = _damage_compatible_edge(
                reads,
                molecules,
                source.sequence,
                candidate.sequence,
                candidate.offset,
                profile,
                raw_candidates,
            )
            classifications[reason] += 1
            if reason != "damage_compatible":
                continue
            target = (candidate.read_index, candidate.reverse)
            source_length = len(source.sequence)
            target_length = len(candidate.sequence)
            left = max(0, -candidate.offset)
            right = max(0, candidate.offset + target_length - source_length)
            if right and not left:
                edge = ExactEdge(
                    (source_index, False),
                    target,
                    candidate.offset,
                    source_length - candidate.offset,
                    corrections,
                )
            elif left and not right:
                edge = ExactEdge(
                    target,
                    (source_index, False),
                    -candidate.offset,
                    target_length + candidate.offset,
                    tuple(
                        (position - candidate.offset, base)
                        for position, base in corrections
                    ),
                )
            else:
                classifications["non_dovetail_geometry"] += 1
                continue
            if (
                edge.overlap >= len(exact.topology_sequences[edge.source[0]])
                or edge.overlap >= len(exact.topology_sequences[edge.target[0]])
            ):
                classifications["containment_geometry"] += 1
                continue
            if _add_damage_edge(edges, edge, exact.topology_sequences, ambiguous):
                accepted_damage_edges += 1

    reduced = _transitive_reduction(edges, exact.topology_sequences)
    graph = ExactStringGraph(
        exact.topology_sequences,
        tuple(sorted(edges.values(), key=lambda edge: (edge.source, edge.target))),
        tuple(sorted(reduced.values(), key=lambda edge: (edge.source, edge.target))),
        _unitigs(exact.topology_sequences, reduced),
        exact.evidence_placements,
        exact.contained_groups,
        exact.ambiguous_evidence_reads,
    )
    return graph, {
        "candidate_classifications": dict(classifications),
        "candidate_search": asdict(candidate_diagnostics),
        "raw_candidate_search": asdict(raw_diagnostics),
        "damage_directed_edges": accepted_damage_edges,
        "ambiguous_physical_edges": len(ambiguous),
    }


def evidence_pool(
    reads: list[Read],
    graph: ExactStringGraph,
    molecule_ids: list[int] | None = None,
) -> tuple[UnitigEvidencePool, dict]:
    """Compose raw placements onto unitigs without arbitrary allocation."""
    original = ProgressiveSequencePool.from_reads(reads, molecule_ids)
    by_topology: dict[int, list[EvidencePlacement]] = defaultdict(list)
    for placement in graph.evidence_placements:
        by_topology[placement.topology_node].append(placement)

    proposed: dict[int, set[tuple[int, bool, int]]] = defaultdict(set)
    for unitig_index, unitig in enumerate(graph.unitigs):
        for topology_node, topology_reverse, topology_offset in unitig.layout:
            topology_length = len(graph.topology_sequences[topology_node])
            for placement in by_topology[topology_node]:
                raw_length = len(reads[placement.raw_index].sequence)
                if topology_reverse:
                    offset = (
                        topology_offset
                        + topology_length
                        - placement.offset
                        - raw_length
                    )
                    reverse = not placement.reverse
                else:
                    offset = topology_offset + placement.offset
                    reverse = placement.reverse
                proposed[placement.raw_index].add(
                    (unitig_index, reverse, offset)
                )

    unique = {
        raw_index: next(iter(placements))
        for raw_index, placements in proposed.items()
        if len(placements) == 1
    }
    placements_by_unitig: dict[int, list[RawPlacement]] = defaultdict(list)
    for raw_index, (unitig_index, reverse, offset) in unique.items():
        placements_by_unitig[unitig_index].append(RawPlacement(
            raw_index,
            offset,
            reverse,
            0,
            len(reads[raw_index].sequence),
        ))

    derived = []
    for unitig_index, unitig in enumerate(graph.unitigs):
        placements = tuple(sorted(
            placements_by_unitig.get(unitig_index, ()),
            key=lambda item: (item.read_index, item.offset, item.reverse),
        ))
        contributing = frozenset(
            original.records[placement.read_index].molecule_id
            for placement in placements
        )
        contig = Read(f"exact_unitig_{unitig_index + 1}", unitig.sequence)
        derived.append(SequenceRecord(
            index=unitig_index,
            molecule_id=-(unitig_index + 1),
            raw=contig,
            current=contig,
            state=SequenceState.EXTENDED_CONTIG,
            generation=1,
            contributing_molecules=contributing,
            raw_placements=placements,
        ))
    ambiguous = sorted(
        raw_index for raw_index, placements in proposed.items()
        if len(placements) != 1
    )
    missing = sorted(set(range(len(reads))) - set(proposed))
    return UnitigEvidencePool(original.records, tuple(derived)), {
        "placed_raw_reads": len(unique),
        "ambiguous_raw_reads": len(ambiguous),
        "missing_raw_reads": len(missing),
        "ambiguous_raw_indices": ambiguous,
        "missing_raw_indices": missing,
    }
