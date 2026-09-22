"""Truth-blind phase evidence from a complete accepted overlap graph."""

from collections import defaultdict, deque
from dataclasses import dataclass

from anvaya.overlap_progressive import ProgressiveSequencePool
from anvaya.phase_blocks import ComponentPhaseEvidence, phase_placements
from anvaya.raw_consensus import DamageProfile, RawPlacement
from anvaya.sequences import reverse_complement


@dataclass(frozen=True, slots=True)
class GraphNodePlacement:
    node_index: int
    reverse: bool
    offset: int


@dataclass(frozen=True, slots=True)
class GraphCoordinateConflict:
    source: tuple[int, bool]
    target: tuple[int, bool]
    assigned_offset: int
    proposed_offset: int


@dataclass(frozen=True, slots=True)
class GraphComponentCoordinates:
    physical_nodes: tuple[int, ...]
    placements: tuple[GraphNodePlacement, ...]
    span: int
    conflicts: tuple[GraphCoordinateConflict, ...]
    orientation_conflicts: tuple[int, ...]

    @property
    def consistent(self):
        return not self.conflicts and not self.orientation_conflicts


@dataclass(frozen=True, slots=True)
class GraphComponentPhaseEvidence:
    coordinates: GraphComponentCoordinates
    phase: ComponentPhaseEvidence | None


def graph_component_coordinates(sequences, edges):
    """Solve every accepted overlap constraint without selecting graph paths.

    Node and target coordinates are 0-based. An edge constrains
    ``target_offset = source_offset + shift``. Bidirected mirror components are
    reported once. Any inconsistent cycle or use of both orientations of one
    physical node makes the component unresolved.
    """
    adjacency = defaultdict(list)
    for edge in edges.values():
        adjacency[edge.source].append((edge.target, edge.shift, edge.source, edge.target))
        adjacency[edge.target].append((edge.source, -edge.shift, edge.source, edge.target))

    unseen = set(adjacency)
    oriented_components = []
    while unseen:
        root = min(unseen)
        component = set()
        queue = [root]
        while queue:
            node = queue.pop()
            if node in component:
                continue
            component.add(node)
            unseen.discard(node)
            queue.extend(target for target, *_ in adjacency[node])
        oriented_components.append(component)

    represented = {node[0] for node in adjacency}
    oriented_components.extend({(index, False)} for index in range(len(sequences))
                               if index not in represented)

    results = []
    seen_orientations = set()
    for component in sorted(oriented_components, key=lambda value: sorted(value)):
        signature = frozenset(component)
        if signature in seen_orientations:
            continue
        seen_orientations.add(signature)
        seen_orientations.add(frozenset(
            (index, not reverse) for index, reverse in component
        ))
        physical = tuple(sorted({node[0] for node in component}))
        root = min(component)
        offsets = {root: 0}
        conflicts = set()
        queue = deque((root,))
        while queue:
            source = queue.popleft()
            for target, delta, edge_source, edge_target in adjacency[source]:
                if target not in component:
                    continue
                proposed = offsets[source] + delta
                if target not in offsets:
                    offsets[target] = proposed
                    queue.append(target)
                elif offsets[target] != proposed:
                    conflicts.add(GraphCoordinateConflict(
                        edge_source, edge_target, offsets[target], proposed,
                    ))

        minimum = min(offsets.values())
        placements = tuple(sorted(
            (GraphNodePlacement(index, reverse, offset - minimum)
             for (index, reverse), offset in offsets.items()),
            key=lambda item: (item.offset, item.node_index, item.reverse),
        ))
        span = max(item.offset + len(sequences[item.node_index]) for item in placements)
        orientations = defaultdict(set)
        for placement in placements:
            orientations[placement.node_index].add(placement.reverse)
        orientation_conflicts = tuple(sorted(
            index for index, values in orientations.items() if len(values) > 1
        ))
        results.append(GraphComponentCoordinates(
            physical,
            placements,
            span,
            tuple(sorted(conflicts, key=lambda item: (
                item.source, item.target, item.assigned_offset, item.proposed_offset,
            ))),
            orientation_conflicts,
        ))
    return tuple(results)
def phase_graph_components(reads, molecule_ids, sequences, groups, edges, profile):
    """Project all coordinate-consistent graph alternatives into phase evidence.

    ``groups`` maps each sequence to ``(raw_read_index, raw_reverse)`` pairs, as
    produced before string-graph path selection. Coordinate-conflicting graph
    components return ``phase=None`` and retain their conflicts for audit.
    """
    if len(reads) != len(molecule_ids):
        raise ValueError("molecule IDs must align one-to-one with reads")
    pool = ProgressiveSequencePool.from_reads(reads, molecule_ids)
    results = []
    for coordinates in graph_component_coordinates(sequences, edges):
        if not coordinates.consistent:
            results.append(GraphComponentPhaseEvidence(coordinates, None))
            continue
        placements = []
        for node in coordinates.placements:
            for raw_index, raw_reverse in groups[sequences[node.node_index]]:
                read_length = len(reads[raw_index].sequence)
                placements.append(RawPlacement(
                    raw_index,
                    node.offset,
                    node.reverse != raw_reverse,
                    0,
                    read_length,
                ))
                if sequences[node.node_index] == reverse_complement(sequences[node.node_index]):
                    placements.append(RawPlacement(
                        raw_index,
                        node.offset,
                        node.reverse == raw_reverse,
                        0,
                        read_length,
                    ))
        evidence = phase_placements(pool, coordinates.span, placements, profile)
        results.append(GraphComponentPhaseEvidence(coordinates, evidence))
    return tuple(results)
