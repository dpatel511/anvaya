"""Orientation-aware de Bruijn graph construction and unitig extraction."""

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field

from anvaya.kmers import kmers
from anvaya.metrics import GraphSummary
from anvaya.sequences import reverse_complement

Handle = int
PhysicalEdge = tuple[Handle, Handle]


@dataclass(slots=True)
class StrandSupport:
    """Observation counts relative to a canonical k-mer."""

    forward: int = 0
    reverse: int = 0
    palindromic: int = 0

    @property
    def total(self) -> int:
        return self.forward + self.reverse + self.palindromic


@dataclass(slots=True)
class BidirectedDeBruijnGraph:
    """A graph of canonical nodes traversed through oriented handles."""

    sequences: list[str] = field(default_factory=list)
    adjacency: dict[Handle, Counter[Handle]] = field(default_factory=dict)
    edge_support: dict[PhysicalEdge, StrandSupport] = field(default_factory=dict)


def flip_handle(handle: Handle) -> Handle:
    """Return the opposite orientation of a handle."""
    return handle ^ 1


def physical_edge(source: Handle, target: Handle) -> PhysicalEdge:
    """Return one identity for both traversals of a physical edge."""
    reverse_edge = (flip_handle(target), flip_handle(source))
    return min((source, target), reverse_edge)


def _canonical_orientation(sequence: str) -> tuple[str, bool]:
    reverse = reverse_complement(sequence)
    if sequence <= reverse:
        return sequence, False
    return reverse, True


def oriented_sequence(graph: BidirectedDeBruijnGraph, handle: Handle) -> str:
    """Return the node sequence spelled by an oriented handle."""
    sequence = graph.sequences[handle >> 1]
    if handle & 1:
        return reverse_complement(sequence)
    return sequence


def _node_handle(
    sequence: str,
    node_ids: dict[str, int],
    graph: BidirectedDeBruijnGraph,
) -> Handle:
    canonical, is_reverse = _canonical_orientation(sequence)
    node_id = node_ids.get(canonical)
    if node_id is None:
        node_id = len(graph.sequences)
        node_ids[canonical] = node_id
        graph.sequences.append(canonical)
        graph.adjacency[node_id << 1] = Counter()
        graph.adjacency[(node_id << 1) | 1] = Counter()
    return (node_id << 1) | is_reverse


def build_bidirected_dbg(
    reads: Sequence[str], k: int, min_count: int = 1
) -> BidirectedDeBruijnGraph:
    """Build a canonical graph while retaining strand-specific support."""
    if not reads:
        raise ValueError("reads must not be empty")
    if not isinstance(k, int) or isinstance(k, bool):
        raise TypeError("k must be an integer")
    if k < 2:
        raise ValueError("k must be at least 2 for graph construction")
    if not isinstance(min_count, int) or isinstance(min_count, bool):
        raise TypeError("min_count must be an integer")
    if min_count < 1:
        raise ValueError("min_count must be at least 1")

    observed: dict[str, StrandSupport] = {}
    for read in reads:
        if len(read) < k:
            raise ValueError("each read must be at least k bases long")
        for kmer in kmers(read, k):
            if "N" in kmer:
                continue
            reverse = reverse_complement(kmer)
            canonical = min(kmer, reverse)
            support = observed.setdefault(canonical, StrandSupport())
            if kmer == reverse:
                support.palindromic += 1
            elif kmer == canonical:
                support.forward += 1
            else:
                support.reverse += 1

    graph = BidirectedDeBruijnGraph()
    node_ids: dict[str, int] = {}

    while observed:
        kmer, support = observed.popitem()
        if support.total < min_count:
            continue

        source = _node_handle(kmer[:-1], node_ids, graph)
        target = _node_handle(kmer[1:], node_ids, graph)
        reverse_source = flip_handle(target)
        reverse_target = flip_handle(source)
        edge = physical_edge(source, target)

        graph.adjacency[source][target] += support.total
        if (source, target) != (reverse_source, reverse_target):
            graph.adjacency[reverse_source][reverse_target] += support.total
        graph.edge_support[edge] = support

    return graph


def bidirected_node_degrees(
    graph: BidirectedDeBruijnGraph,
) -> tuple[dict[Handle, int], dict[Handle, int]]:
    """Return incoming and outgoing degrees for every oriented handle."""
    in_degrees = dict.fromkeys(graph.adjacency, 0)
    out_degrees = {
        handle: len(successors) for handle, successors in graph.adjacency.items()
    }
    for successors in graph.adjacency.values():
        for successor in successors:
            in_degrees[successor] += 1
    return in_degrees, out_degrees


def _spell_path(graph: BidirectedDeBruijnGraph, path: list[Handle]) -> str:
    sequence = oriented_sequence(graph, path[0])
    return sequence + "".join(
        oriented_sequence(graph, handle)[-1] for handle in path[1:]
    )


def extract_bidirected_unitigs(graph: BidirectedDeBruijnGraph) -> list[str]:
    """Return one sequence for each maximal non-branching physical path."""
    in_degrees, out_degrees = bidirected_node_degrees(graph)
    visited: set[PhysicalEdge] = set()
    unitigs: list[str] = []

    for start in sorted(graph.adjacency):
        if in_degrees[start] == 1 and out_degrees[start] == 1:
            continue

        for successor in sorted(graph.adjacency[start]):
            edge = physical_edge(start, successor)
            if edge in visited:
                continue

            path = [start, successor]
            visited.add(edge)
            current = successor

            while in_degrees[current] == 1 and out_degrees[current] == 1:
                following = next(iter(graph.adjacency[current]))
                edge = physical_edge(current, following)
                if edge in visited:
                    break
                path.append(following)
                visited.add(edge)
                current = following

            unitigs.append(_spell_path(graph, path))

    for start in sorted(graph.adjacency):
        for successor in sorted(graph.adjacency[start]):
            edge = physical_edge(start, successor)
            if edge in visited:
                continue

            path = [start, successor]
            visited.add(edge)
            current = successor

            while current != start:
                following = next(iter(graph.adjacency[current]))
                edge = physical_edge(current, following)
                if edge in visited:
                    break
                path.append(following)
                visited.add(edge)
                current = following

            unitigs.append(_spell_path(graph, path))

    return unitigs


def summarize_bidirected_graph(
    graph: BidirectedDeBruijnGraph, unitigs: list[str] | None = None
) -> GraphSummary:
    """Return topology and support counts for a bidirected graph."""
    in_degrees, out_degrees = bidirected_node_degrees(graph)
    if unitigs is None:
        unitigs = extract_bidirected_unitigs(graph)

    source_nodes = {
        handle >> 1
        for handle in graph.adjacency
        if in_degrees[handle] == 0 and out_degrees[handle] > 0
    }
    sink_nodes = {
        handle >> 1
        for handle in graph.adjacency
        if in_degrees[handle] > 0 and out_degrees[handle] == 0
    }
    branching_nodes = {
        handle >> 1
        for handle in graph.adjacency
        if in_degrees[handle] > 1 or out_degrees[handle] > 1
    }

    return GraphSummary(
        nodes=len(graph.sequences),
        edges=len(graph.edge_support),
        observations=sum(support.total for support in graph.edge_support.values()),
        sources=len(source_nodes),
        sinks=len(sink_nodes),
        branching_nodes=len(branching_nodes),
        unitigs=len(unitigs),
    )
