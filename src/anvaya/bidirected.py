"""Compact orientation-aware de Bruijn graph construction and traversal."""

from array import array
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field

from anvaya.metrics import GraphSummary
from anvaya.sequences import normalize_dna

Handle = int

_BASE_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}
_COUNT_BITS = 32
_COUNT_MASK = (1 << _COUNT_BITS) - 1
_REVERSE_INCREMENT = 1 << _COUNT_BITS
_PALINDROMIC_INCREMENT = 1 << (2 * _COUNT_BITS)
_MAX_OUT_DEGREE = 4


@dataclass(slots=True, frozen=True)
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
    """Canonical nodes and physical edges stored in compact integer arrays."""

    node_length: int
    node_codes: list[int] = field(default_factory=list)
    edge_sources: array = field(default_factory=lambda: array("Q"))
    edge_targets: array = field(default_factory=lambda: array("Q"))
    forward_support: array = field(default_factory=lambda: array("I"))
    reverse_support: array = field(default_factory=lambda: array("I"))
    palindromic_support: array = field(default_factory=lambda: array("I"))
    out_edges: array = field(default_factory=lambda: array("i"))
    out_degrees: bytearray = field(default_factory=bytearray)
    observations: int = 0

    @property
    def node_count(self) -> int:
        return len(self.node_codes)

    @property
    def edge_count(self) -> int:
        return len(self.edge_sources)

    def strand_support(self, edge_id: int) -> StrandSupport:
        """Return strand-specific observation counts for one physical edge."""
        return StrandSupport(
            self.forward_support[edge_id],
            self.reverse_support[edge_id],
            self.palindromic_support[edge_id],
        )


def flip_handle(handle: Handle) -> Handle:
    """Return the opposite orientation of a handle."""
    return handle ^ 1


def _reverse_complement_code(code: int, length: int) -> int:
    reverse = 0
    for _ in range(length):
        reverse = (reverse << 2) | (3 - (code & 3))
        code >>= 2
    return reverse


def _decode_dna(code: int, length: int) -> str:
    bases = "ACGT"
    sequence = ["A"] * length
    for index in range(length - 1, -1, -1):
        sequence[index] = bases[code & 3]
        code >>= 2
    return "".join(sequence)


def _encoded_kmers(sequence: str, k: int) -> Iterator[tuple[int, int]]:
    """Yield canonical k-mer codes and their strand classification."""
    normalized = normalize_dna(sequence)
    mask = (1 << (2 * k)) - 1
    reverse_shift = 2 * (k - 1)
    forward = 0
    reverse = 0
    valid_bases = 0

    for base in normalized:
        if base == "N":
            forward = 0
            reverse = 0
            valid_bases = 0
            continue

        bits = _BASE_BITS[base]
        forward = ((forward << 2) | bits) & mask
        reverse = (reverse >> 2) | ((3 - bits) << reverse_shift)
        valid_bases += 1
        if valid_bases < k:
            continue

        if forward < reverse:
            yield forward, 0
        elif reverse < forward:
            yield reverse, 1
        else:
            yield forward, 2


def _unpack_support(packed: int) -> StrandSupport:
    return StrandSupport(
        packed & _COUNT_MASK,
        (packed >> _COUNT_BITS) & _COUNT_MASK,
        packed >> (2 * _COUNT_BITS),
    )


def _node_handle(
    code: int,
    reverse_code: int,
    node_ids: dict[int, int],
    graph: BidirectedDeBruijnGraph,
) -> Handle:
    is_reverse = reverse_code < code
    canonical = reverse_code if is_reverse else code
    node_id = node_ids.get(canonical)
    if node_id is None:
        node_id = len(graph.node_codes)
        node_ids[canonical] = node_id
        graph.node_codes.append(canonical)
    return (node_id << 1) | is_reverse


def _append_arc(graph: BidirectedDeBruijnGraph, source: Handle, edge_id: int) -> None:
    degree = graph.out_degrees[source]
    if degree >= _MAX_OUT_DEGREE:
        raise ValueError("oriented de Bruijn graph node has more than four outgoing edges")
    graph.out_edges[source * _MAX_OUT_DEGREE + degree] = edge_id
    graph.out_degrees[source] = degree + 1


def build_bidirected_dbg(
    reads: Sequence[str], k: int, min_count: int = 1
) -> BidirectedDeBruijnGraph:
    """Build a compact canonical graph with strand-specific edge support."""
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

    observed: dict[int, int] = {}
    for read in reads:
        if len(read) < k:
            raise ValueError("each read must be at least k bases long")
        for canonical, orientation in _encoded_kmers(read, k):
            if orientation == 0:
                increment = 1
            elif orientation == 1:
                increment = _REVERSE_INCREMENT
            else:
                increment = _PALINDROMIC_INCREMENT
            observed[canonical] = observed.get(canonical, 0) + increment

    graph = BidirectedDeBruijnGraph(node_length=k - 1)
    node_ids: dict[int, int] = {}
    node_mask = (1 << (2 * (k - 1))) - 1

    while observed:
        kmer_code, packed = observed.popitem()
        support = _unpack_support(packed)
        if support.total < min_count:
            continue

        reverse_code = _reverse_complement_code(kmer_code, k)
        prefix = kmer_code >> 2
        suffix = kmer_code & node_mask
        source = _node_handle(prefix, reverse_code & node_mask, node_ids, graph)
        target = _node_handle(suffix, reverse_code >> 2, node_ids, graph)

        graph.edge_sources.append(source)
        graph.edge_targets.append(target)
        graph.forward_support.append(support.forward)
        graph.reverse_support.append(support.reverse)
        graph.palindromic_support.append(support.palindromic)
        graph.observations += support.total

    handle_count = graph.node_count * 2
    graph.out_degrees = bytearray(handle_count)
    graph.out_edges = array("i", [-1]) * (handle_count * _MAX_OUT_DEGREE)

    for edge_id in range(graph.edge_count):
        source = graph.edge_sources[edge_id]
        target = graph.edge_targets[edge_id]
        _append_arc(graph, source, edge_id)
        reverse_source = flip_handle(target)
        if reverse_source != source:
            _append_arc(graph, reverse_source, edge_id)

    return graph


def oriented_sequence(graph: BidirectedDeBruijnGraph, handle: Handle) -> str:
    """Return the node sequence spelled by an oriented handle."""
    code = graph.node_codes[handle >> 1]
    if handle & 1:
        code = _reverse_complement_code(code, graph.node_length)
    return _decode_dna(code, graph.node_length)


def _successor(
    graph: BidirectedDeBruijnGraph, handle: Handle, edge_id: int
) -> Handle:
    source = graph.edge_sources[edge_id]
    target = graph.edge_targets[edge_id]
    if handle == source:
        return target
    return flip_handle(source)


def _outgoing_edges(
    graph: BidirectedDeBruijnGraph, handle: Handle
) -> Iterator[int]:
    offset = handle * _MAX_OUT_DEGREE
    for index in range(graph.out_degrees[handle]):
        yield graph.out_edges[offset + index]


def _spell_path(graph: BidirectedDeBruijnGraph, path: list[Handle]) -> str:
    sequence = oriented_sequence(graph, path[0])
    return sequence + "".join(
        oriented_sequence(graph, handle)[-1] for handle in path[1:]
    )


def extract_bidirected_unitigs(graph: BidirectedDeBruijnGraph) -> list[str]:
    """Return one sequence for each maximal non-branching physical path."""
    visited = bytearray(graph.edge_count)
    unitigs: list[str] = []
    handle_count = graph.node_count * 2

    for start in range(handle_count):
        in_degree = graph.out_degrees[flip_handle(start)]
        out_degree = graph.out_degrees[start]
        if in_degree == 1 and out_degree == 1:
            continue

        for edge_id in _outgoing_edges(graph, start):
            if visited[edge_id]:
                continue

            successor = _successor(graph, start, edge_id)
            path = [start, successor]
            visited[edge_id] = 1
            current = successor

            while (
                graph.out_degrees[flip_handle(current)] == 1
                and graph.out_degrees[current] == 1
            ):
                following_edge = graph.out_edges[current * _MAX_OUT_DEGREE]
                if visited[following_edge]:
                    break
                current = _successor(graph, current, following_edge)
                path.append(current)
                visited[following_edge] = 1

            unitigs.append(_spell_path(graph, path))

    for start in range(handle_count):
        for edge_id in _outgoing_edges(graph, start):
            if visited[edge_id]:
                continue

            successor = _successor(graph, start, edge_id)
            path = [start, successor]
            visited[edge_id] = 1
            current = successor

            while current != start:
                following_edge = graph.out_edges[current * _MAX_OUT_DEGREE]
                if visited[following_edge]:
                    break
                current = _successor(graph, current, following_edge)
                path.append(current)
                visited[following_edge] = 1

            unitigs.append(_spell_path(graph, path))

    return unitigs


def summarize_bidirected_graph(
    graph: BidirectedDeBruijnGraph, unitigs: list[str] | None = None
) -> GraphSummary:
    """Return topology and support counts for a bidirected graph."""
    if unitigs is None:
        unitigs = extract_bidirected_unitigs(graph)

    source_nodes = bytearray(graph.node_count)
    sink_nodes = bytearray(graph.node_count)
    branching_nodes = bytearray(graph.node_count)
    for handle in range(graph.node_count * 2):
        in_degree = graph.out_degrees[flip_handle(handle)]
        out_degree = graph.out_degrees[handle]
        node_id = handle >> 1
        if in_degree == 0 and out_degree > 0:
            source_nodes[node_id] = 1
        if in_degree > 0 and out_degree == 0:
            sink_nodes[node_id] = 1
        if in_degree > 1 or out_degree > 1:
            branching_nodes[node_id] = 1

    return GraphSummary(
        nodes=graph.node_count,
        edges=graph.edge_count,
        observations=graph.observations,
        sources=sum(source_nodes),
        sinks=sum(sink_nodes),
        branching_nodes=sum(branching_nodes),
        unitigs=len(unitigs),
    )
