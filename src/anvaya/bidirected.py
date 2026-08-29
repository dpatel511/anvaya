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


@dataclass(slots=True, frozen=True)
class EndSupport:
    """Read-end evidence normalized to the canonical edge orientation."""

    left: tuple[int, ...]
    right: tuple[int, ...]
    internal: int = 0
    ambiguous_terminal: int = 0


@dataclass(slots=True, frozen=True)
class MoleculeEndLink:
    """One terminal edge observation linked to its read and molecule."""

    read_index: int
    molecule_id: int
    end: str
    distance: int
    orientation: int
    base_quality: int | None = None


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
    end_window: int = 0
    left_end_support: array = field(default_factory=lambda: array("I"))
    right_end_support: array = field(default_factory=lambda: array("I"))
    internal_support: array = field(default_factory=lambda: array("I"))
    ambiguous_end_support: array = field(default_factory=lambda: array("I"))
    molecule_links_collected: bool = False
    molecule_link_offsets: array = field(default_factory=lambda: array("Q"))
    molecule_link_tokens: array = field(default_factory=lambda: array("Q"))
    molecule_link_qualities: bytearray = field(default_factory=bytearray)
    read_molecule_ids: array = field(default_factory=lambda: array("Q"))
    out_edges: array = field(default_factory=lambda: array("i"))
    out_degrees: bytearray = field(default_factory=bytearray)
    active_edge_count: int = 0
    observations: int = 0
    terminal_observations: int = 0

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

    def end_support(self, edge_id: int) -> EndSupport:
        """Return terminal and internal support for one physical edge."""
        if self.end_window == 0:
            raise ValueError("read-end evidence was not collected")
        offset = edge_id * self.end_window
        end = offset + self.end_window
        return EndSupport(
            left=tuple(self.left_end_support[offset:end]),
            right=tuple(self.right_end_support[offset:end]),
            internal=self.internal_support[edge_id],
            ambiguous_terminal=self.ambiguous_end_support[edge_id],
        )

    def molecule_end_links(
        self,
        edge_id: int,
    ) -> tuple[MoleculeEndLink, ...]:
        """Return terminal observations linked to source-read indices."""
        if not self.molecule_links_collected:
            raise ValueError("molecule links were not collected")
        start = self.molecule_link_offsets[edge_id]
        end = self.molecule_link_offsets[edge_id + 1]
        links: list[MoleculeEndLink] = []
        for index in range(start, end):
            token = self.molecule_link_tokens[index]
            read_side_orientation, distance = divmod(
                token,
                self.end_window,
            )
            read_index, side_orientation = divmod(
                read_side_orientation,
                4,
            )
            side, orientation = divmod(side_orientation, 2)
            quality = self.molecule_link_qualities[index]
            links.append(
                MoleculeEndLink(
                    read_index=read_index,
                    molecule_id=self.read_molecule_ids[read_index],
                    end="left" if side == 0 else "right",
                    distance=distance,
                    orientation=orientation,
                    base_quality=None if quality == 255 else quality,
                )
            )
        return tuple(links)

    def molecule_end_distances(
        self,
        edge_id: int,
        end_name: str,
    ) -> Iterator[int]:
        """Yield retained terminal distances without allocating link objects."""
        if not self.molecule_links_collected:
            raise ValueError("molecule links were not collected")
        if end_name not in {"left", "right"}:
            raise ValueError("end_name must be 'left' or 'right'")
        expected_side = 0 if end_name == "left" else 1
        start = self.molecule_link_offsets[edge_id]
        end = self.molecule_link_offsets[edge_id + 1]
        for token in self.molecule_link_tokens[start:end]:
            read_side_orientation, distance = divmod(
                token,
                self.end_window,
            )
            side_orientation = read_side_orientation % 4
            if side_orientation // 2 == expected_side:
                yield distance


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
    for _, canonical, orientation in _positioned_encoded_kmers(sequence, k):
        yield canonical, orientation


def _positioned_encoded_kmers(
    sequence: str, k: int
) -> Iterator[tuple[int, int, int]]:
    """Yield read offsets, canonical k-mers, and strand classifications."""
    normalized = normalize_dna(sequence)
    mask = (1 << (2 * k)) - 1
    reverse_shift = 2 * (k - 1)
    forward = 0
    reverse = 0
    valid_bases = 0

    for index, base in enumerate(normalized):
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
            yield index - k + 1, forward, 0
        elif reverse < forward:
            yield index - k + 1, reverse, 1
        else:
            yield index - k + 1, forward, 2


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
    reads: Sequence[str],
    k: int,
    min_count: int = 1,
    end_window: int = 0,
    track_molecule_links: bool = False,
    read_qualities: Sequence[Sequence[int] | None] | None = None,
    read_molecule_ids: Sequence[int] | None = None,
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
    if not isinstance(end_window, int) or isinstance(end_window, bool):
        raise TypeError("end_window must be an integer")
    if end_window < 0:
        raise ValueError("end_window must not be negative")
    if not isinstance(track_molecule_links, bool):
        raise TypeError("track_molecule_links must be a boolean")
    if track_molecule_links and end_window == 0:
        raise ValueError("molecule links require read-end evidence")
    if read_qualities is not None and len(read_qualities) != len(reads):
        raise ValueError("read qualities must align one-to-one with reads")
    if read_molecule_ids is not None and len(read_molecule_ids) != len(reads):
        raise ValueError("molecule IDs must align one-to-one with reads")
    if read_molecule_ids is not None and any(
        not isinstance(molecule_id, int)
        or isinstance(molecule_id, bool)
        or molecule_id < 0
        for molecule_id in read_molecule_ids
    ):
        raise ValueError("molecule IDs must be non-negative integers")
    if read_qualities is not None:
        for read, qualities in zip(reads, read_qualities, strict=True):
            if qualities is not None and len(qualities) != len(read):
                raise ValueError("read quality lengths must match read lengths")
            if qualities is not None and any(
                not isinstance(quality, int) or isinstance(quality, bool)
                for quality in qualities
            ):
                raise TypeError("read qualities must be integers")
            if qualities is not None and any(
                quality < 0 or quality > 254 for quality in qualities
            ):
                raise ValueError("read qualities must be between 0 and 254")

    observed: dict[int, int] = {}
    terminal: dict[int, list[int]] = {}
    for read in reads:
        if len(read) < k:
            raise ValueError("each read must be at least k bases long")
        last_start = len(read) - k
        for start, canonical, orientation in _positioned_encoded_kmers(read, k):
            if orientation == 0:
                increment = 1
            elif orientation == 1:
                increment = _REVERSE_INCREMENT
            else:
                increment = _PALINDROMIC_INCREMENT
            observed[canonical] = observed.get(canonical, 0) + increment

            if end_window == 0:
                continue
            left_distance = start
            right_distance = last_start - start
            if orientation == 1:
                left_distance, right_distance = right_distance, left_distance
            is_terminal = (
                left_distance < end_window or right_distance < end_window
            )
            if not is_terminal:
                continue

            evidence = terminal.get(canonical)
            if evidence is None:
                evidence = [0] * (2 * end_window + 2)
                terminal[canonical] = evidence
            if orientation == 2:
                evidence[-1] += 1
                continue
            if left_distance < end_window:
                evidence[left_distance] += 1
            if right_distance < end_window:
                evidence[end_window + right_distance] += 1
            evidence[-2] += 1

    graph = BidirectedDeBruijnGraph(node_length=k - 1, end_window=end_window)
    node_ids: dict[int, int] = {}
    edge_ids_by_code: dict[int, int] = {}
    node_mask = (1 << (2 * (k - 1))) - 1

    while observed:
        kmer_code, packed = observed.popitem()
        support = _unpack_support(packed)
        evidence = terminal.pop(kmer_code, None)
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
        if track_molecule_links:
            edge_ids_by_code[kmer_code] = graph.edge_count - 1
        if end_window:
            if evidence is None:
                evidence = [0] * (2 * end_window + 2)
            graph.left_end_support.extend(evidence[:end_window])
            graph.right_end_support.extend(evidence[end_window : 2 * end_window])
            terminal_support = evidence[-2] + evidence[-1]
            graph.internal_support.append(support.total - terminal_support)
            graph.ambiguous_end_support.append(evidence[-1])
            graph.terminal_observations += terminal_support

    if track_molecule_links:
        graph.read_molecule_ids = array(
            "Q",
            range(len(reads)) if read_molecule_ids is None else read_molecule_ids,
        )
        _collect_molecule_links(
            graph,
            reads,
            k,
            edge_ids_by_code,
            read_qualities,
        )

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

    graph.active_edge_count = graph.edge_count

    return graph


def _collect_molecule_links(
    graph: BidirectedDeBruijnGraph,
    reads: Sequence[str],
    k: int,
    edge_ids_by_code: dict[int, int],
    read_qualities: Sequence[Sequence[int] | None] | None,
) -> None:
    """Store retained terminal observations in compact edge-indexed arrays."""
    counts = array("Q", [0]) * graph.edge_count

    def observations():
        for read_index, read in enumerate(reads):
            qualities = (
                None if read_qualities is None else read_qualities[read_index]
            )
            last_start = len(read) - k
            for start, canonical, orientation in _positioned_encoded_kmers(
                read,
                k,
            ):
                edge_id = edge_ids_by_code.get(canonical)
                if edge_id is None or orientation == 2:
                    continue
                left_distance = start
                right_distance = last_start - start
                if orientation == 1:
                    left_distance, right_distance = (
                        right_distance,
                        left_distance,
                    )
                left_base_index = start if orientation == 1 else start + k - 1
                right_base_index = start + k - 1 if orientation == 1 else start
                if left_distance < graph.end_window:
                    quality = (
                        None
                        if qualities is None
                        else qualities[left_base_index]
                    )
                    yield (
                        edge_id,
                        read_index,
                        0,
                        orientation,
                        left_distance,
                        quality,
                    )
                if right_distance < graph.end_window:
                    quality = (
                        None
                        if qualities is None
                        else qualities[right_base_index]
                    )
                    yield (
                        edge_id,
                        read_index,
                        1,
                        orientation,
                        right_distance,
                        quality,
                    )

    for edge_id, _, _, _, _, _ in observations():
        counts[edge_id] += 1

    graph.molecule_link_offsets.append(0)
    for count in counts:
        graph.molecule_link_offsets.append(
            graph.molecule_link_offsets[-1] + count
        )
    graph.molecule_link_tokens = array(
        "Q",
        [0],
    ) * graph.molecule_link_offsets[-1]
    graph.molecule_link_qualities = bytearray(
        [255]
    ) * graph.molecule_link_offsets[-1]
    cursors = graph.molecule_link_offsets[:-1]
    for (
        edge_id,
        read_index,
        side,
        orientation,
        distance,
        quality,
    ) in observations():
        token = (
            (read_index * 4 + side * 2 + orientation)
            * graph.end_window
            + distance
        )
        index = cursors[edge_id]
        graph.molecule_link_tokens[index] = token
        graph.molecule_link_qualities[index] = (
            255 if quality is None else quality
        )
        cursors[edge_id] += 1
    graph.molecule_links_collected = True


def collect_target_edge_links(
    graph: BidirectedDeBruijnGraph,
    reads: Sequence[str],
    edge_ids: Sequence[int],
    read_qualities: Sequence[Sequence[int] | None] | None = None,
    read_molecule_ids: Sequence[int] | None = None,
) -> dict[int, dict[int, int | None]]:
    """Aggregate selected edge observations once per source molecule."""
    if read_qualities is not None and len(read_qualities) != len(reads):
        raise ValueError("read qualities must align one-to-one with reads")
    if read_molecule_ids is not None and len(read_molecule_ids) != len(reads):
        raise ValueError("molecule IDs must align one-to-one with reads")

    k = graph.node_length + 1
    code_to_edge: dict[int, int] = {}
    for edge_id in set(edge_ids):
        source = graph.edge_sources[edge_id]
        sequence = oriented_sequence(graph, source) + _oriented_last_base(
            graph, graph.edge_targets[edge_id]
        )
        _, code, orientation = next(_positioned_encoded_kmers(sequence, k))
        if orientation != 2:
            code_to_edge[code] = edge_id

    links: dict[int, dict[int, int | None]] = {
        edge_id: {} for edge_id in code_to_edge.values()
    }
    for read_index, read in enumerate(reads):
        qualities = None if read_qualities is None else read_qualities[read_index]
        molecule_id = (
            read_index
            if read_molecule_ids is None
            else read_molecule_ids[read_index]
        )
        for start, code, orientation in _positioned_encoded_kmers(read, k):
            edge_id = code_to_edge.get(code)
            if edge_id is None or orientation == 2:
                continue
            base_index = start + k - 1 if orientation == 0 else start
            quality = None if qualities is None else qualities[base_index]
            previous = links[edge_id].get(molecule_id)
            if previous is None or (quality is not None and quality > previous):
                links[edge_id][molecule_id] = quality
    return links


def oriented_sequence(graph: BidirectedDeBruijnGraph, handle: Handle) -> str:
    """Return the node sequence spelled by an oriented handle."""
    code = graph.node_codes[handle >> 1]
    if handle & 1:
        code = _reverse_complement_code(code, graph.node_length)
    return _decode_dna(code, graph.node_length)


def _oriented_last_base(
    graph: BidirectedDeBruijnGraph, handle: Handle
) -> str:
    """Return the final base of an oriented node without decoding it."""
    code = graph.node_codes[handle >> 1]
    if handle & 1:
        bits = (code >> (2 * (graph.node_length - 1))) & 3
        return "ACGT"[3 - bits]
    return "ACGT"[code & 3]


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


def _edge_support_total(graph: BidirectedDeBruijnGraph, edge_id: int) -> int:
    return (
        graph.forward_support[edge_id]
        + graph.reverse_support[edge_id]
        + graph.palindromic_support[edge_id]
    )


def _remove_arc(
    graph: BidirectedDeBruijnGraph, source: Handle, edge_id: int
) -> None:
    offset = source * _MAX_OUT_DEGREE
    degree = graph.out_degrees[source]
    for index in range(degree):
        if graph.out_edges[offset + index] != edge_id:
            continue
        last_index = offset + degree - 1
        graph.out_edges[offset + index] = graph.out_edges[last_index]
        graph.out_edges[last_index] = -1
        graph.out_degrees[source] = degree - 1
        return
    raise ValueError("edge is not active from the requested handle")


def _deactivate_edge(graph: BidirectedDeBruijnGraph, edge_id: int) -> int:
    source = graph.edge_sources[edge_id]
    target = graph.edge_targets[edge_id]
    _remove_arc(graph, source, edge_id)
    reverse_source = flip_handle(target)
    if reverse_source != source:
        _remove_arc(graph, reverse_source, edge_id)

    support = _edge_support_total(graph, edge_id)
    graph.active_edge_count -= 1
    graph.observations -= support
    if graph.end_window:
        graph.terminal_observations -= (
            support - graph.internal_support[edge_id]
        )
    return support


def _spell_path(graph: BidirectedDeBruijnGraph, path: list[Handle]) -> str:
    sequence = oriented_sequence(graph, path[0])
    return sequence + "".join(
        _oriented_last_base(graph, handle) for handle in path[1:]
    )


def spell_edge_path(
    graph: BidirectedDeBruijnGraph,
    start: Handle,
    edge_ids: Sequence[int],
) -> str:
    """Spell an oriented sequence from a start handle and physical edges."""
    sequence = [oriented_sequence(graph, start)]
    current = start
    for edge_id in edge_ids:
        current = _successor(graph, current, edge_id)
        sequence.append(_oriented_last_base(graph, current))
    return "".join(sequence)


def iter_bidirected_unitig_paths(
    graph: BidirectedDeBruijnGraph,
) -> Iterator[tuple[Handle, tuple[int, ...]]]:
    """Yield each maximal physical path as a start handle and edge IDs."""
    visited = bytearray(graph.edge_count)
    handle_count = graph.node_count * 2

    for start in range(handle_count):
        in_degree = graph.out_degrees[flip_handle(start)]
        out_degree = graph.out_degrees[start]
        if in_degree == 1 and out_degree == 1:
            continue

        for edge_id in _outgoing_edges(graph, start):
            if visited[edge_id]:
                continue

            path = [edge_id]
            visited[edge_id] = 1
            current = _successor(graph, start, edge_id)

            while (
                graph.out_degrees[flip_handle(current)] == 1
                and graph.out_degrees[current] == 1
            ):
                following_edge = graph.out_edges[current * _MAX_OUT_DEGREE]
                if visited[following_edge]:
                    break
                current = _successor(graph, current, following_edge)
                path.append(following_edge)
                visited[following_edge] = 1

            yield start, tuple(path)

    for start in range(handle_count):
        for edge_id in _outgoing_edges(graph, start):
            if visited[edge_id]:
                continue

            path = [edge_id]
            visited[edge_id] = 1
            current = _successor(graph, start, edge_id)

            while current != start:
                following_edge = graph.out_edges[current * _MAX_OUT_DEGREE]
                if visited[following_edge]:
                    break
                current = _successor(graph, current, following_edge)
                path.append(following_edge)
                visited[following_edge] = 1

            yield start, tuple(path)


def extract_bidirected_unitigs(graph: BidirectedDeBruijnGraph) -> list[str]:
    """Return one sequence for each maximal non-branching physical path."""
    return [
        spell_edge_path(graph, start, edge_ids)
        for start, edge_ids in iter_bidirected_unitig_paths(graph)
    ]


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

    active_nodes = sum(
        graph.out_degrees[node_id << 1] > 0
        or graph.out_degrees[(node_id << 1) | 1] > 0
        for node_id in range(graph.node_count)
    )

    return GraphSummary(
        nodes=active_nodes,
        edges=graph.active_edge_count,
        observations=graph.observations,
        sources=sum(source_nodes),
        sinks=sum(sink_nodes),
        branching_nodes=sum(branching_nodes),
        unitigs=len(unitigs),
    )
