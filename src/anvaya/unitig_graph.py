"""Compacted unitig-level representation of an orientation-aware graph."""

from array import array
from dataclasses import dataclass, field

from anvaya.bidirected import (
    BidirectedDeBruijnGraph,
    Handle,
    _edge_support_total,
    _oriented_last_base,
    _outgoing_edges,
    _successor,
    flip_handle,
    iter_bidirected_unitig_paths,
    oriented_sequence,
)
from anvaya.sequences import reverse_complement

_MAX_LINK_DEGREE = 4


@dataclass(slots=True, frozen=True)
class UnitigSupport:
    """Coverage and read-end evidence aggregated across one unitig."""

    edge_count: int
    observations: int
    minimum_edge_support: int
    mean_edge_support: float
    forward_observations: int
    reverse_observations: int
    palindromic_observations: int
    left_terminal_observations: int
    right_terminal_observations: int
    ambiguous_terminal_observations: int
    terminal_observations: int
    internal_observations: int


@dataclass(slots=True)
class CompactedUnitigGraph:
    """Physical unitigs with links addressable through oriented handles."""

    k: int
    sequences: list[str] = field(default_factory=list)
    start_handles: array = field(default_factory=lambda: array("Q"))
    end_handles: array = field(default_factory=lambda: array("Q"))
    first_edge_ids: array = field(default_factory=lambda: array("I"))
    last_edge_ids: array = field(default_factory=lambda: array("I"))
    edge_counts: array = field(default_factory=lambda: array("I"))
    observations: array = field(default_factory=lambda: array("Q"))
    minimum_edge_support: array = field(default_factory=lambda: array("I"))
    forward_observations: array = field(default_factory=lambda: array("Q"))
    reverse_observations: array = field(default_factory=lambda: array("Q"))
    palindromic_observations: array = field(default_factory=lambda: array("Q"))
    left_terminal_observations: array = field(default_factory=lambda: array("Q"))
    right_terminal_observations: array = field(default_factory=lambda: array("Q"))
    ambiguous_terminal_observations: array = field(
        default_factory=lambda: array("Q")
    )
    terminal_observations: array = field(default_factory=lambda: array("Q"))
    internal_observations: array = field(default_factory=lambda: array("Q"))
    edge_handles: array = field(default_factory=lambda: array("q"))
    out_links: array = field(default_factory=lambda: array("q"))
    out_degrees: bytearray = field(default_factory=bytearray)

    @property
    def unitig_count(self) -> int:
        return len(self.sequences)

    @property
    def oriented_link_count(self) -> int:
        return sum(self.out_degrees)

    def oriented_sequence(self, handle: Handle) -> str:
        """Return a unitig sequence in the requested orientation."""
        sequence = self.sequences[handle >> 1]
        return reverse_complement(sequence) if handle & 1 else sequence

    def outgoing(self, handle: Handle) -> tuple[Handle, ...]:
        """Return oriented unitigs reachable from one oriented unitig."""
        offset = handle * _MAX_LINK_DEGREE
        return tuple(
            self.out_links[offset + index]
            for index in range(self.out_degrees[handle])
        )

    def support(self, unitig_id: int) -> UnitigSupport:
        """Return aggregate support for one physical unitig."""
        edge_count = self.edge_counts[unitig_id]
        observations = self.observations[unitig_id]
        return UnitigSupport(
            edge_count=edge_count,
            observations=observations,
            minimum_edge_support=self.minimum_edge_support[unitig_id],
            mean_edge_support=observations / edge_count,
            forward_observations=self.forward_observations[unitig_id],
            reverse_observations=self.reverse_observations[unitig_id],
            palindromic_observations=self.palindromic_observations[unitig_id],
            left_terminal_observations=(
                self.left_terminal_observations[unitig_id]
            ),
            right_terminal_observations=(
                self.right_terminal_observations[unitig_id]
            ),
            ambiguous_terminal_observations=(
                self.ambiguous_terminal_observations[unitig_id]
            ),
            terminal_observations=self.terminal_observations[unitig_id],
            internal_observations=self.internal_observations[unitig_id],
        )

    def oriented_support(self, handle: Handle) -> UnitigSupport:
        """Return evidence oriented consistently with a unitig handle."""
        support = self.support(handle >> 1)
        if not handle & 1:
            return support
        return UnitigSupport(
            edge_count=support.edge_count,
            observations=support.observations,
            minimum_edge_support=support.minimum_edge_support,
            mean_edge_support=support.mean_edge_support,
            forward_observations=support.reverse_observations,
            reverse_observations=support.forward_observations,
            palindromic_observations=support.palindromic_observations,
            left_terminal_observations=support.right_terminal_observations,
            right_terminal_observations=support.left_terminal_observations,
            ambiguous_terminal_observations=(
                support.ambiguous_terminal_observations
            ),
            terminal_observations=support.terminal_observations,
            internal_observations=support.internal_observations,
        )


def _append_link(
    graph: CompactedUnitigGraph, source: Handle, target: Handle
) -> None:
    degree = graph.out_degrees[source]
    if degree >= _MAX_LINK_DEGREE:
        raise ValueError("oriented unitig has more than four outgoing links")
    graph.out_links[source * _MAX_LINK_DEGREE + degree] = target
    graph.out_degrees[source] = degree + 1


def build_compacted_unitig_graph(
    source_graph: BidirectedDeBruijnGraph,
    *,
    track_edge_handles: bool = False,
) -> CompactedUnitigGraph:
    """Compact active non-branching paths while retaining links and evidence."""
    if not isinstance(track_edge_handles, bool):
        raise TypeError("track_edge_handles must be a boolean")
    graph = CompactedUnitigGraph(k=source_graph.node_length + 1)
    if track_edge_handles:
        graph.edge_handles = array("q", [-1]) * source_graph.edge_count
    entries: dict[tuple[Handle, int], Handle] = {}
    internal_exits: set[tuple[Handle, int]] = set()

    for start, edge_ids in iter_bidirected_unitig_paths(source_graph):
        unitig_id = graph.unitig_count
        current = start
        observations = 0
        minimum_support: int | None = None
        forward = 0
        reverse = 0
        palindromic = 0
        left_terminal = 0
        right_terminal = 0
        ambiguous_terminal = 0
        terminal = 0
        internal = 0
        sequence = [oriented_sequence(source_graph, start)]

        for edge_id in edge_ids:
            traverses_forward = current == source_graph.edge_sources[edge_id]
            if track_edge_handles:
                graph.edge_handles[edge_id] = (
                    unitig_id << 1
                    if traverses_forward
                    else (unitig_id << 1) | 1
                )
            support = _edge_support_total(source_graph, edge_id)
            observations += support
            minimum_support = (
                support if minimum_support is None else min(minimum_support, support)
            )
            edge_forward = source_graph.forward_support[edge_id]
            edge_reverse = source_graph.reverse_support[edge_id]
            if not traverses_forward:
                edge_forward, edge_reverse = edge_reverse, edge_forward
            forward += edge_forward
            reverse += edge_reverse
            palindromic += source_graph.palindromic_support[edge_id]
            if source_graph.end_window:
                edge_internal = source_graph.internal_support[edge_id]
                internal += edge_internal
                terminal += support - edge_internal
                offset = edge_id * source_graph.end_window
                end = offset + source_graph.end_window
                edge_left = sum(source_graph.left_end_support[offset:end])
                edge_right = sum(source_graph.right_end_support[offset:end])
                if not traverses_forward:
                    edge_left, edge_right = edge_right, edge_left
                left_terminal += edge_left
                right_terminal += edge_right
                ambiguous_terminal += source_graph.ambiguous_end_support[edge_id]
            else:
                internal += support
            current = _successor(source_graph, current, edge_id)
            sequence.append(_oriented_last_base(source_graph, current))

        graph.sequences.append("".join(sequence))
        graph.start_handles.append(start)
        graph.end_handles.append(current)
        graph.first_edge_ids.append(edge_ids[0])
        graph.last_edge_ids.append(edge_ids[-1])
        graph.edge_counts.append(len(edge_ids))
        graph.observations.append(observations)
        graph.minimum_edge_support.append(minimum_support or 0)
        graph.forward_observations.append(forward)
        graph.reverse_observations.append(reverse)
        graph.palindromic_observations.append(palindromic)
        graph.left_terminal_observations.append(left_terminal)
        graph.right_terminal_observations.append(right_terminal)
        graph.ambiguous_terminal_observations.append(ambiguous_terminal)
        graph.terminal_observations.append(terminal)
        graph.internal_observations.append(internal)

        entries[(start, edge_ids[0])] = unitig_id << 1
        reverse_entry = (flip_handle(current), edge_ids[-1])
        entries.setdefault(reverse_entry, (unitig_id << 1) | 1)
        path_edges = set(edge_ids)
        for edge_id in _outgoing_edges(source_graph, current):
            if edge_id in path_edges:
                internal_exits.add((unitig_id << 1, edge_id))
        reverse_end = flip_handle(start)
        for edge_id in _outgoing_edges(source_graph, reverse_end):
            if edge_id in path_edges:
                internal_exits.add(((unitig_id << 1) | 1, edge_id))

    handle_count = graph.unitig_count * 2
    graph.out_degrees = bytearray(handle_count)
    graph.out_links = array("q", [-1]) * (handle_count * _MAX_LINK_DEGREE)

    for unitig_id in range(graph.unitig_count):
        for orientation in (0, 1):
            source = (unitig_id << 1) | orientation
            end = (
                graph.end_handles[unitig_id]
                if orientation == 0
                else flip_handle(graph.start_handles[unitig_id])
            )
            for edge_id in _outgoing_edges(source_graph, end):
                target = entries.get((end, edge_id))
                if target is None:
                    if (source, edge_id) in internal_exits:
                        continue
                    raise RuntimeError("unitig entry is missing for an active graph edge")
                _append_link(graph, source, target)

    return graph
