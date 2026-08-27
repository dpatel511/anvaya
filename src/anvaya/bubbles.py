"""Non-destructive detection of simple bubbles in bidirected graphs."""

from dataclasses import dataclass

from anvaya.bidirected import (
    BidirectedDeBruijnGraph,
    _edge_support_total,
    _outgoing_edges,
    _successor,
    flip_handle,
)


@dataclass(slots=True, frozen=True)
class BubblePath:
    """One alternative path through a bubble."""

    edge_ids: tuple[int, ...]
    observations: int
    minimum_edge_support: int


@dataclass(slots=True, frozen=True)
class Bubble:
    """Alternative paths that leave and rejoin the same oriented handles."""

    start: int
    end: int
    paths: tuple[BubblePath, ...]


def trace_linear_branch(
    graph: BidirectedDeBruijnGraph,
    start: int,
    first_edge: int,
    max_edges: int,
) -> tuple[int, tuple[int, ...]] | None:
    """Follow one branch until its first non-linear handle."""
    path = [first_edge]
    current = _successor(graph, start, first_edge)
    visited = {start}

    while (
        current not in visited
        and graph.out_degrees[flip_handle(current)] == 1
        and graph.out_degrees[current] == 1
    ):
        if len(path) == max_edges:
            return None
        visited.add(current)
        edge_id = next(_outgoing_edges(graph, current))
        path.append(edge_id)
        current = _successor(graph, current, edge_id)

    if current in visited:
        return None
    return current, tuple(path)


def find_simple_bubbles(
    graph: BidirectedDeBruijnGraph,
    max_edges: int = 64,
) -> list[Bubble]:
    """Report bounded, edge-disjoint linear alternatives without graph changes."""
    if not isinstance(max_edges, int) or isinstance(max_edges, bool):
        raise TypeError("max_edges must be an integer")
    if max_edges < 1:
        raise ValueError("max_edges must be at least 1")

    bubbles: list[Bubble] = []
    seen: set[tuple[tuple[int, ...], ...]] = set()

    for start in range(graph.node_count * 2):
        if graph.out_degrees[start] < 2:
            continue

        branches = [
            trace_linear_branch(graph, start, edge_id, max_edges)
            for edge_id in sorted(_outgoing_edges(graph, start))
        ]
        if any(branch is None for branch in branches):
            continue

        completed = [branch for branch in branches if branch is not None]
        ends = {end for end, _ in completed}
        if len(ends) != 1:
            continue
        end = completed[0][0]
        if end == start or graph.out_degrees[flip_handle(end)] < 2:
            continue

        edge_paths = [edge_ids for _, edge_ids in completed]
        used_edges: set[int] = set()
        paths_are_disjoint = True
        for edge_ids in edge_paths:
            if used_edges.intersection(edge_ids):
                paths_are_disjoint = False
                break
            used_edges.update(edge_ids)
        if not paths_are_disjoint:
            continue

        signature = tuple(sorted(tuple(sorted(path)) for path in edge_paths))
        if signature in seen:
            continue
        seen.add(signature)

        paths = tuple(
            BubblePath(
                edge_ids=edge_ids,
                observations=sum(
                    _edge_support_total(graph, edge_id) for edge_id in edge_ids
                ),
                minimum_edge_support=min(
                    _edge_support_total(graph, edge_id) for edge_id in edge_ids
                ),
            )
            for edge_ids in edge_paths
        )
        bubbles.append(Bubble(start=start, end=end, paths=paths))

    return bubbles
