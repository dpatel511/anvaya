"""Conservative simplification of compact orientation-aware graphs."""

from dataclasses import dataclass

from anvaya.bidirected import (
    BidirectedDeBruijnGraph,
    _MAX_OUT_DEGREE,
    _deactivate_edge,
    _edge_support_total,
    _outgoing_edges,
    _successor,
    flip_handle,
)

_SUPPORT_RATIO_DENOMINATOR = 5
_MAX_CLEANING_ROUNDS = 5


@dataclass(slots=True, frozen=True)
class TipCleaningSummary:
    """Counts describing one conservative tip-cleaning run."""

    tips_removed: int = 0
    edges_removed: int = 0
    observations_removed: int = 0
    rounds: int = 0


@dataclass(slots=True, frozen=True)
class TipCandidate:
    """A weak dead-end path reported before graph modification."""

    start: int
    end: int
    edge_ids: tuple[int, ...]


def _tip_candidate(
    graph: BidirectedDeBruijnGraph,
    start: int,
    max_tip_length: int,
) -> TipCandidate | None:
    path: list[int] = []
    current = start
    maximum_support = 0

    while graph.out_degrees[current] == 1:
        edge_id = graph.out_edges[current * _MAX_OUT_DEGREE]
        path.append(edge_id)
        if graph.node_length + len(path) > max_tip_length:
            return None

        maximum_support = max(maximum_support, _edge_support_total(graph, edge_id))
        current = _successor(graph, current, edge_id)
        if not (
            graph.out_degrees[flip_handle(current)] == 1
            and graph.out_degrees[current] == 1
        ):
            break

    if graph.out_degrees[flip_handle(current)] <= 1:
        return None

    last_edge = path[-1]
    competing_support = max(
        (
            _edge_support_total(graph, edge_id)
            for edge_id in _outgoing_edges(graph, flip_handle(current))
            if edge_id != last_edge
        ),
        default=0,
    )
    if maximum_support * _SUPPORT_RATIO_DENOMINATOR > competing_support:
        return None
    return TipCandidate(start=start, end=current, edge_ids=tuple(path))


def find_weak_tip_candidates(
    graph: BidirectedDeBruijnGraph,
    max_tip_length: int | None = None,
) -> list[TipCandidate]:
    """Report weak tips without changing the graph."""
    if max_tip_length is None:
        max_tip_length = 2 * (graph.node_length + 1)
    if max_tip_length < graph.node_length + 1:
        raise ValueError("max_tip_length must be at least k")

    candidates: list[TipCandidate] = []
    seen: set[tuple[int, ...]] = set()
    for start in range(graph.node_count * 2):
        if (
            graph.out_degrees[flip_handle(start)] != 0
            or graph.out_degrees[start] != 1
        ):
            continue
        candidate = _tip_candidate(graph, start, max_tip_length)
        if candidate is None:
            continue
        signature = tuple(sorted(candidate.edge_ids))
        if signature in seen:
            continue
        seen.add(signature)
        candidates.append(candidate)
    return candidates


def remove_weak_tips(
    graph: BidirectedDeBruijnGraph,
    max_tip_length: int | None = None,
) -> TipCleaningSummary:
    """Remove short dead ends supported at no more than 20% of a competitor."""
    if max_tip_length is None:
        max_tip_length = 2 * (graph.node_length + 1)
    if max_tip_length < graph.node_length + 1:
        raise ValueError("max_tip_length must be at least k")

    tips_removed = 0
    edges_removed = 0
    observations_removed = 0
    rounds = 0

    for _ in range(_MAX_CLEANING_ROUNDS):
        removed_this_round = 0
        for start in range(graph.node_count * 2):
            if (
                graph.out_degrees[flip_handle(start)] != 0
                or graph.out_degrees[start] != 1
            ):
                continue

            candidate = _tip_candidate(graph, start, max_tip_length)
            if candidate is None:
                continue

            for edge_id in candidate.edge_ids:
                observations_removed += _deactivate_edge(graph, edge_id)
            tips_removed += 1
            edges_removed += len(candidate.edge_ids)
            removed_this_round += 1

        if removed_this_round == 0:
            break
        rounds += 1

    return TipCleaningSummary(
        tips_removed=tips_removed,
        edges_removed=edges_removed,
        observations_removed=observations_removed,
        rounds=rounds,
    )
