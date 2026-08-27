"""Non-destructive detection and matching of incomplete graph branches."""

from dataclasses import dataclass

from anvaya.bidirected import (
    BidirectedDeBruijnGraph,
    _edge_support_total,
    _outgoing_edges,
)
from anvaya.bubbles import trace_linear_branch
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.tip_matching import TipBackboneMatch, match_branch_to_backbone


@dataclass(slots=True, frozen=True)
class IncompleteBranchCandidate:
    """A weak bounded branch not already represented as a tip or bubble."""

    start: int
    end: int
    edge_ids: tuple[int, ...]
    observations: int
    minimum_edge_support: int
    competing_minimum_edge_support: int


def _path_support(
    graph: BidirectedDeBruijnGraph,
    edge_ids: tuple[int, ...],
) -> tuple[int, int, int]:
    supports = tuple(
        _edge_support_total(graph, edge_id) for edge_id in edge_ids
    )
    return sum(supports), min(supports), max(supports)


def find_incomplete_branch_candidates(
    graph: BidirectedDeBruijnGraph,
    max_edges: int = 64,
    maximum_relative_support: float = 0.50,
) -> list[IncompleteBranchCandidate]:
    """Report weak non-rejoining branches bounded by graph decisions."""
    if not isinstance(max_edges, int) or isinstance(max_edges, bool):
        raise TypeError("max_edges must be an integer")
    if max_edges < 1:
        raise ValueError("max_edges must be at least 1")
    if (
        not isinstance(maximum_relative_support, (int, float))
        or isinstance(maximum_relative_support, bool)
    ):
        raise TypeError("maximum_relative_support must be numeric")
    if not 0.0 <= maximum_relative_support <= 1.0:
        raise ValueError("maximum_relative_support must be between 0 and 1")

    candidates: list[IncompleteBranchCandidate] = []
    seen: set[tuple[int, ...]] = set()
    tip_signatures = {
        tuple(sorted(tip.edge_ids))
        for tip in find_weak_tip_candidates(graph)
    }

    for start in range(graph.node_count * 2):
        if graph.out_degrees[start] < 2:
            continue

        traced = [
            trace_linear_branch(graph, start, edge_id, max_edges)
            for edge_id in sorted(_outgoing_edges(graph, start))
        ]
        if any(branch is None for branch in traced):
            continue
        branches = [branch for branch in traced if branch is not None]
        if len({end for end, _ in branches}) == 1:
            continue

        supported = [
            (end, edge_ids, *_path_support(graph, edge_ids))
            for end, edge_ids in branches
        ]
        strongest = max(
            supported,
            key=lambda branch: (branch[3], branch[2], len(branch[1])),
        )
        strongest_end, strongest_edges, _, strongest_minimum, _ = strongest

        for end, edge_ids, observations, minimum, maximum in supported:
            if edge_ids == strongest_edges or end == strongest_end:
                continue
            if maximum > strongest_minimum * maximum_relative_support:
                continue
            signature = tuple(sorted(edge_ids))
            if signature in tip_signatures:
                continue
            match = match_branch_to_backbone(graph, start, edge_ids)
            if match is None:
                continue
            if (
                maximum
                > match.backbone_minimum_edge_support
                * maximum_relative_support
            ):
                continue

            if signature in seen:
                continue
            seen.add(signature)
            candidates.append(
                IncompleteBranchCandidate(
                    start=start,
                    end=end,
                    edge_ids=edge_ids,
                    observations=observations,
                    minimum_edge_support=minimum,
                    competing_minimum_edge_support=(
                        match.backbone_minimum_edge_support
                    ),
                )
            )

    return candidates


def match_incomplete_branch_to_backbone(
    graph: BidirectedDeBruijnGraph,
    candidate: IncompleteBranchCandidate,
) -> TipBackboneMatch | None:
    """Match a detected incomplete branch to its local backbone."""
    return match_branch_to_backbone(
        graph,
        candidate.start,
        candidate.edge_ids,
    )
