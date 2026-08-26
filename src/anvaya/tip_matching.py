"""Non-destructive matching of weak tips to local backbone paths."""

from dataclasses import dataclass

from anvaya.bidirected import (
    BidirectedDeBruijnGraph,
    _edge_support_total,
    _outgoing_edges,
    _successor,
    flip_handle,
    spell_edge_path,
)
from anvaya.classification import Substitution, substitutions
from anvaya.cleaning import TipCandidate

_RY_TRANSLATION = str.maketrans({"A": "R", "G": "R", "C": "Y", "T": "Y"})


@dataclass(slots=True, frozen=True)
class TipBackboneMatch:
    """A weak tip aligned to one locally competing graph continuation."""

    branch_handle: int
    tip_edge_ids: tuple[int, ...]
    backbone_edge_ids: tuple[int, ...]
    tip_sequence: str
    backbone_sequence: str
    tip_observations: int
    backbone_observations: int
    backbone_minimum_edge_support: int
    relative_coverage: float
    sequence_identity: float
    ry_identity: float
    substitutions: tuple[Substitution, ...]
    damage_compatible: bool


def _identity(left: str, right: str) -> float:
    if len(left) != len(right) or not left:
        return 0.0
    return sum(
        left_base == right_base
        for left_base, right_base in zip(left, right, strict=True)
    ) / len(left)


def _ry_sequence(sequence: str) -> str:
    return sequence.translate(_RY_TRANSLATION)


def _linear_continuation(
    graph: BidirectedDeBruijnGraph,
    start: int,
    first_edge: int,
    edge_count: int,
    excluded_edges: set[int],
) -> tuple[int, ...] | None:
    path = [first_edge]
    current = _successor(graph, start, first_edge)

    while len(path) < edge_count:
        if not (
            graph.out_degrees[current] == 1
            and graph.out_degrees[flip_handle(current)] == 1
        ):
            return None
        following = next(_outgoing_edges(graph, current))
        if following in excluded_edges or following in path:
            return None
        path.append(following)
        current = _successor(graph, current, following)

    return tuple(path)


def _oriented_terminal_counts(
    graph: BidirectedDeBruijnGraph,
    start: int,
    edge_ids: tuple[int, ...],
    edge_index: int,
) -> tuple[int, int]:
    current = start
    for index, edge_id in enumerate(edge_ids):
        evidence = graph.end_support(edge_id)
        left = sum(evidence.left)
        right = sum(evidence.right)
        if current != graph.edge_sources[edge_id]:
            left, right = right, left
        if index == edge_index:
            return left, right
        current = _successor(graph, current, edge_id)
    return 0, 0


def _damage_compatible(
    graph: BidirectedDeBruijnGraph,
    branch_handle: int,
    tip_edge_ids: tuple[int, ...],
    changes: tuple[Substitution, ...],
) -> bool:
    if graph.end_window == 0 or not changes:
        return False
    for change in changes:
        edge_index = change.position - graph.node_length - 1
        if not 0 <= edge_index < len(tip_edge_ids):
            return False
        left, right = _oriented_terminal_counts(
            graph,
            branch_handle,
            tip_edge_ids,
            edge_index,
        )
        if change.reference == "C" and change.alternative == "T" and left:
            continue
        if change.reference == "G" and change.alternative == "A" and right:
            continue
        return False
    return True


def match_tip_to_backbone(
    graph: BidirectedDeBruijnGraph,
    tip: TipCandidate,
) -> TipBackboneMatch | None:
    """Match a weak tip to its best equal-length linear competitor."""
    if not tip.edge_ids:
        return None

    branch_handle = flip_handle(tip.end)
    tip_edge_ids = tuple(reversed(tip.edge_ids))
    tip_sequence = spell_edge_path(graph, branch_handle, tip_edge_ids)
    tip_observations = sum(
        _edge_support_total(graph, edge_id) for edge_id in tip_edge_ids
    )
    excluded = set(tip.edge_ids)
    best: TipBackboneMatch | None = None
    best_key: tuple[float, float, int, int] | None = None

    for edge_id in sorted(_outgoing_edges(graph, branch_handle)):
        if edge_id in excluded:
            continue
        backbone_edge_ids = _linear_continuation(
            graph,
            branch_handle,
            edge_id,
            len(tip_edge_ids),
            excluded,
        )
        if backbone_edge_ids is None:
            continue

        backbone_sequence = spell_edge_path(
            graph,
            branch_handle,
            backbone_edge_ids,
        )
        sequence_identity = _identity(backbone_sequence, tip_sequence)
        ry_identity = _identity(
            _ry_sequence(backbone_sequence),
            _ry_sequence(tip_sequence),
        )
        backbone_supports = tuple(
            _edge_support_total(graph, path_edge)
            for path_edge in backbone_edge_ids
        )
        backbone_observations = sum(backbone_supports)
        minimum_support = min(backbone_supports)
        changes = substitutions(backbone_sequence, tip_sequence)
        candidate = TipBackboneMatch(
            branch_handle=branch_handle,
            tip_edge_ids=tip_edge_ids,
            backbone_edge_ids=backbone_edge_ids,
            tip_sequence=tip_sequence,
            backbone_sequence=backbone_sequence,
            tip_observations=tip_observations,
            backbone_observations=backbone_observations,
            backbone_minimum_edge_support=minimum_support,
            relative_coverage=(
                tip_observations / backbone_observations
                if backbone_observations
                else 0.0
            ),
            sequence_identity=sequence_identity,
            ry_identity=ry_identity,
            substitutions=changes,
            damage_compatible=_damage_compatible(
                graph,
                branch_handle,
                tip_edge_ids,
                changes,
            ),
        )
        key = (
            ry_identity,
            sequence_identity,
            minimum_support,
            backbone_observations,
        )
        if best_key is None or key > best_key:
            best = candidate
            best_key = key

    return best
