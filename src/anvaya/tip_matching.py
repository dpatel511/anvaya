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
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    current = start
    for index, edge_id in enumerate(edge_ids):
        evidence = graph.end_support(edge_id)
        left = evidence.left
        right = evidence.right
        if current != graph.edge_sources[edge_id]:
            left, right = right, left
        if index == edge_index:
            return left, right
        current = _successor(graph, current, edge_id)
    return (), ()


def _exact_end_count(
    graph: BidirectedDeBruijnGraph,
    counts: tuple[int, ...],
    expected_end: str,
) -> int:
    """Count observations where the appended base is inside the end window."""
    offset = graph.node_length if expected_end == "left" else 0
    return sum(
        count
        for distance, count in enumerate(counts)
        if distance + offset < graph.end_window
    )


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
        if (
            change.reference == "C"
            and change.alternative == "T"
            and _exact_end_count(graph, left, "left")
        ):
            continue
        if (
            change.reference == "G"
            and change.alternative == "A"
            and _exact_end_count(graph, right, "right")
        ):
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
    return match_branch_to_backbone(graph, branch_handle, tip_edge_ids)


def match_branch_to_backbone(
    graph: BidirectedDeBruijnGraph,
    branch_handle: int,
    branch_edge_ids: tuple[int, ...],
) -> TipBackboneMatch | None:
    """Match one oriented branch to its best equal-length linear competitor."""
    if not branch_edge_ids:
        return None

    tip_edge_ids = branch_edge_ids
    excluded = set(tip_edge_ids)
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

        candidate = match_competing_paths(
            graph,
            branch_handle,
            tip_edge_ids,
            backbone_edge_ids,
            require_same_end=False,
        )
        assert candidate is not None
        key = (
            candidate.ry_identity,
            candidate.sequence_identity,
            candidate.backbone_minimum_edge_support,
            candidate.backbone_observations,
        )
        if best_key is None or key > best_key:
            best = candidate
            best_key = key

    return best


def match_competing_paths(
    graph: BidirectedDeBruijnGraph,
    branch_handle: int,
    alternative_edge_ids: tuple[int, ...],
    reference_edge_ids: tuple[int, ...],
    *,
    require_same_end: bool = True,
) -> TipBackboneMatch | None:
    """Compare two explicit equal-length paths leaving the same handle."""
    if (
        not alternative_edge_ids
        or len(alternative_edge_ids) != len(reference_edge_ids)
    ):
        return None

    def path_end(edge_ids: tuple[int, ...]) -> int:
        current = branch_handle
        for edge_id in edge_ids:
            current = _successor(graph, current, edge_id)
        return current

    if (
        require_same_end
        and path_end(alternative_edge_ids) != path_end(reference_edge_ids)
    ):
        return None

    alternative_sequence = spell_edge_path(
        graph,
        branch_handle,
        alternative_edge_ids,
    )
    reference_sequence = spell_edge_path(
        graph,
        branch_handle,
        reference_edge_ids,
    )
    if len(alternative_sequence) != len(reference_sequence):
        return None

    alternative_observations = sum(
        _edge_support_total(graph, edge_id) for edge_id in alternative_edge_ids
    )
    reference_supports = tuple(
        _edge_support_total(graph, edge_id) for edge_id in reference_edge_ids
    )
    reference_observations = sum(reference_supports)
    changes = substitutions(reference_sequence, alternative_sequence)
    return TipBackboneMatch(
        branch_handle=branch_handle,
        tip_edge_ids=alternative_edge_ids,
        backbone_edge_ids=reference_edge_ids,
        tip_sequence=alternative_sequence,
        backbone_sequence=reference_sequence,
        tip_observations=alternative_observations,
        backbone_observations=reference_observations,
        backbone_minimum_edge_support=min(reference_supports),
        relative_coverage=(
            alternative_observations / reference_observations
            if reference_observations
            else 0.0
        ),
        sequence_identity=_identity(reference_sequence, alternative_sequence),
        ry_identity=_identity(
            _ry_sequence(reference_sequence),
            _ry_sequence(alternative_sequence),
        ),
        substitutions=changes,
        damage_compatible=_damage_compatible(
            graph,
            branch_handle,
            alternative_edge_ids,
            changes,
        ),
    )
