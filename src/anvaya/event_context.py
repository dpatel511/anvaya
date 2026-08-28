"""Local graph context for calibrating matched event likelihoods."""

from dataclasses import dataclass

from anvaya.bidirected import BidirectedDeBruijnGraph, _edge_support_total
from anvaya.tip_classification import classify_tip_match
from anvaya.tip_matching import TipBackboneMatch


@dataclass(slots=True, frozen=True)
class EventGraphContext:
    """Report-only topology, multiplicity, and molecule evidence."""

    event_type: str
    path_edge_count: int
    substitution_count: int
    alternative_path_observations: int
    reference_path_observations: int
    alternative_minimum_edge_support: int
    reference_minimum_edge_support: int
    relative_coverage: float
    sequence_identity: float
    ry_identity: float
    damage_compatible: bool
    branch_out_degree: int
    terminal_fraction: float
    strand_balance: float | None
    terminal_enrichment: float
    substitution_terminal_fraction: float
    mean_damage_distance: float | None
    joint_molecule_fraction: float | None
    mean_base_quality: float | None


def extract_event_graph_context(
    graph: BidirectedDeBruijnGraph,
    match: TipBackboneMatch,
    event_type: str,
) -> EventGraphContext:
    """Extract existing evidence without modifying or reclassifying the graph."""
    evidence = classify_tip_match(graph, match).evidence
    return EventGraphContext(
        event_type=event_type,
        path_edge_count=len(match.tip_edge_ids),
        substitution_count=len(match.substitutions),
        alternative_path_observations=match.tip_observations,
        reference_path_observations=match.backbone_observations,
        alternative_minimum_edge_support=min(
            _edge_support_total(graph, edge_id)
            for edge_id in match.tip_edge_ids
        ),
        reference_minimum_edge_support=match.backbone_minimum_edge_support,
        relative_coverage=match.relative_coverage,
        sequence_identity=match.sequence_identity,
        ry_identity=match.ry_identity,
        damage_compatible=match.damage_compatible,
        branch_out_degree=graph.out_degrees[match.branch_handle],
        terminal_fraction=evidence.tip_terminal_fraction,
        strand_balance=evidence.tip_strand_balance,
        terminal_enrichment=evidence.terminal_enrichment,
        substitution_terminal_fraction=(
            evidence.substitution_terminal_fraction
        ),
        mean_damage_distance=evidence.mean_damage_distance,
        joint_molecule_fraction=evidence.joint_molecule_fraction,
        mean_base_quality=evidence.mean_base_quality,
    )
