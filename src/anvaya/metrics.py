"""Graph-level assembly metrics."""

from dataclasses import dataclass

from anvaya.graph import DeBruijnGraph, node_degrees
from anvaya.unitigs import extract_unitigs


@dataclass(frozen=True)
class GraphSummary:
    """Counts describing the topology and support of a graph."""

    nodes: int
    edges: int
    observations: int
    sources: int
    sinks: int
    branching_nodes: int
    unitigs: int


def summarize_graph(
    graph: DeBruijnGraph, unitigs: list[str] | None = None
) -> GraphSummary:
    """Return basic topology and support counts for *graph*."""
    in_degrees, out_degrees = node_degrees(graph)
    if unitigs is None:
        unitigs = extract_unitigs(graph)

    return GraphSummary(
        nodes=len(graph),
        edges=sum(len(successors) for successors in graph.values()),
        observations=sum(sum(successors.values()) for successors in graph.values()),
        sources=sum(
            in_degrees[node] == 0 and out_degrees[node] > 0 for node in graph
        ),
        sinks=sum(
            in_degrees[node] > 0 and out_degrees[node] == 0 for node in graph
        ),
        branching_nodes=sum(
            in_degrees[node] > 1 or out_degrees[node] > 1 for node in graph
        ),
        unitigs=len(unitigs),
    )
