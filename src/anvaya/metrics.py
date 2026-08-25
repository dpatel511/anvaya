"""Summary metrics for de Bruijn graphs."""

from dataclasses import dataclass

from anvaya.graph import branching_nodes, sink_nodes, source_nodes, DeBruijnGraph
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


def summarize_graph(graph: DeBruijnGraph) -> GraphSummary:
    """Return basic topology and support counts for *graph*."""
    return GraphSummary(
        nodes=len(graph),
        edges=sum(len(successors) for successors in graph.values()),
        observations=sum(sum(successors.values()) for successors in graph.values()),
        sources=len(source_nodes(graph)),
        sinks=len(sink_nodes(graph)),
        branching_nodes=len(branching_nodes(graph)),
        unitigs=len(extract_unitigs(graph)),
    )
