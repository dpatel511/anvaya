"""Basic de Bruijn graph construction and topology."""

from collections import Counter
from collections.abc import Sequence

from anvaya.kmers import kmers

DeBruijnGraph = dict[str, Counter[str]]


def build_dbg(reads: Sequence[str], k: int) -> DeBruijnGraph:
    """Build a directed de Bruijn graph with edge multiplicities.

    Nodes are (k-1)-mers. Each observed k-mer contributes an edge from its
    prefix to its suffix. Repeated observations increase the edge count.
    """
    if not reads:
        raise ValueError("reads must not be empty")
    if not isinstance(k, int) or isinstance(k, bool):
        raise TypeError("k must be an integer")
    if k < 2:
        raise ValueError("k must be at least 2 for graph construction")

    graph: DeBruijnGraph = {}

    for read in reads:
        if len(read) < k:
            raise ValueError("each read must be at least k bases long")

        for kmer in kmers(read, k):
            prefix = kmer[:-1]
            suffix = kmer[1:]
            graph.setdefault(prefix, Counter())[suffix] += 1
            graph.setdefault(suffix, Counter())

    return graph


def out_degree(graph: DeBruijnGraph, node: str) -> int:
    """Return the number of distinct outgoing edges from *node*."""
    return len(graph[node])


def in_degree(graph: DeBruijnGraph, node: str) -> int:
    """Return the number of distinct incoming edges to *node*."""
    if node not in graph:
        raise KeyError(node)
    return sum(node in successors for successors in graph.values())


def source_nodes(graph: DeBruijnGraph) -> set[str]:
    """Return nodes with no incoming edges and at least one outgoing edge."""
    return {
        node
        for node in graph
        if in_degree(graph, node) == 0 and out_degree(graph, node) > 0
    }


def sink_nodes(graph: DeBruijnGraph) -> set[str]:
    """Return nodes with at least one incoming edge and no outgoing edges."""
    return {
        node
        for node in graph
        if in_degree(graph, node) > 0 and out_degree(graph, node) == 0
    }


def branching_nodes(graph: DeBruijnGraph) -> set[str]:
    """Return nodes with multiple incoming or outgoing edges."""
    return {
        node
        for node in graph
        if in_degree(graph, node) > 1 or out_degree(graph, node) > 1
    }
