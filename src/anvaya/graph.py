"""Basic de Bruijn graph construction and topology."""

from collections import Counter
from collections.abc import Sequence

from anvaya.kmers import kmers

DeBruijnGraph = dict[str, Counter[str]]
DegreeMap = dict[str, int]


def build_dbg(reads: Sequence[str], k: int) -> DeBruijnGraph:
    """Build a directed de Bruijn graph with edge multiplicities.

    Nodes are (k-1)-mers. Each unambiguous observed k-mer contributes an edge
    from its prefix to its suffix. Repeated observations increase the count.
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
            if "N" in kmer:
                continue
            prefix = kmer[:-1]
            suffix = kmer[1:]
            graph.setdefault(prefix, Counter())[suffix] += 1
            graph.setdefault(suffix, Counter())

    return graph


def out_degree(graph: DeBruijnGraph, node: str) -> int:
    """Return the number of distinct outgoing edges from *node*."""
    return len(graph[node])


def node_degrees(graph: DeBruijnGraph) -> tuple[DegreeMap, DegreeMap]:
    """Return distinct incoming and outgoing degrees for every node."""
    in_degrees = dict.fromkeys(graph, 0)
    out_degrees = {node: len(successors) for node, successors in graph.items()}

    for successors in graph.values():
        for successor in successors:
            in_degrees[successor] += 1

    return in_degrees, out_degrees


def in_degree(graph: DeBruijnGraph, node: str) -> int:
    """Return the number of distinct incoming edges to *node*."""
    if node not in graph:
        raise KeyError(node)
    return sum(node in successors for successors in graph.values())


def source_nodes(graph: DeBruijnGraph) -> set[str]:
    """Return nodes with no incoming edges and at least one outgoing edge."""
    in_degrees, out_degrees = node_degrees(graph)
    return {
        node
        for node in graph
        if in_degrees[node] == 0 and out_degrees[node] > 0
    }


def sink_nodes(graph: DeBruijnGraph) -> set[str]:
    """Return nodes with at least one incoming edge and no outgoing edges."""
    in_degrees, out_degrees = node_degrees(graph)
    return {
        node
        for node in graph
        if in_degrees[node] > 0 and out_degrees[node] == 0
    }


def branching_nodes(graph: DeBruijnGraph) -> set[str]:
    """Return nodes with multiple incoming or outgoing edges."""
    in_degrees, out_degrees = node_degrees(graph)
    return {
        node
        for node in graph
        if in_degrees[node] > 1 or out_degrees[node] > 1
    }
