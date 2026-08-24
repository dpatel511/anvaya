"""Basic de Bruijn graph construction."""

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
