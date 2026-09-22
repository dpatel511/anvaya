"""Conservative opt-in substitution guard; no haplotype or strain inference."""

from collections import Counter

from anvaya.mixture_diagnostics import robust_alleles


def variable_neighbours(columns):
    """Select neighbours independently of the focal base, from immutable columns."""
    neighbours = {}
    for position, column in enumerate(columns):
        robust = robust_alleles(column)
        if sum(count >= 2 for count in Counter(robust.values()).values()) >= 2:
            neighbours[position] = robust
    return neighbours


def blocking_neighbour(position, column, before, proposed, neighbours):
    """Return a 0-based witness or None, requiring a complete two-group partition.

    All Q20 focal molecules must have robust neighbour observations. Their joint
    table must contain exactly two cells, different at both sites, each supported
    by >=2 molecules. Each focal allele needs >=1 robust observation as well.
    No mixture-score threshold, statistical test or read reassignment is applied.
    """
    if before == proposed:
        return None
    focal = {m: o[0] for m, o in column.items() if o is not None and o[1] >= 20}
    if set(focal.values()) != {before, proposed}:
        return None
    counts = Counter(focal.values())
    if min(counts.values()) < 2:
        return None
    if not {before, proposed} <= set(robust_alleles(column).values()):
        return None
    for neighbour, robust in neighbours.items():
        if neighbour == position or not focal.keys() <= robust.keys():
            continue
        joint = Counter((a, robust[m]) for m, a in focal.items())
        if len(joint) == 2 and len({b for _, b in joint}) == 2 and min(joint.values()) >= 2:
            return neighbour
    return None
