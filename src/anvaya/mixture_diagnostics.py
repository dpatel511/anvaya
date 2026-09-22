"""Diagnostic model evidence and allele linkage on fixed molecule observations.

No calls or placements are changed. Columns use the consensus caller's existing
(oriented allele index, Phred quality, four likelihoods) observations, deduplicated
by molecule. Mixture fractions have a discrete uniform prior, not an MLE fit.
"""

import csv
from collections import Counter
from itertools import combinations
from math import exp, log

BASES = "ACGT"
FRACTIONS = tuple(i / 20 for i in range(1, 20))
PAIRS = tuple(combinations(range(4), 2))


def log_mean_exp(values):
    maximum = max(values)
    return maximum + log(sum(exp(v - maximum) for v in values) / len(values))


def mixture_evidence(observations):
    """Natural-log BF: six pairs x 19 fractions versus four single bases.

    Each pair and fraction has equal prior mass within the mixture model; each
    base has equal mass within the single model. No between-model prior is set.
    Fractions span .05 through .95, so this does not model arbitrarily rare alleles.
    Empty observations give BF=1 (log BF=0), i.e. no evidence.
    """
    vectors = Counter(tuple(observation[2]) for observation in observations)
    single = [sum(count * log(max(vector[a], 1e-300)) for vector, count in vectors.items()) for a in range(4)]
    pair_logs = []
    for a, b in PAIRS:
        fractions = [sum(count * log(max(f * vector[a] + (1 - f) * vector[b], 1e-300))
                         for vector, count in vectors.items()) for f in FRACTIONS]
        pair_logs.append(log_mean_exp(fractions))
    best_pair = PAIRS[max(range(len(PAIRS)), key=pair_logs.__getitem__)]
    return {
        "log_bf_mixture_vs_single": log_mean_exp(pair_logs) - log_mean_exp(single),
        "best_single": BASES[max(range(4), key=single.__getitem__)],
        "best_pair": "".join(BASES[a] for a in best_pair),
    }


def robust_alleles(column):
    """Neighbour evidence uses the existing Q20 / alternative likelihood <.01 rule."""
    return {molecule: observation[0] for molecule, observation in column.items()
            if observation is not None and observation[1] >= 20
            and max(p for a, p in enumerate(observation[2]) if a != observation[0]) < 0.01}


def diagnostic_writers(mixture_report, linkage_report):
    mixture = csv.writer(mixture_report, delimiter="\t", lineterminator="\n") if mixture_report is not None else None
    linkage = csv.writer(linkage_report, delimiter="\t", lineterminator="\n") if linkage_report is not None else None
    if mixture is not None:
        mixture.writerow(["contig_id", "position_0based", "before", "after", "molecules",
                          "observed_counts_ACGT", "robust_counts_ACGT", "best_single", "best_pair",
                          "log_bf_mixture_vs_single", "support_status", "linked_neighbour_sites"])
    if linkage is not None:
        linkage.writerow(["contig_id", "focal_position_0based", "neighbour_position_0based",
                          "focal_observed_allele", "neighbour_robust_allele", "molecules"])
    return mixture, linkage


def write_diagnostics(columns, contig_id, before, after, mixture_writer, linkage_writer):
    """Report every column with >=2 usable observed alleles, including retained sites.

    Neighbours are selected independently of the focal allele: at least two robust
    molecules per allele for at least two alleles. Co-occurrences require a Q20
    focal observation and a robust neighbour on the same molecule. No phasing,
    significance test, read reassignment or claim of independent neighbouring sites.
    """
    neighbours = {}
    for position, column in enumerate(columns):
        robust = robust_alleles(column)
        if sum(count >= 2 for count in Counter(robust.values()).values()) >= 2:
            neighbours[position] = robust
    for position, column in enumerate(columns):
        observations = [o for o in column.values() if o is not None]
        counts = Counter(o[0] for o in observations)
        if len(counts) < 2:
            continue
        evidence = mixture_evidence(observations)
        focal = {m: o[0] for m, o in column.items() if o is not None and o[1] >= 20}
        linked_sites = 0
        for neighbour, robust in neighbours.items():
            if neighbour == position:
                continue
            joint = Counter((allele, robust[m]) for m, allele in focal.items() if m in robust)
            # Require observed variation at both sites among the shared molecules.
            if len({a for a, _ in joint}) < 2 or len({b for _, b in joint}) < 2:
                continue
            linked_sites += 1
            if linkage_writer is not None:
                for (a, b), count in sorted(joint.items()):
                    linkage_writer.writerow([contig_id, position, neighbour, BASES[a], BASES[b], count])
        if mixture_writer is not None:
            robust_counts = Counter(robust_alleles(column).values())
            mixture_writer.writerow([contig_id, position, before[position], after[position], len(observations),
                ",".join(str(counts[a]) for a in range(4)),
                ",".join(str(robust_counts[a]) for a in range(4)),
                evidence["best_single"], evidence["best_pair"], f'{evidence["log_bf_mixture_vs_single"]:.10g}',
                "insufficient_support" if len(observations) < 3 else "diagnostic_only", linked_sites])
