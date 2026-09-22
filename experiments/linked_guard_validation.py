"""Fixed-layout synthetic validation, not an end-to-end assembly benchmark.

Known source haplotype 0 is the target. Haplotype 1 is competing recruitment.
Coordinates are 0-based, half-open. All fixtures are generated in memory.
"""

import argparse
import csv
import hashlib
import itertools
import json
import random
from collections import Counter
from dataclasses import replace
from pathlib import Path

from anvaya.overlap_progressive import ProgressiveSequencePool
from anvaya.raw_consensus import DamageProfile, RawPlacement, project_raw_consensus
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


def fixture(seed, depth, competitor_fraction, damage, profile_mode, separation, quality):
    rng = random.Random(seed)
    truth = list("".join(rng.choice("ACGT") for _ in range(120)))
    truth[40], truth[40 + separation] = "T", "C"
    truth = "".join(truth)
    competitor = list(truth)
    competitor[40], competitor[40 + separation] = "C", "T"
    competitor = "".join(competitor)
    baseline = list(truth)
    # Independent baseline errors prevent a guard that simply preserves everything
    # from looking optimal. These are not simulated earlier assembly decisions.
    for p, base in enumerate(baseline):
        if rng.random() < .02:
            baseline[p] = rng.choice([b for b in "ACGT" if b != base])
    reads, placements, covered = [], [], set()
    competitors = 0
    for i in range(depth):
        is_competitor = rng.random() < competitor_fraction
        competitors += is_competitor
        haplotype = competitor if is_competitor else truth
        start = rng.randrange(15, 41)
        stop = start + rng.randrange(31, 66)
        reverse = bool(rng.randrange(2))
        raw = haplotype[start:stop]
        if reverse:
            raw = reverse_complement(raw)
        sequence = list(raw)
        for p, base in enumerate(sequence):
            distance = p if base == "C" else len(raw) - 1 - p
            if base in "CG" and distance < 5 and rng.random() < damage * (.65 ** distance):
                base = "T" if base == "C" else "A"
            if rng.random() < 10 ** (-quality / 10):
                base = rng.choice([b for b in "ACGT" if b != base])
            sequence[p] = base
        reads.append(Read(str(i), "".join(sequence), (quality,) * len(raw)))
        placements.append(RawPlacement(i, start, reverse, 0, len(raw)))
        covered.update(range(start, stop))
    supplied = damage if profile_mode == "matched" else min(.9, damage + .3)
    profile = DamageProfile(tuple(supplied * .65 ** p for p in range(5)),
                            tuple(supplied * .65 ** p for p in range(5)))
    pool = ProgressiveSequencePool.from_reads(reads)
    center = replace(pool.records[0].corrected(Read("synthetic", "".join(baseline))),
                     raw_placements=tuple(placements))
    return pool.replace_record(center), profile, truth, covered, competitors


def evaluate(pool, profile, truth, covered, variants=()):
    before = pool.records[0].current.sequence
    unguarded, _ = project_raw_consensus(pool, profile)
    guarded, diagnostics = project_raw_consensus(pool, profile, linked_allele_guard=True)
    old, new = unguarded[0].sequence, guarded[0].sequence
    if len(old) != len(truth) or len(new) != len(truth):
        raise AssertionError("layout length changed")
    changed = [p for p in covered if old[p] != new[p]]
    if len(changed) != diagnostics.linked_allele_rejections or any(new[p] != before[p] for p in changed):
        raise AssertionError("guard output/decision mismatch")
    return {
        "covered_positions": len(covered),
        "baseline_errors": sum(before[p] != truth[p] for p in covered),
        "unguarded_errors": sum(old[p] != truth[p] for p in covered),
        "guarded_errors": sum(new[p] != truth[p] for p in covered),
        "blocked_changes": len(changed),
        "prevented_errors": sum(old[p] != truth[p] and new[p] == truth[p] for p in changed),
        "blocked_beneficial_corrections": sum(old[p] == truth[p] and new[p] != truth[p] for p in changed),
        "both_wrong": sum(old[p] != truth[p] and new[p] != truth[p] for p in changed),
        "source_variant_sites": len(variants),
        "unguarded_source_alleles": sum(old[p] == truth[p] for p in variants),
        "guarded_source_alleles": sum(new[p] == truth[p] for p in variants),
    }


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    sources = [root / "src/anvaya/linked_allele_guard.py", root / "src/anvaya/raw_consensus.py", Path(__file__)]
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    totals, strata = Counter(), {}
    with (output / "cases.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = None
        for params in itertools.product(range(73001, 73009), (6, 12, 24), (0, .2, .5, .8),
                                        (0, .2, .4), ("matched", "overestimated"), (6, 35), (20, 30)):
            seed, depth, fraction, damage, mode, separation, quality = params
            pool, profile, truth, covered, realized = fixture(*params)
            # Retention refers to the designated source allele, not recovery of both haplotypes.
            variants = [p for p in (40, 40 + separation) if p in covered and fraction > 0]
            metrics = evaluate(pool, profile, truth, covered, variants)
            row = dict(seed=seed, depth=depth, competitor_fraction=fraction, realized_competitors=realized,
                       damage=damage, profile=mode, variant_separation=separation, quality=quality, **metrics)
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
                writer.writeheader()
            writer.writerow(row)
            totals.update(metrics)
            totals["cases"] += 1
            key = f"competitor={fraction},profile={mode},separation={separation}"
            strata.setdefault(key, Counter()).update(metrics)
    if any(hashlib.sha256(p.read_bytes()).hexdigest() != hashes[str(p.relative_to(root))] for p in sources):
        raise AssertionError("benchmark sources changed during evaluation")
    summary = {"source_sha256": hashes, "totals": dict(totals), "strata": {k: dict(v) for k, v in strata.items()}}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(totals), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    run(parser.parse_args().output)
