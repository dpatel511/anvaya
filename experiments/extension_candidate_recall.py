"""Small synthetic truth-placement recall benchmark; 0-based ungapped offsets."""

import argparse
import csv
import hashlib
import itertools
import json
import random
from collections import Counter
from pathlib import Path

from anvaya.overlap_assembly import _candidate_alignments
from anvaya.overlap_index import _anchor_index
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


def case(seed, length, overlap, damage, error, reverse):
    rng = random.Random(seed)
    truth = "".join(rng.choice("ACGT") for _ in range(2 * length - overlap))
    def corrupt(sequence):
        result = []
        for p, base in enumerate(sequence):
            distance = p if base == "C" else len(sequence) - 1 - p
            if base in "CG" and distance < 5 and rng.random() < damage * .65 ** distance:
                base = "T" if base == "C" else "A"
            if rng.random() < error:
                base = rng.choice([b for b in "ACGT" if b != base])
            result.append(base)
        return "".join(result)
    target = corrupt(truth[:length])
    partner = truth[length-overlap:]
    if reverse:
        partner = reverse_complement(partner)
    reads = [Read("partner", corrupt(partner))]
    reads.extend(Read(f"decoy_{i}", "".join(rng.choice("ACGT") for _ in range(length))) for i in range(3))
    anchors = _anchor_index(reads, 15, 0, 8, 100)
    trace = {}
    candidates = _candidate_alignments(target, reads, list(range(4)), anchors, set(),
        anchor_k=15, anchors_per_read=8, maximum_anchor_occurrences=100,
        minimum_anchor_matches=1, minimum_overlap=30, minimum_identity=.90,
        minimum_ry_identity=.99, position_bits=length.bit_length(), target_window=length,
        stage_trace=trace)
    truth_key = (0, reverse, length-overlap)
    selected = {(c.read_index, c.reverse, c.offset) for c in candidates}
    return {"stage": trace.get(truth_key, "no_true_offset_anchor_vote"),
            "truth_selected": int(truth_key in selected),
            "other_placements_selected": len(selected - {truth_key})}


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    totals, strata = Counter(), {}
    case_count = 0
    with (output / "cases.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = None
        for params in itertools.product(range(91001, 91021), (50, 100), (30, 40), (0, .4), (0, .01), (False, True)):
            seed, length, overlap, damage, error, reverse = params
            result = case(*params)
            case_count += 1
            row = dict(seed=seed, length=length, overlap=overlap, damage=damage, error=error, reverse=reverse, **result)
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
                writer.writeheader()
            writer.writerow(row)
            totals[result["stage"]] += 1
            key = f"damage={damage},error={error}"
            strata.setdefault(key, Counter())[result["stage"]] += 1
            totals["other_placements_selected"] += result["other_placements_selected"]
    root = Path(__file__).resolve().parents[1]
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in
              (Path(__file__), root / "src/anvaya/overlap_assembly.py", root / "src/anvaya/overlap_index.py")}
    if sum(v for k, v in totals.items() if k != "other_placements_selected") != case_count:
        raise AssertionError("case/stage accounting mismatch")
    summary = dict(cases=case_count, counts=dict(totals), strata={k:dict(v) for k,v in strata.items()}, source_sha256=hashes)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
