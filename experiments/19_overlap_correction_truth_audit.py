"""Classify experimental overlap corrections against a known reference."""

import argparse
import csv
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


def _read_events(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        events = list(csv.DictReader(handle, delimiter="\t"))
    for event in events:
        sequence_after = list(event["sequence_before"])
        position = int(event["position"])
        sequence_after[position] = event["corrected_base"]
        event["sequence_after"] = "".join(sequence_after)
    return events


def _write_queries(events: list[dict[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for index, event in enumerate(events):
            handle.write(f">event_{index}_before\n{event['sequence_before']}\n")
            handle.write(f">event_{index}_after\n{event['sequence_after']}\n")


def _parse_paf(path: Path) -> dict[str, list[dict[str, object]]]:
    alignments: dict[str, list[dict[str, object]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            tags = {
                field[:2]: field[5:]
                for field in fields[12:]
                if len(field) >= 6 and field[2:5] in {":i:", ":Z:"}
            }
            query_length = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            if query_start != 0 or query_end != query_length or "NM" not in tags:
                continue
            alignments[fields[0]].append({
                "target": fields[5],
                "strand": fields[4],
                "target_start": int(fields[7]),
                "target_end": int(fields[8]),
                "matches": int(fields[9]),
                "block_length": int(fields[10]),
                "edit_distance": int(tags["NM"]),
            })
    return alignments


def _unique_best(
    alignments: list[dict[str, object]],
) -> dict[str, object] | None:
    if not alignments:
        return None
    score = lambda alignment: (
        int(alignment["matches"]),
        -int(alignment["edit_distance"]),
        -int(alignment["block_length"]),
    )
    best_score = max(score(alignment) for alignment in alignments)
    best = [alignment for alignment in alignments if score(alignment) == best_score]
    loci = {
        (
            alignment["target"],
            alignment["strand"],
            alignment["target_start"],
            alignment["target_end"],
        )
        for alignment in best
    }
    return best[0] if len(loci) == 1 else None


def _classify(
    before_alignments: list[dict[str, object]],
    after_alignments: list[dict[str, object]],
) -> tuple[str, int | None, int | None]:
    before = _unique_best(before_alignments)
    after = _unique_best(after_alignments)
    if before is None or after is None:
        label = "unstable_or_unmapped" if before_alignments or after_alignments else "unmapped"
        return label, None, None
    before_locus = (
        before["target"], before["strand"],
        before["target_start"], before["target_end"],
    )
    after_locus = (
        after["target"], after["strand"],
        after["target_start"], after["target_end"],
    )
    before_edits = int(before["edit_distance"])
    after_edits = int(after["edit_distance"])
    if before_locus != after_locus:
        return "unstable_or_unmapped", before_edits, after_edits
    if after_edits < before_edits:
        return "fixes_mismatch", before_edits, after_edits
    if after_edits > before_edits:
        return "introduces_mismatch", before_edits, after_edits
    return "neutral", before_edits, after_edits


def run_audit(
    correction_report: Path,
    reference: Path,
    output_dir: Path,
    minimap2: str = "minimap2",
) -> Counter[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    events = _read_events(correction_report)
    queries = output_dir / "correction-queries.fasta"
    paf = output_dir / "correction-alignments.paf"
    _write_queries(events, queries)
    with paf.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                minimap2, "-x", "sr", "-k", "9", "-w", "5", "-N", "5",
                "-c", str(reference), str(queries),
            ],
            check=True,
            stdout=handle,
        )
    alignments = _parse_paf(paf)
    counts: Counter[str] = Counter()
    fieldnames = list(events[0]) if events else []
    fieldnames += ["classification", "before_edit_distance", "after_edit_distance"]
    with (output_dir / "correction-truth-audit.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for index, event in enumerate(events):
            label, before_edits, after_edits = _classify(
                alignments.get(f"event_{index}_before", []),
                alignments.get(f"event_{index}_after", []),
            )
            counts[label] += 1
            writer.writerow({
                **event,
                "classification": label,
                "before_edit_distance": "" if before_edits is None else before_edits,
                "after_edit_distance": "" if after_edits is None else after_edits,
            })
    resolved = counts["fixes_mismatch"] + counts["introduces_mismatch"]
    summary = {
        "events": len(events),
        **dict(sorted(counts.items())),
        "resolved_precision": (
            counts["fixes_mismatch"] / resolved if resolved else None
        ),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return counts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--correction-report", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--minimap2", default="minimap2")
    arguments = parser.parse_args()
    counts = run_audit(
        arguments.correction_report,
        arguments.reference,
        arguments.output_dir,
        arguments.minimap2,
    )
    print(f"events={sum(counts.values())}")
    for label, count in sorted(counts.items()):
        print(f"{label}={count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
