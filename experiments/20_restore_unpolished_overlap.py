"""Restore the exact pre-polish assembly from an overlap correction report."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from anvaya.output import write_fasta
from anvaya.reads import load_reads


def restore_unpolished(
    polished_contigs: Path,
    correction_report: Path,
    output: Path,
) -> int:
    contigs = load_reads(polished_contigs)
    corrections: dict[int, list[tuple[int, str, str]]] = defaultdict(list)
    with correction_report.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            corrections[int(row["cluster"]) - 1].append((
                int(row["position"]),
                row["original_base"],
                row["corrected_base"],
            ))

    restored = [list(contig.sequence) for contig in contigs]
    for contig_index, events in corrections.items():
        if not 0 <= contig_index < len(restored):
            raise ValueError("correction report does not match contig order")
        for position, original_base, corrected_base in events:
            if not 0 <= position < len(restored[contig_index]):
                raise ValueError("correction position is outside its contig")
            if restored[contig_index][position] != corrected_base:
                raise ValueError("polished contig does not contain the reported correction")
            restored[contig_index][position] = original_base

    write_fasta(("".join(sequence) for sequence in restored), output)
    return sum(len(events) for events in corrections.values())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--polished-contigs", required=True, type=Path)
    parser.add_argument("--correction-report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    restored = restore_unpolished(
        arguments.polished_contigs,
        arguments.correction_report,
        arguments.output,
    )
    print(f"restored_bases={restored}")
    print(f"output={arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
