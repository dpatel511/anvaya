#!/usr/bin/env python3
"""Split assembled records from unresolved fragments using raw molecule counts."""

import argparse
import csv
import hashlib
import json
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _n50(lengths: list[int]) -> int:
    threshold = (sum(lengths) + 1) // 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            return length
    return 0


def _molecule_counts(path: Path) -> dict[str, int]:
    counts = {}
    closed = set()
    current = None
    molecules = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"contig_id", "molecule_id"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("placements require contig_id and molecule_id columns")
        for row in reader:
            contig = row["contig_id"]
            molecule = row["molecule_id"]
            if not contig or not molecule:
                raise ValueError("placement contig and molecule IDs must not be empty")
            if contig != current:
                if current is not None:
                    counts[current] = len(molecules)
                    closed.add(current)
                if contig in closed:
                    raise ValueError(f"placements are not grouped by contig: {contig}")
                current = contig
                molecules = set()
            molecules.add(molecule)
    if current is not None:
        counts[current] = len(molecules)
    if not counts:
        raise ValueError("placements must contain at least one record")
    return counts


def split_assembled_output(
    fasta: Path,
    placements: Path,
    assembled: Path,
    unresolved: Path,
    report: Path,
    minimum_molecules: int = 2,
) -> dict:
    if minimum_molecules < 2:
        raise ValueError("minimum molecules must be at least 2")
    paths = [fasta.resolve(), placements.resolve(), assembled.resolve(),
             unresolved.resolve(), report.resolve()]
    if len(paths) != len(set(paths)):
        raise ValueError("inputs and outputs must be distinct")
    if not fasta.is_file() or not placements.is_file():
        raise ValueError("FASTA and placements inputs must exist")
    for output in (assembled, unresolved, report):
        if output.exists():
            raise ValueError(f"refusing to overwrite {output}")

    counts = _molecule_counts(placements)
    seen = set()
    assembled_lengths = []
    unresolved_lengths = []
    current_id = None
    current_lines = []
    current_length = 0

    def emit(assembled_handle, unresolved_handle) -> None:
        nonlocal current_id, current_lines, current_length
        if current_id is None:
            return
        if current_id in seen:
            raise ValueError(f"duplicate FASTA identifier: {current_id}")
        if current_id not in counts:
            raise ValueError(f"FASTA record lacks placement evidence: {current_id}")
        seen.add(current_id)
        if counts[current_id] >= minimum_molecules:
            assembled_handle.writelines(current_lines)
            assembled_lengths.append(current_length)
        else:
            unresolved_handle.writelines(current_lines)
            unresolved_lengths.append(current_length)

    for output in (assembled, unresolved, report):
        output.parent.mkdir(parents=True, exist_ok=True)
    with fasta.open(encoding="utf-8") as source, assembled.open(
        "w", encoding="utf-8"
    ) as assembled_handle, unresolved.open("w", encoding="utf-8") as unresolved_handle:
        for raw_line in source:
            if raw_line.startswith(">"):
                emit(assembled_handle, unresolved_handle)
                current_id = raw_line[1:].split()[0]
                if not current_id:
                    raise ValueError("FASTA identifiers must not be empty")
                current_lines = [raw_line]
                current_length = 0
            else:
                if current_id is None:
                    raise ValueError("sequence occurs before the first FASTA header")
                current_lines.append(raw_line)
                current_length += len(raw_line.strip())
        emit(assembled_handle, unresolved_handle)

    missing = set(counts) - seen
    if missing:
        raise ValueError(f"placements reference {len(missing)} absent FASTA records")
    result = {
        "minimum_distinct_molecules": minimum_molecules,
        "coordinate_convention": "placement coordinates are 0-based, half-open; unused",
        "input": str(fasta.resolve()),
        "input_sha256": _sha256(fasta),
        "placements": str(placements.resolve()),
        "placements_sha256": _sha256(placements),
        "input_records": len(seen),
        "input_bases": sum(assembled_lengths) + sum(unresolved_lengths),
        "assembled": {
            "path": str(assembled.resolve()),
            "sha256": _sha256(assembled),
            "records": len(assembled_lengths),
            "bases": sum(assembled_lengths),
            "n50": _n50(assembled_lengths),
            "longest": max(assembled_lengths, default=0),
        },
        "unresolved": {
            "path": str(unresolved.resolve()),
            "sha256": _sha256(unresolved),
            "records": len(unresolved_lengths),
            "bases": sum(unresolved_lengths),
            "n50": _n50(unresolved_lengths),
            "longest": max(unresolved_lengths, default=0),
        },
    }
    report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--placements", required=True, type=Path)
    parser.add_argument("--assembled", required=True, type=Path)
    parser.add_argument("--unresolved", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--minimum-molecules", type=int, default=2)
    args = parser.parse_args()
    result = split_assembled_output(
        args.fasta, args.placements, args.assembled, args.unresolved,
        args.report, args.minimum_molecules,
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
