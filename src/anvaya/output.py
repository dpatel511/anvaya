"""Assembler output utilities."""

from collections.abc import Sequence
from pathlib import Path

from anvaya.sequences import normalize_dna


def write_fasta(
    sequences: Sequence[str],
    path: str | Path,
    line_width: int = 80,
) -> None:
    """Write sequences to a FASTA file with deterministic unitig names."""
    if not isinstance(line_width, int) or isinstance(line_width, bool):
        raise TypeError("line_width must be an integer")
    if line_width < 1:
        raise ValueError("line_width must be at least 1")

    output_path = Path(path)
    with output_path.open(mode="w", encoding="utf-8") as handle:
        for index, sequence in enumerate(sequences, start=1):
            normalized = normalize_dna(sequence)
            handle.write(f">unitig_{index} length={len(normalized)}\n")
            for start in range(0, len(normalized), line_width):
                handle.write(normalized[start : start + line_width] + "\n")
