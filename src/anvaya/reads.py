"""FASTA and FASTQ input."""

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from anvaya.sequences import normalize_dna


@dataclass(frozen=True)
class Read:
    """A sequencing read and its optional Phred quality scores."""

    name: str
    sequence: str
    qualities: tuple[int, ...] | None = None


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _load_fasta(first_header: str, handle: TextIO) -> list[Read]:
    reads: list[Read] = []
    name = first_header[1:].strip()
    sequence_parts: list[str] = []

    if not name:
        raise ValueError("FASTA header must contain a read name")

    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            reads.append(Read(name, normalize_dna("".join(sequence_parts))))
            name = line[1:].strip()
            sequence_parts = []
            if not name:
                raise ValueError("FASTA header must contain a read name")
        else:
            sequence_parts.append(line)

    reads.append(Read(name, normalize_dna("".join(sequence_parts))))
    return reads


def _load_fastq(first_header: str, handle: TextIO) -> list[Read]:
    reads: list[Read] = []
    header = first_header

    while header:
        if not header.startswith("@"):
            raise ValueError("FASTQ record must start with '@'")

        name = header[1:].strip()
        sequence_line = handle.readline()
        separator = handle.readline()
        quality_line = handle.readline()
        if not name or not sequence_line or not separator or not quality_line:
            raise ValueError("incomplete FASTQ record")
        if not separator.startswith("+"):
            raise ValueError("FASTQ separator must start with '+'")

        sequence = normalize_dna(sequence_line.strip())
        quality = quality_line.rstrip("\r\n")
        if len(quality) != len(sequence):
            raise ValueError("FASTQ sequence and quality lengths must match")

        scores = tuple(ord(character) - 33 for character in quality)
        if any(score < 0 for score in scores):
            raise ValueError("FASTQ contains an invalid quality character")
        reads.append(Read(name, sequence, scores))

        header = handle.readline()
        while header and not header.strip():
            header = handle.readline()
        header = header.rstrip("\r\n")

    return reads


def load_reads(path: str | Path) -> list[Read]:
    """Load reads from a plain or gzipped FASTA/FASTQ file."""
    input_path = Path(path)
    with _open_text(input_path) as handle:
        first_line = handle.readline().rstrip("\r\n")
        if not first_line:
            raise ValueError("sequence file must not be empty")
        if first_line.startswith(">"):
            return _load_fasta(first_line, handle)
        if first_line.startswith("@"):
            return _load_fastq(first_line, handle)
        raise ValueError("input must be FASTA or FASTQ")
