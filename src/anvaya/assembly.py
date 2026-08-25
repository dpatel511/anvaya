"""End-to-end assembly operations."""

from collections.abc import Sequence
from pathlib import Path

from anvaya.graph import build_dbg
from anvaya.reads import load_reads
from anvaya.unitigs import extract_unitigs


def assemble(reads: Sequence[str], k: int, min_count: int = 1) -> list[str]:
    """Assemble reads into maximal non-branching unitig sequences."""
    return extract_unitigs(build_dbg(reads, k, min_count))


def assemble_file(
    path: str | Path, k: int, min_count: int = 1
) -> list[str]:
    """Load a FASTA/FASTQ file and assemble its observed read orientations."""
    return assemble(
        [read.sequence for read in load_reads(path)],
        k,
        min_count,
    )
