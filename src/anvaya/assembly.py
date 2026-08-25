"""End-to-end assembly operations."""

from collections.abc import Sequence

from anvaya.graph import build_dbg
from anvaya.unitigs import extract_unitigs


def assemble(reads: Sequence[str], k: int) -> list[str]:
    """Assemble reads into maximal non-branching unitig sequences."""
    return extract_unitigs(build_dbg(reads, k))
