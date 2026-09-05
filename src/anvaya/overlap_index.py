"""Canonical anchor sketches and sequence identity for overlap discovery."""

import heapq
from collections import defaultdict

from anvaya.reads import Read
from anvaya.sequences import reverse_complement


def _ry(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "RYRYN"))


def _identity(first: str, second: str) -> float:
    if not first:
        return 0.0
    return sum(a == b for a, b in zip(first, second, strict=True)) / len(first)


def _encoded_kmer(sequence: str) -> int | None:
    value = 0
    for base in sequence:
        code = "ACGT".find(base)
        if code < 0:
            return None
        value = (value << 2) | code
    return value


def _sketch_anchors(
    sequence: str, anchor_k: int, end_window: int, anchors_per_read: int
) -> list[tuple[int, int, bool]]:
    """Return a bounded bottom sketch of canonical internal anchors."""
    start = end_window
    stop = len(sequence) - end_window - anchor_k + 1
    if start >= stop:
        return []
    initial = sequence[start : start + anchor_k]
    forward = _encoded_kmer(initial)
    reverse = _encoded_kmer(reverse_complement(initial))
    if "N" in sequence[start : len(sequence) - end_window]:
        # Ambiguous bases are uncommon; keep the simpler safe path for them.
        candidates = []
        seen: set[int] = set()
        for position in range(start, stop):
            anchor = sequence[position : position + anchor_k]
            encoded = _encoded_kmer(anchor)
            encoded_reverse = _encoded_kmer(reverse_complement(anchor))
            if encoded is None or encoded_reverse is None:
                continue
            canonical = min(encoded, encoded_reverse)
            if canonical in seen:
                continue
            seen.add(canonical)
            candidates.append(
                (canonical, position, encoded_reverse < encoded)
            )
        return heapq.nsmallest(anchors_per_read, candidates)

    assert forward is not None and reverse is not None

    candidates: list[tuple[int, int, bool]] = []
    seen: set[int] = set()
    mask = (1 << (2 * anchor_k)) - 1
    high_shift = 2 * (anchor_k - 1)
    for position in range(start, stop):
        canonical = min(forward, reverse)
        if canonical in seen:
            pass
        else:
            seen.add(canonical)
            candidates.append((canonical, position, reverse < forward))
        next_position = position + anchor_k
        if next_position < len(sequence) - end_window:
            code = "ACGT".find(sequence[next_position])
            forward = ((forward << 2) & mask) | code
            reverse = ((3 - code) << high_shift) | (reverse >> 2)
    return heapq.nsmallest(anchors_per_read, candidates)


def _anchor_index(
    reads: list[Read],
    anchor_k: int,
    end_window: int,
    anchors_per_read: int,
    maximum_anchor_occurrences: int,
) -> dict[int, list[int]]:
    """Index a fixed number of packed canonical anchors per read."""
    index: dict[int, list[int]] = defaultdict(list)
    position_bits = max(1, max(len(read.sequence) for read in reads).bit_length())
    payload_bits = position_bits + 1
    for read_index, read in enumerate(reads):
        for anchor, position, reverse in _sketch_anchors(
            read.sequence, anchor_k, end_window, anchors_per_read
        ):
            occurrences = index[anchor]
            if len(occurrences) <= maximum_anchor_occurrences:
                occurrences.append(
                    (read_index << payload_bits) | (position << 1) | reverse
                )
    return dict(index)
