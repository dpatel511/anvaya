"""Conservative report-only solid-k-mer read-correction candidates."""

import csv
from dataclasses import dataclass
from pathlib import Path

from anvaya.bidirected import _positioned_encoded_kmers
from anvaya.reads import Read


@dataclass(slots=True, frozen=True)
class CorrectionSummary:
    """Counts from a non-destructive correction-candidate audit."""

    reads: int = 0
    audited_reads: int = 0
    reads_without_quality: int = 0
    low_quality_bases: int = 0
    would_correct: int = 0
    protected_damage: int = 0
    ambiguous: int = 0
    incomplete_rescue: int = 0


def _code(sequence: str) -> int | None:
    encoded = _positioned_encoded_kmers(sequence, len(sequence))
    return next(encoded, (0, None, 0))[1]


def _supports(
    sequence: str,
    position: int,
    k: int,
    counts: dict[int, int],
) -> list[int]:
    first = max(0, position - k + 1)
    last = min(position, len(sequence) - k)
    supports = []
    for start in range(first, last + 1):
        code = _code(sequence[start : start + k])
        supports.append(0 if code is None else counts.get(code, 0))
    return supports


def _damage_reversal(
    observed: str,
    proposed: str,
    position: int,
    read_length: int,
    end_window: int,
) -> bool:
    return (
        observed == "T" and proposed == "C" and position < end_window
    ) or (
        observed == "A"
        and proposed == "G"
        and read_length - 1 - position < end_window
    )


def write_correction_report(
    reads: list[Read],
    path: Path,
    *,
    k: int,
    min_count: int,
    maximum_quality: int = 20,
    end_window: int = 5,
    maximum_reads: int = 0,
) -> CorrectionSummary:
    """Write auditable substitution candidates without changing reads."""
    if k < 2:
        raise ValueError("k must be at least 2")
    if min_count < 2:
        raise ValueError("correction audit requires min_count of at least 2")
    if not 0 <= maximum_quality <= 254:
        raise ValueError("maximum_quality must be between 0 and 254")
    if maximum_reads < 0:
        raise ValueError("maximum_reads must not be negative")
    counts: dict[int, int] = {}
    for read in reads:
        if len(read.sequence) < k:
            continue
        for _, code, _ in _positioned_encoded_kmers(read.sequence, k):
            counts[code] = counts.get(code, 0) + 1

    rows: list[dict[str, object]] = []
    reads_without_quality = 0
    low_quality_bases = 0
    decisions = {"would_correct": 0, "protect_damage": 0, "ambiguous": 0, "incomplete_rescue": 0}
    if maximum_reads and maximum_reads < len(reads):
        selected = [
            (
                index * len(reads) // maximum_reads,
                reads[index * len(reads) // maximum_reads],
            )
            for index in range(maximum_reads)
        ]
    else:
        selected = list(enumerate(reads))
    for read_index, read in selected:
        if read.qualities is None:
            reads_without_quality += 1
            continue
        if len(read.sequence) < k:
            continue
        for position, quality in enumerate(read.qualities):
            if quality > maximum_quality:
                continue
            low_quality_bases += 1
            original_supports = _supports(read.sequence, position, k, counts)
            original_solid = sum(value >= min_count for value in original_supports)
            if original_solid == len(original_supports):
                continue
            candidates: list[tuple[tuple[int, int, int], str, list[int]]] = []
            for base in "ACGT":
                if base == read.sequence[position]:
                    continue
                sequence = (
                    read.sequence[:position] + base + read.sequence[position + 1 :]
                )
                supports = _supports(sequence, position, k, counts)
                score = (
                    sum(value >= min_count for value in supports),
                    min(supports, default=0),
                    sum(supports),
                )
                candidates.append((score, base, supports))
            best_score = max(candidate[0] for candidate in candidates)
            best = [candidate for candidate in candidates if candidate[0] == best_score]
            if best_score[0] <= original_solid:
                continue
            proposed = best[0][1] if len(best) == 1 else ""
            if len(best) > 1:
                decision = "ambiguous"
                reason = "multiple_equal_solid_kmer_rescues"
            elif _damage_reversal(
                read.sequence[position],
                proposed,
                position,
                len(read.sequence),
                end_window,
            ):
                decision = "protect_damage"
                reason = "terminal_damage_compatible_reversal"
            elif best_score[0] == len(best[0][2]):
                decision = "would_correct"
                reason = "unique_complete_solid_kmer_rescue"
            else:
                decision = "incomplete_rescue"
                reason = "unique_but_not_all_kmers_rescued"
            decisions[decision] += 1
            rows.append(
                {
                    "read_index": read_index,
                    "read_name": read.name,
                    "position": position,
                    "observed_base": read.sequence[position],
                    "proposed_base": proposed,
                    "base_quality": quality,
                    "decision": decision,
                    "reason": reason,
                    "overlapping_kmers": len(original_supports),
                    "original_solid_kmers": original_solid,
                    "proposed_solid_kmers": best_score[0],
                    "original_min_support": min(original_supports, default=0),
                    "proposed_min_support": best_score[1],
                }
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "read_index",
        "read_name",
        "position",
        "observed_base",
        "proposed_base",
        "base_quality",
        "decision",
        "reason",
        "overlapping_kmers",
        "original_solid_kmers",
        "proposed_solid_kmers",
        "original_min_support",
        "proposed_min_support",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return CorrectionSummary(
        reads=len(reads),
        audited_reads=len(selected),
        reads_without_quality=reads_without_quality,
        low_quality_bases=low_quality_bases,
        would_correct=decisions["would_correct"],
        protected_damage=decisions["protect_damage"],
        ambiguous=decisions["ambiguous"],
        incomplete_rescue=decisions["incomplete_rescue"],
    )
