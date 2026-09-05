"""Conservative direct merging of overlapping paired-end reads."""

from dataclasses import dataclass

from anvaya.damage_consensus import _identity
from anvaya.overlap_assembly import _ry
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


@dataclass(frozen=True, slots=True)
class PairedMergeDiagnostics:
    """Counts from splitting paired reads into merged and unmerged evidence."""

    input_pairs: int = 0
    merged_pairs: int = 0
    unmerged_pairs: int = 0
    ambiguous_pairs: int = 0
    merged_bases: int = 0


def _reverse_qualities(read: Read) -> tuple[int, ...] | None:
    return None if read.qualities is None else tuple(reversed(read.qualities))


def _merge_at_offset(left: Read, right: Read, offset: int) -> Read:
    right_sequence = reverse_complement(right.sequence)
    right_qualities = _reverse_qualities(right)
    start = min(0, offset)
    end = max(len(left.sequence), offset + len(right_sequence))
    sequence: list[str] = []
    qualities: list[int] | None = (
        [] if left.qualities is not None and right_qualities is not None else None
    )
    for position in range(start, end):
        left_index = position
        right_index = position - offset
        left_base = (
            left.sequence[left_index]
            if 0 <= left_index < len(left.sequence)
            else None
        )
        right_base = (
            right_sequence[right_index]
            if 0 <= right_index < len(right_sequence)
            else None
        )
        left_quality = (
            left.qualities[left_index]
            if left_base is not None and left.qualities is not None
            else -1
        )
        right_quality = (
            right_qualities[right_index]
            if right_base is not None and right_qualities is not None
            else -1
        )
        if left_base is None:
            base, quality = right_base, right_quality
        elif right_base is None:
            base, quality = left_base, left_quality
        elif left_base == right_base:
            base, quality = left_base, max(left_quality, right_quality)
        elif right_quality > left_quality:
            base, quality = right_base, right_quality
        else:
            base, quality = left_base, left_quality
        assert base is not None
        sequence.append(base)
        if qualities is not None:
            qualities.append(quality)
    return Read(
        f"{left.name}|merged",
        "".join(sequence),
        None if qualities is None else tuple(qualities),
    )


def merge_overlapping_pairs(
    left_reads: list[Read],
    right_reads: list[Read],
    *,
    minimum_overlap: int = 20,
    minimum_identity: float = 0.95,
    minimum_ry_identity: float = 0.99,
) -> tuple[list[Read], list[int], PairedMergeDiagnostics]:
    """Merge uniquely best R1/R2 overlaps and retain both mates otherwise."""
    if len(left_reads) != len(right_reads):
        raise ValueError("paired-end inputs must contain the same number of reads")
    if minimum_overlap < 1:
        raise ValueError("minimum_overlap must be at least 1")
    if not 0.0 <= minimum_identity <= 1.0:
        raise ValueError("minimum_identity must be between 0 and 1")
    if not 0.0 <= minimum_ry_identity <= 1.0:
        raise ValueError("minimum_ry_identity must be between 0 and 1")

    output: list[Read] = []
    molecule_ids: list[int] = []
    merged_pairs = ambiguous_pairs = merged_bases = 0
    for molecule, (left, right) in enumerate(zip(left_reads, right_reads)):
        reverse_right = reverse_complement(right.sequence)
        candidates: list[tuple[tuple[int, float, float], int]] = []
        minimum_offset = minimum_overlap - len(reverse_right)
        maximum_offset = len(left.sequence) - minimum_overlap
        for offset in range(minimum_offset, maximum_offset + 1):
            left_start = max(0, offset)
            right_start = max(0, -offset)
            overlap = min(
                len(left.sequence) - left_start,
                len(reverse_right) - right_start,
            )
            if overlap < minimum_overlap:
                continue
            left_overlap = left.sequence[left_start:left_start + overlap]
            right_overlap = reverse_right[right_start:right_start + overlap]
            dna_identity = _identity(left_overlap, right_overlap)
            ry_identity = _identity(_ry(left_overlap), _ry(right_overlap))
            if (
                dna_identity >= minimum_identity
                and ry_identity >= minimum_ry_identity
            ):
                candidates.append(
                    ((overlap, ry_identity, dna_identity), offset)
                )
        if candidates:
            best_score = max(score for score, _ in candidates)
            best_offsets = [
                offset for score, offset in candidates if score == best_score
            ]
        else:
            best_offsets = []
        if len(best_offsets) == 1:
            merged = _merge_at_offset(left, right, best_offsets[0])
            output.append(merged)
            molecule_ids.append(molecule)
            merged_pairs += 1
            merged_bases += len(merged.sequence)
        else:
            output.extend((left, right))
            molecule_ids.extend((molecule, molecule))
            ambiguous_pairs += int(len(best_offsets) > 1)

    return output, molecule_ids, PairedMergeDiagnostics(
        input_pairs=len(left_reads),
        merged_pairs=merged_pairs,
        unmerged_pairs=len(left_reads) - merged_pairs,
        ambiguous_pairs=ambiguous_pairs,
        merged_bases=merged_bases,
    )
