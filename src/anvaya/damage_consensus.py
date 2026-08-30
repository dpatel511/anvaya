"""Conservative damage-aware consensus correction from overlapping reads."""

import csv
import heapq
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from anvaya.reads import Read
from anvaya.sequences import reverse_complement


@dataclass(slots=True, frozen=True)
class DamageConsensusSummary:
    """Counts from damage-aware read consensus."""

    reads: int = 0
    candidate_bases: int = 0
    corrected_bases: int = 0
    corrected_reads: int = 0
    insufficient_support: int = 0
    ambiguous: int = 0


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


def _alignments(
    target_index: int,
    reads: list[Read],
    molecule_ids: list[int],
    anchors: dict[int, list[int]],
    *,
    anchor_k: int,
    end_window: int,
    minimum_anchor_matches: int,
    maximum_anchor_occurrences: int,
    minimum_overlap: int,
    minimum_identity: float,
    minimum_ry_identity: float,
    anchors_per_read: int,
    position_bits: int,
) -> list[tuple[str, int, tuple[int, ...] | None]]:
    target = reads[target_index].sequence
    payload_bits = position_bits + 1
    position_mask = (1 << position_bits) - 1
    votes: dict[tuple[int, bool, int], int] = defaultdict(int)
    for anchor, target_position, target_reverse in _sketch_anchors(
        target, anchor_k, end_window, anchors_per_read
    ):
        occurrences = anchors.get(anchor, ())
        if not occurrences or len(occurrences) > maximum_anchor_occurrences:
            continue
        for packed in occurrences:
            member_index = packed >> payload_bits
            if molecule_ids[member_index] == molecule_ids[target_index]:
                continue
            member_position = (packed >> 1) & position_mask
            reverse = target_reverse != bool(packed & 1)
            if reverse:
                member_position = len(reads[member_index].sequence) - member_position - anchor_k
            votes[(member_index, reverse, target_position - member_position)] += 1

    by_molecule: dict[int, list[tuple[int, int, bool, int]]] = defaultdict(list)
    for (member_index, reverse, offset), support in votes.items():
        if support < minimum_anchor_matches:
            continue
        read = reads[member_index]
        member = reverse_complement(read.sequence) if reverse else read.sequence
        start = max(0, offset)
        stop = min(len(target), offset + len(member))
        if stop - start < minimum_overlap:
            continue
        member_start = start - offset
        target_overlap = target[start:stop]
        member_overlap = member[member_start : member_start + stop - start]
        if _identity(target_overlap, member_overlap) < minimum_identity:
            continue
        if _identity(_ry(target_overlap), _ry(member_overlap)) < minimum_ry_identity:
            continue
        molecule = molecule_ids[member_index]
        by_molecule[molecule].append((support, member_index, reverse, offset))
    selected: list[tuple[str, int, tuple[int, ...] | None]] = []
    for candidates in by_molecule.values():
        best_support = max(item[0] for item in candidates)
        best = [item for item in candidates if item[0] == best_support]
        if len(best) == 1:
            _, member_index, reverse, offset = best[0]
            read = reads[member_index]
            member = reverse_complement(read.sequence) if reverse else read.sequence
            qualities = read.qualities
            if reverse and qualities is not None:
                qualities = tuple(reversed(qualities))
            selected.append((member, offset, qualities))
    return selected


def apply_damage_aware_consensus(
    reads: list[Read],
    *,
    molecule_ids: list[int] | None = None,
    report_path: Path | None = None,
    end_window: int = 5,
    anchor_k: int = 15,
    minimum_anchor_matches: int = 2,
    maximum_anchor_occurrences: int = 100,
    minimum_overlap: int = 30,
    minimum_identity: float = 0.90,
    minimum_ry_identity: float = 0.99,
    minimum_support: int = 3,
    dominance_ratio: float = 4.0,
    minimum_base_quality: int = 20,
    anchors_per_read: int = 8,
) -> tuple[list[Read], DamageConsensusSummary]:
    """Correct strongly supported terminal deaminations and otherwise abstain."""
    if not reads:
        return [], DamageConsensusSummary()
    if end_window < 1:
        raise ValueError("end_window must be at least 1")
    if anchor_k < 3:
        raise ValueError("anchor_k must be at least 3")
    if minimum_support < 2:
        raise ValueError("minimum_support must be at least 2")
    if anchors_per_read < minimum_anchor_matches:
        raise ValueError("anchors_per_read must cover minimum_anchor_matches")
    if dominance_ratio < 1:
        raise ValueError("dominance_ratio must be at least 1")
    if not 0 <= minimum_identity <= 1 or not 0 <= minimum_ry_identity <= 1:
        raise ValueError("identity thresholds must be between 0 and 1")
    if not 0 <= minimum_base_quality <= 254:
        raise ValueError("minimum_base_quality must be between 0 and 254")
    molecules = list(range(len(reads))) if molecule_ids is None else molecule_ids
    if len(molecules) != len(reads):
        raise ValueError("molecule IDs must align one-to-one with reads")

    position_bits = max(1, max(len(read.sequence) for read in reads).bit_length())
    anchors = _anchor_index(
        reads,
        anchor_k,
        end_window,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    corrected: list[Read] = []
    fieldnames = [
        "read_index",
        "read_name",
        "position",
        "observed_base",
        "proposed_base",
        "decision",
        "support_A",
        "support_C",
        "support_G",
        "support_T",
        "aligned_molecules",
    ]
    report_handle = None
    writer = None
    if report_path is not None:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_handle = report_path.open("w", encoding="utf-8", newline="")
        writer = csv.DictWriter(
            report_handle, fieldnames=fieldnames, delimiter="\t"
        )
        writer.writeheader()
    corrected_bases = 0
    corrected_reads = 0
    insufficient = 0
    ambiguous = 0
    candidate_bases = 0
    for read_index, read in enumerate(reads):
        positions = [
            (position, "C")
            for position in range(min(end_window, len(read.sequence)))
            if read.sequence[position] == "T"
        ] + [
            (position, "G")
            for position in range(
                max(0, len(read.sequence) - end_window), len(read.sequence)
            )
            if read.sequence[position] == "A"
        ]
        if not positions:
            corrected.append(read)
            continue
        alignments = _alignments(
            read_index,
            reads,
            molecules,
            anchors,
            anchor_k=anchor_k,
            end_window=end_window,
            minimum_anchor_matches=minimum_anchor_matches,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_overlap=minimum_overlap,
            minimum_identity=minimum_identity,
            minimum_ry_identity=minimum_ry_identity,
            anchors_per_read=anchors_per_read,
            position_bits=position_bits,
        )
        bases = list(read.sequence)
        changed = False
        for position, proposed in positions:
            candidate_bases += 1
            counts = {base: 0 for base in "ACGT"}
            for member, offset, qualities in alignments:
                member_position = position - offset
                if not end_window <= member_position < len(member) - end_window:
                    continue
                if (
                    qualities is None
                    or qualities[member_position] < minimum_base_quality
                ):
                    continue
                base = member[member_position]
                if base in counts:
                    counts[base] += 1
            runner_up = max(
                (count for base, count in counts.items() if base != proposed),
                default=0,
            )
            if counts[proposed] < minimum_support:
                decision = "insufficient_support"
                insufficient += 1
            elif runner_up and counts[proposed] <= dominance_ratio * runner_up:
                decision = "ambiguous"
                ambiguous += 1
            else:
                decision = "correct"
                bases[position] = proposed
                corrected_bases += 1
                changed = True
            if writer is not None:
                writer.writerow({
                    "read_index": read_index,
                    "read_name": read.name,
                    "position": position,
                    "observed_base": read.sequence[position],
                    "proposed_base": proposed,
                    "decision": decision,
                    "support_A": counts["A"],
                    "support_C": counts["C"],
                    "support_G": counts["G"],
                    "support_T": counts["T"],
                    "aligned_molecules": len(alignments),
                })
        if changed:
            corrected_reads += 1
            corrected.append(Read(read.name, "".join(bases), read.qualities))
        else:
            corrected.append(read)

    if report_handle is not None:
        report_handle.close()
    return corrected, DamageConsensusSummary(
        reads=len(reads),
        candidate_bases=candidate_bases,
        corrected_bases=corrected_bases,
        corrected_reads=corrected_reads,
        insufficient_support=insufficient,
        ambiguous=ambiguous,
    )
