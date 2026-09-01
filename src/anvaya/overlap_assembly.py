"""Bounded direct overlap-layout-consensus assembly for merged fragments."""

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from anvaya.damage_consensus import _anchor_index, _identity, _ry, _sketch_anchors
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


@dataclass(slots=True, frozen=True)
class OverlapAssemblySummary:
    """Audit counts for conservative direct overlap assembly."""

    input_reads: int = 0
    output_contigs: int = 0
    clusters: int = 0
    clustered_reads: int = 0
    extension_rounds: int = 0
    extended_contigs: int = 0
    added_bases: int = 0
    contig_iterations: int = 0
    contig_merges: int = 0
    ambiguous_extensions: int = 0
    corrected_bases: int = 0
    correction_candidates: int = 0
    correction_terminal_only: int = 0
    correction_insufficient_support: int = 0
    correction_ambiguous: int = 0


@dataclass(slots=True, frozen=True)
class OverlapCorrectionEvent:
    """One accepted correction retained for reference-assisted auditing."""

    cluster: int
    round: int
    center_name: str
    position: int
    original_base: str
    corrected_base: str
    internal_support: int
    competing_support: int
    sequence_before: str
    sequence_after: str


@dataclass(slots=True, frozen=True)
class _Alignment:
    read_index: int
    sequence: str
    offset: int
    support: int


@dataclass(slots=True, frozen=True)
class _Consensus:
    sequence: str
    left_extension: int
    right_extension: int
    corrected_bases: int
    correction_candidates: int
    correction_terminal_only: int
    correction_insufficient_support: int
    correction_ambiguous: int
    corrections: tuple[tuple[int, str, str, int, int], ...]


def _target_anchors(
    target: str,
    anchor_k: int,
    anchors_per_end: int,
    window: int,
) -> list[tuple[int, int, bool]]:
    """Sketch only contig ends, where an unused read can extend the layout."""
    if len(target) <= 2 * window:
        return _sketch_anchors(target, anchor_k, 0, anchors_per_end * 2)
    anchors = _sketch_anchors(target[:window], anchor_k, 0, anchors_per_end)
    suffix_start = len(target) - window
    anchors.extend(
        (anchor, position + suffix_start, reverse)
        for anchor, position, reverse in _sketch_anchors(
            target[suffix_start:], anchor_k, 0, anchors_per_end
        )
    )
    return anchors


def _candidate_alignments(
    target: str,
    reads: list[Read],
    molecule_ids: list[int],
    anchors: dict[int, list[int]],
    unavailable: list[bool] | set[int],
    *,
    anchor_k: int,
    anchors_per_read: int,
    maximum_anchor_occurrences: int,
    minimum_anchor_matches: int,
    minimum_overlap: int,
    minimum_identity: float,
    minimum_ry_identity: float,
    position_bits: int,
    target_window: int,
) -> list[_Alignment]:
    payload_bits = position_bits + 1
    position_mask = (1 << position_bits) - 1
    votes: dict[tuple[int, bool, int], int] = defaultdict(int)
    for anchor, target_position, target_reverse in _target_anchors(
        target, anchor_k, anchors_per_read, target_window
    ):
        occurrences = anchors.get(anchor, ())
        if not occurrences or len(occurrences) > maximum_anchor_occurrences:
            continue
        for packed in occurrences:
            member_index = packed >> payload_bits
            if (
                member_index in unavailable
                if isinstance(unavailable, set)
                else unavailable[member_index]
            ):
                continue
            member_position = (packed >> 1) & position_mask
            reverse = target_reverse != bool(packed & 1)
            if reverse:
                member_position = (
                    len(reads[member_index].sequence)
                    - member_position
                    - anchor_k
                )
            votes[(member_index, reverse, target_position - member_position)] += 1

    by_molecule: dict[int, list[_Alignment]] = defaultdict(list)
    for (member_index, reverse, offset), support in votes.items():
        if support < minimum_anchor_matches:
            continue
        member = reads[member_index].sequence
        if reverse:
            member = reverse_complement(member)
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
        by_molecule[molecule_ids[member_index]].append(
            _Alignment(member_index, member, offset, support)
        )

    selected: list[_Alignment] = []
    for candidates in by_molecule.values():
        best_support = max(candidate.support for candidate in candidates)
        best = [candidate for candidate in candidates if candidate.support == best_support]
        if len(best) == 1:
            selected.append(best[0])
    return selected


def _consensus(
    target: str,
    members: list[_Alignment],
    *,
    minimum_extension_support: int,
    dominance_ratio: float,
    minimum_correction_support: int = 3,
    correction_dominance_ratio: float = 2.0,
    damage_end_window: int = 0,
    allow_internal_consensus: bool = True,
) -> _Consensus:
    observations: dict[int, dict[str, int]] = defaultdict(
        lambda: {base: 0 for base in "ACGT"}
    )
    internal_observations: dict[int, dict[str, int]] = defaultdict(
        lambda: {base: 0 for base in "ACGT"}
    )
    for position, base in enumerate(target):
        if base in "ACGT":
            observations[position][base] += 1
    for member in members:
        for member_position, base in enumerate(member.sequence):
            if base in "ACGT":
                target_position = member.offset + member_position
                observations[target_position][base] += 1
                if (
                    member_position >= damage_end_window
                    and member_position < len(member.sequence) - damage_end_window
                ):
                    internal_observations[target_position][base] += 1

    accepted_extensions: dict[int, str] = {}
    accepted_internal: dict[int, str] = {}
    accepted_corrections: dict[int, str] = {}
    correction_evidence: dict[int, tuple[int, int]] = {}
    correction_candidates = correction_terminal_only = 0
    correction_insufficient_support = correction_ambiguous = 0
    for position, counts in observations.items():
        ranked = sorted(counts.items(), key=lambda item: item[1], reverse=True)
        base, support = ranked[0]
        runner_up = ranked[1][1]
        if position < 0 or position >= len(target):
            if support >= minimum_extension_support and (
                runner_up == 0 or support >= dominance_ratio * runner_up
            ):
                accepted_extensions[position] = base
            continue
        if allow_internal_consensus and base != target[position] and (
            support >= minimum_extension_support
            and (runner_up == 0 or support >= dominance_ratio * runner_up)
        ):
            accepted_internal[position] = base
        expected: str | None = None
        if position < damage_end_window and target[position] == "T":
            expected = "C"
        elif (
            position >= len(target) - damage_end_window
            and target[position] == "A"
        ):
            expected = "G"
        if expected is None or counts[expected] == 0:
            continue
        correction_candidates += 1
        correction_counts = internal_observations[position]
        correction_support = correction_counts[expected]
        if correction_support == 0:
            correction_terminal_only += 1
            continue
        if correction_support < minimum_correction_support:
            correction_insufficient_support += 1
            continue
        competing_support = max(
            (1 if base == target[position] else 0) + correction_counts[base]
            for base in "ACGT"
            if base != expected
        )
        if correction_support < correction_dominance_ratio * competing_support:
            correction_ambiguous += 1
            continue
        accepted_corrections[position] = expected
        correction_evidence[position] = (correction_support, competing_support)

    left = 0
    while left - 1 in accepted_extensions:
        left -= 1
    right = len(target)
    while right in accepted_extensions:
        right += 1
    sequence = list(target)
    for position, base in accepted_internal.items():
        sequence[position] = base
    for position, base in accepted_corrections.items():
        sequence[position] = base
    assembled = (
        "".join(accepted_extensions[position] for position in range(left, 0))
        + "".join(sequence)
        + "".join(
            accepted_extensions[position] for position in range(len(target), right)
        )
    )
    return _Consensus(
        assembled,
        -left,
        right - len(target),
        len(accepted_corrections),
        correction_candidates,
        correction_terminal_only,
        correction_insufficient_support,
        correction_ambiguous,
        tuple(
            (
                position,
                target[position],
                corrected_base,
                correction_evidence[position][0],
                correction_evidence[position][1],
            )
            for position, corrected_base in sorted(accepted_corrections.items())
        ),
    )


def write_overlap_correction_report(
    events: list[OverlapCorrectionEvent],
    path: str | Path,
) -> None:
    """Write accepted overlap corrections for an external truth audit."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [field.name for field in OverlapCorrectionEvent.__dataclass_fields__.values()]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for event in events:
            writer.writerow({field: getattr(event, field) for field in fields})


def _unique_best_extensions(
    target: str,
    candidates: list[_Alignment],
    *,
    minimum_overlap_margin: int,
) -> tuple[list[_Alignment], set[int], int]:
    """Select at most one clearly best proper extension on each side."""
    sides: dict[str, list[tuple[tuple[int, int, int], _Alignment]]] = {
        "left": [],
        "right": [],
    }
    contained: set[int] = set()
    for candidate in candidates:
        start = max(0, candidate.offset)
        stop = min(len(target), candidate.offset + len(candidate.sequence))
        overlap = stop - start
        left_extension = max(0, -candidate.offset)
        right_extension = max(
            0, candidate.offset + len(candidate.sequence) - len(target)
        )
        if not left_extension and not right_extension:
            contained.add(candidate.read_index)
            continue
        if left_extension:
            sides["left"].append(
                ((overlap, candidate.support, left_extension), candidate)
            )
        if right_extension:
            sides["right"].append(
                ((overlap, candidate.support, right_extension), candidate)
            )

    selected: dict[int, _Alignment] = {}
    ambiguous = 0
    for ranked in sides.values():
        ranked.sort(key=lambda item: item[0], reverse=True)
        if not ranked:
            continue
        best_score, best = ranked[0]
        if len(ranked) > 1:
            second_score = ranked[1][0]
            if best_score[0] - second_score[0] < minimum_overlap_margin:
                ambiguous += 1
                continue
        selected[best.read_index] = best
    return list(selected.values()), contained, ambiguous


def _merge_contig_pass(
    contigs: list[Read],
    *,
    anchor_k: int,
    anchors_per_read: int,
    maximum_anchor_occurrences: int,
    minimum_anchor_matches: int,
    minimum_overlap: int,
    minimum_identity: float,
    minimum_ry_identity: float,
    maximum_rounds: int,
    minimum_overlap_margin: int,
) -> tuple[list[Read], int, int, int, int]:
    """Merge unambiguous contig overlaps while preserving every other contig."""
    if not contigs:
        return [], 0, 0, 0, 0
    position_bits = max(
        1, max(len(contig.sequence) for contig in contigs).bit_length()
    )
    target_window = max(len(contig.sequence) for contig in contigs)
    anchors = _anchor_index(
        contigs, anchor_k, 0, anchors_per_read, maximum_anchor_occurrences
    )
    molecule_ids = list(range(len(contigs)))
    assigned = [False] * len(contigs)
    output: list[Read] = []
    merges = rounds = added_bases = ambiguous = 0
    order = sorted(
        range(len(contigs)),
        key=lambda index: (-len(contigs[index].sequence), index),
    )
    for center_index in order:
        if assigned[center_index]:
            continue
        assigned[center_index] = True
        target = contigs[center_index].sequence
        seed_length = len(target)
        for _ in range(maximum_rounds):
            candidates = _candidate_alignments(
                target,
                contigs,
                molecule_ids,
                anchors,
                assigned,
                anchor_k=anchor_k,
                anchors_per_read=anchors_per_read,
                maximum_anchor_occurrences=maximum_anchor_occurrences,
                minimum_anchor_matches=minimum_anchor_matches,
                minimum_overlap=minimum_overlap,
                minimum_identity=minimum_identity,
                minimum_ry_identity=minimum_ry_identity,
                position_bits=position_bits,
                target_window=target_window,
            )
            selected, contained, branch_count = _unique_best_extensions(
                target,
                candidates,
                minimum_overlap_margin=minimum_overlap_margin,
            )
            ambiguous += branch_count
            for read_index in contained:
                assigned[read_index] = True
            if not selected:
                break
            for candidate in selected:
                assigned[candidate.read_index] = True
            result = _consensus(
                target,
                selected,
                minimum_extension_support=1,
                dominance_ratio=4.0,
            )
            if result.sequence == target:
                break
            target = result.sequence
            merges += len(selected)
            rounds += 1
        added_bases += len(target) - seed_length
        output.append(Read(f"merged_contig_{len(output) + 1}", target))
    return output, merges, rounds, added_bases, ambiguous


def _polish_frozen_contigs(
    contigs: list[Read],
    reads: list[Read],
    molecule_ids: list[int],
    anchors: dict[int, list[int]],
    *,
    anchor_k: int,
    anchors_per_read: int,
    maximum_anchor_occurrences: int,
    minimum_anchor_matches: int,
    minimum_overlap: int,
    minimum_identity: float,
    minimum_ry_identity: float,
    position_bits: int,
    target_window: int,
    minimum_correction_support: int,
    correction_dominance_ratio: float,
    damage_end_window: int,
    correction_events: list[OverlapCorrectionEvent] | None,
) -> tuple[list[Read], int, int, int, int, int]:
    """Polish final contig ends without changing frozen layout or membership."""
    polished: list[Read] = []
    corrected_bases = correction_candidates = correction_terminal_only = 0
    correction_insufficient_support = correction_ambiguous = 0
    for contig_index, contig in enumerate(contigs):
        candidates = _candidate_alignments(
            contig.sequence,
            reads,
            molecule_ids,
            anchors,
            set(),
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=minimum_overlap,
            minimum_identity=minimum_identity,
            minimum_ry_identity=minimum_ry_identity,
            position_bits=position_bits,
            target_window=target_window,
        )
        result = _consensus(
            contig.sequence,
            candidates,
            minimum_extension_support=len(candidates) + 2,
            dominance_ratio=4.0,
            minimum_correction_support=minimum_correction_support,
            correction_dominance_ratio=correction_dominance_ratio,
            damage_end_window=damage_end_window,
            allow_internal_consensus=False,
        )
        corrected_bases += result.corrected_bases
        correction_candidates += result.correction_candidates
        correction_terminal_only += result.correction_terminal_only
        correction_insufficient_support += result.correction_insufficient_support
        correction_ambiguous += result.correction_ambiguous
        if correction_events is not None:
            for (
                position,
                original_base,
                corrected_base,
                internal_support,
                competing_support,
            ) in result.corrections:
                corrected_sequence = list(contig.sequence)
                corrected_sequence[position] = corrected_base
                correction_events.append(OverlapCorrectionEvent(
                    cluster=contig_index + 1,
                    round=1,
                    center_name=contig.name,
                    position=position,
                    original_base=original_base,
                    corrected_base=corrected_base,
                    internal_support=internal_support,
                    competing_support=competing_support,
                    sequence_before=contig.sequence,
                    sequence_after="".join(corrected_sequence),
                ))
        polished.append(Read(contig.name, result.sequence))
    return (
        polished,
        corrected_bases,
        correction_candidates,
        correction_terminal_only,
        correction_insufficient_support,
        correction_ambiguous,
    )


def assemble_overlap_contigs(
    reads: list[Read],
    *,
    molecule_ids: list[int] | None = None,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
    minimum_identity: float = 0.90,
    minimum_ry_identity: float = 0.99,
    minimum_cluster_size: int = 5,
    minimum_consensus_support: int = 5,
    minimum_correction_support: int = 3,
    dominance_ratio: float = 4.0,
    correction_dominance_ratio: float = 2.0,
    damage_end_window: int = 5,
    maximum_rounds: int = 3,
    maximum_contig_iterations: int = 3,
    minimum_overlap_margin: int = 3,
    minimum_output_length: int = 0,
    correction_events: list[OverlapCorrectionEvent] | None = None,
) -> tuple[list[Read], OverlapAssemblySummary]:
    """Build contigs directly by bounded iterative overlap extension."""
    if not reads:
        return [], OverlapAssemblySummary()
    if minimum_cluster_size < 2:
        raise ValueError("minimum_cluster_size must be at least 2")
    if minimum_consensus_support < 2:
        raise ValueError("minimum_consensus_support must be at least 2")
    if minimum_correction_support < 2:
        raise ValueError("minimum_correction_support must be at least 2")
    if damage_end_window < 0:
        raise ValueError("damage_end_window must not be negative")
    if maximum_rounds < 1:
        raise ValueError("maximum_rounds must be at least 1")
    if maximum_contig_iterations < 0:
        raise ValueError("maximum_contig_iterations must not be negative")
    if minimum_overlap_margin < 1:
        raise ValueError("minimum_overlap_margin must be at least 1")
    if minimum_output_length < 0:
        raise ValueError("minimum_output_length must not be negative")
    if anchors_per_read < minimum_anchor_matches:
        raise ValueError("anchors_per_read must cover minimum_anchor_matches")
    molecules = list(range(len(reads))) if molecule_ids is None else molecule_ids
    if len(molecules) != len(reads):
        raise ValueError("molecule IDs must align one-to-one with reads")

    position_bits = max(1, max(len(read.sequence) for read in reads).bit_length())
    target_window = max(len(read.sequence) for read in reads)
    anchors = _anchor_index(
        reads, anchor_k, 0, anchors_per_read, maximum_anchor_occurrences
    )
    claimed = [False] * len(reads)
    contigs: list[Read] = []
    clusters = clustered_reads = extension_rounds = 0
    extended_contigs = added_bases = corrected_bases = 0
    correction_candidates = correction_terminal_only = 0
    correction_insufficient_support = correction_ambiguous = 0
    order = sorted(
        range(len(reads)),
        key=lambda index: (-len(reads[index].sequence), index),
    )
    for center_index in order:
        if claimed[center_index]:
            continue
        target = reads[center_index].sequence
        seed_length = len(target)
        claimed[center_index] = True
        unavailable = {center_index}
        members: list[_Alignment] = []
        cluster_indices = {center_index}
        rounds = 0
        for _ in range(maximum_rounds):
            candidates = _candidate_alignments(
                target,
                reads,
                molecules,
                anchors,
                unavailable,
                anchor_k=anchor_k,
                anchors_per_read=anchors_per_read,
                maximum_anchor_occurrences=maximum_anchor_occurrences,
                minimum_anchor_matches=minimum_anchor_matches,
                minimum_overlap=minimum_overlap,
                minimum_identity=minimum_identity,
                minimum_ry_identity=minimum_ry_identity,
                position_bits=position_bits,
                target_window=target_window,
            )
            if not candidates:
                break
            if not members and len(candidates) + 1 < minimum_cluster_size:
                break
            tentative_members = members + candidates
            result = _consensus(
                target,
                tentative_members,
                minimum_extension_support=minimum_consensus_support,
                dominance_ratio=dominance_ratio,
                minimum_correction_support=minimum_correction_support,
                correction_dominance_ratio=correction_dominance_ratio,
                damage_end_window=0,
            )
            rounds += 1
            changed = result.sequence != target
            if changed or not members:
                for candidate in candidates:
                    claimed[candidate.read_index] = True
                    unavailable.add(candidate.read_index)
                    cluster_indices.add(candidate.read_index)
                members = tentative_members
            if not changed:
                break
            if result.left_extension:
                members = [
                    _Alignment(
                        member.read_index,
                        member.sequence,
                        member.offset + result.left_extension,
                        member.support,
                    )
                    for member in members
                ]
            target = result.sequence
        if len(cluster_indices) < minimum_cluster_size:
            continue
        clusters += 1
        clustered_reads += len(cluster_indices)
        extension_rounds += rounds
        extension = len(target) - seed_length
        if extension:
            extended_contigs += 1
            added_bases += extension
        contigs.append(Read(f"overlap_contig_{clusters}", target))

    contig_iterations = contig_merges = ambiguous_extensions = 0
    for _ in range(maximum_contig_iterations):
        merged, merges, rounds, added, ambiguous = _merge_contig_pass(
            contigs,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=minimum_overlap,
            minimum_identity=minimum_identity,
            minimum_ry_identity=minimum_ry_identity,
            maximum_rounds=maximum_rounds,
            minimum_overlap_margin=minimum_overlap_margin,
        )
        ambiguous_extensions += ambiguous
        if not merges:
            break
        contigs = merged
        contig_iterations += 1
        contig_merges += merges
        extension_rounds += rounds
        added_bases += added

    if damage_end_window:
        (
            contigs,
            corrected_bases,
            correction_candidates,
            correction_terminal_only,
            correction_insufficient_support,
            correction_ambiguous,
        ) = _polish_frozen_contigs(
            contigs,
            reads,
            molecules,
            anchors,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=minimum_overlap,
            minimum_identity=minimum_identity,
            minimum_ry_identity=minimum_ry_identity,
            position_bits=position_bits,
            target_window=target_window,
            minimum_correction_support=minimum_correction_support,
            correction_dominance_ratio=correction_dominance_ratio,
            damage_end_window=damage_end_window,
            correction_events=correction_events,
        )

    contigs = [
        contig for contig in contigs if len(contig.sequence) >= minimum_output_length
    ]

    return contigs, OverlapAssemblySummary(
        input_reads=len(reads),
        output_contigs=len(contigs),
        clusters=clusters,
        clustered_reads=clustered_reads,
        extension_rounds=extension_rounds,
        extended_contigs=extended_contigs,
        added_bases=added_bases,
        contig_iterations=contig_iterations,
        contig_merges=contig_merges,
        ambiguous_extensions=ambiguous_extensions,
        corrected_bases=corrected_bases,
        correction_candidates=correction_candidates,
        correction_terminal_only=correction_terminal_only,
        correction_insufficient_support=correction_insufficient_support,
        correction_ambiguous=correction_ambiguous,
    )
