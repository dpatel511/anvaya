"""Bounded direct overlap-layout-consensus assembly for merged fragments."""

import csv
from collections import defaultdict
from dataclasses import dataclass, replace
from math import sqrt
from pathlib import Path

from anvaya.damage_consensus import _anchor_index, _identity, _ry, _sketch_anchors
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


_DAMAGE_PAIRS = {frozenset(("C", "T")), frozenset(("G", "A"))}


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
    candidate_offsets: int = 0
    candidate_below_anchor_support: int = 0
    candidate_short_overlap: int = 0
    candidate_dna_rejected: int = 0
    candidate_ry_rejected: int = 0
    candidate_molecule_ambiguous: int = 0
    candidate_alignments: int = 0
    unavailable_anchor_hits: int = 0
    targets_without_candidates: int = 0
    clusters_below_minimum_size: int = 0
    consensus_without_extension: int = 0
    reciprocal_extension_checks: int = 0
    reciprocal_extension_rejections: int = 0
    confidence_ranked_sides: int = 0
    confidence_changed_winners: int = 0
    confidence_ambiguous_extensions: int = 0
    deferred_reads: int = 0
    cross_cluster_reads: int = 0
    deferred_candidate_reads: int = 0
    cross_cluster_candidate_reads: int = 0
    deferred_extension_candidates: int = 0
    cross_cluster_extension_candidates: int = 0
    deferred_ambiguous_assignments: int = 0
    cross_cluster_ambiguous_assignments: int = 0
    deferred_unique_assignments: int = 0
    cross_cluster_unique_assignments: int = 0
    deferred_reciprocal_assignments: int = 0
    cross_cluster_reciprocal_assignments: int = 0
    recruitment_supported_contigs: int = 0
    recruitment_supported_bases: int = 0
    deferred_supported_contigs: int = 0
    deferred_supported_bases: int = 0
    cross_cluster_supported_contigs: int = 0
    cross_cluster_supported_bases: int = 0
    contig_link_candidate_overlaps: int = 0
    contig_link_unique_overlaps: int = 0
    contig_link_reciprocal_overlaps: int = 0
    contig_link_read_supported_overlaps: int = 0
    contig_link_support_at_least_1: int = 0
    contig_link_support_at_least_2: int = 0
    contig_link_support_at_least_3: int = 0
    contig_link_support_at_least_5: int = 0
    contig_link_max_read_support: int = 0
    contig_link_ambiguous_ends: int = 0
    contig_link_cyclic_components: int = 0
    contig_link_linear_chains: int = 0
    contig_link_projected_joins: int = 0
    contig_link_projected_merged_bases: int = 0
    contig_link_projected_longest_contig: int = 0


@dataclass(slots=True)
class _CandidateDiagnostics:
    offsets: int = 0
    below_anchor_support: int = 0
    short_overlap: int = 0
    dna_rejected: int = 0
    ry_rejected: int = 0
    molecule_ambiguous: int = 0
    alignments: int = 0
    unavailable_anchor_hits: int = 0


@dataclass(slots=True)
class _ReciprocalDiagnostics:
    checks: int = 0
    rejections: int = 0


@dataclass(slots=True)
class _RankingDiagnostics:
    ranked_sides: int = 0
    changed_winners: int = 0
    ambiguous_extensions: int = 0


@dataclass(slots=True)
class _RecruitmentDiagnostics:
    deferred_reads: int = 0
    cross_cluster_reads: int = 0
    deferred_candidate_reads: int = 0
    cross_cluster_candidate_reads: int = 0
    deferred_extension_candidates: int = 0
    cross_cluster_extension_candidates: int = 0
    deferred_ambiguous_assignments: int = 0
    cross_cluster_ambiguous_assignments: int = 0
    deferred_unique_assignments: int = 0
    cross_cluster_unique_assignments: int = 0
    deferred_reciprocal_assignments: int = 0
    cross_cluster_reciprocal_assignments: int = 0
    supported_contigs: int = 0
    supported_bases: int = 0
    deferred_supported_contigs: int = 0
    deferred_supported_bases: int = 0
    cross_cluster_supported_contigs: int = 0
    cross_cluster_supported_bases: int = 0


@dataclass(slots=True)
class _ContigLinkDiagnostics:
    candidate_overlaps: int = 0
    unique_overlaps: int = 0
    reciprocal_overlaps: int = 0
    read_supported_overlaps: int = 0
    support_at_least_1: int = 0
    support_at_least_2: int = 0
    support_at_least_3: int = 0
    support_at_least_5: int = 0
    max_read_support: int = 0
    ambiguous_ends: int = 0
    cyclic_components: int = 0
    linear_chains: int = 0
    projected_joins: int = 0
    projected_merged_bases: int = 0
    projected_longest_contig: int = 0


@dataclass(slots=True, frozen=True)
class TwoTierRedundancyDiagnostics:
    """Projected support-three rescue after conservative redundancy removal."""

    rescue_contigs: int = 0
    rescue_bases: int = 0
    contained_by_primary: int = 0
    redundant_with_rescue: int = 0
    extends_primary: int = 0
    ambiguous_primary_extensions: int = 0
    replacement_contigs: int = 0
    replacement_added_bases: int = 0
    replacement_projected_bases: int = 0
    replacement_projected_n50: int = 0
    replacement_projected_longest_contig: int = 0
    novel_contigs: int = 0
    novel_bases: int = 0
    projected_contigs: int = 0
    projected_bases: int = 0
    projected_n50: int = 0
    projected_longest_contig: int = 0


@dataclass(slots=True, frozen=True)
class IterativeReclusteringRound:
    """One output-neutral global reclustering round."""

    iteration: int
    extended_sequences: int
    added_bases: int
    candidate_extensions: int
    reused_candidates: int
    conflicting_candidates: int
    reciprocal_checks: int
    reciprocal_rejections: int
    derived_overlap_rejections: int
    raw_confirmed_mismatch_overlaps: int
    raw_support_rejections: int
    assembled_sequences: int
    redundant_sequences: int
    projected_contigs: int
    projected_bases: int
    projected_n50: int
    projected_longest_contig: int


@dataclass(slots=True, frozen=True)
class IterativeReclusteringDiagnostics:
    """Projection from CarpeDeam-style fixed-pool iterative reclustering."""

    iterations: int = 0
    converged: bool = False
    rollback_rejected: bool = False
    extended_sequences: int = 0
    added_bases: int = 0
    candidate_extensions: int = 0
    reused_candidates: int = 0
    conflicting_candidates: int = 0
    reciprocal_checks: int = 0
    reciprocal_rejections: int = 0
    derived_overlap_rejections: int = 0
    raw_confirmed_mismatch_overlaps: int = 0
    raw_support_rejections: int = 0
    assembled_sequences: int = 0
    redundant_sequences: int = 0
    projected_contigs: int = 0
    projected_bases: int = 0
    projected_n50: int = 0
    projected_longest_contig: int = 0
    rounds: tuple[IterativeReclusteringRound, ...] = ()


@dataclass(slots=True, frozen=True)
class MasterOverlapGraphDiagnostics:
    """Report-only exact string-graph projection of assembled contigs."""

    input_contigs: int = 0
    candidate_dovetails: int = 0
    exact_dovetails: int = 0
    transitive_edges_removed: int = 0
    ambiguous_ends: int = 0
    reciprocal_edges: int = 0
    linear_paths: int = 0
    cyclic_components: int = 0
    merged_contigs: int = 0
    added_bases: int = 0
    projected_contigs: int = 0
    projected_bases: int = 0
    projected_n50: int = 0
    projected_longest_contig: int = 0


@dataclass(slots=True, frozen=True)
class RawConfirmedMasterGraphDiagnostics:
    """Report-only master graph admitting raw-confirmed near-exact edges."""

    input_contigs: int = 0
    candidate_dovetails: int = 0
    exact_dovetails: int = 0
    near_exact_candidates: int = 0
    near_exact_accepted: int = 0
    mismatch_positions: int = 0
    insufficient_raw_support: int = 0
    strain_conflicts: int = 0
    exact_preferred: int = 0
    transitive_edges_removed: int = 0
    ambiguous_ends: int = 0
    reciprocal_edges: int = 0
    linear_paths: int = 0
    cyclic_components: int = 0
    merged_contigs: int = 0
    corrected_overlap_bases: int = 0
    added_bases: int = 0
    projected_contigs: int = 0
    projected_bases: int = 0
    projected_n50: int = 0
    projected_longest_contig: int = 0


@dataclass(slots=True, frozen=True)
class StrainSafeContainmentDiagnostics:
    """Report-only removal of contained contigs lacking raw allele support."""

    input_contigs: int = 0
    candidate_containments: int = 0
    exact_duplicates: int = 0
    near_duplicate_candidates: int = 0
    strain_protected: int = 0
    insufficient_evidence: int = 0
    unsupported_variants: int = 0
    removed_contigs: int = 0
    removed_bases: int = 0
    projected_contigs: int = 0
    projected_bases: int = 0
    projected_n50: int = 0
    projected_longest_contig: int = 0


@dataclass(slots=True, frozen=True)
class _SupportedContigLink:
    first: int
    second: int
    first_end: str
    second_end: str
    overlap: int
    support: int


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
class _MasterOverlapEdge:
    source: tuple[int, bool]
    target: tuple[int, bool]
    shift: int
    overlap: int
    corrections: tuple[tuple[int, str], ...] = ()


@dataclass(slots=True, frozen=True)
class _MasterGraphProjection:
    transitive_edges_removed: int = 0
    ambiguous_ends: int = 0
    reciprocal_edges: int = 0
    linear_paths: int = 0
    cyclic_components: int = 0
    merged_contigs: int = 0
    corrected_overlap_bases: int = 0
    added_bases: int = 0


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
    diagnostics: _CandidateDiagnostics | None = None,
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
                if diagnostics is not None:
                    diagnostics.unavailable_anchor_hits += 1
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
        if diagnostics is not None:
            diagnostics.offsets += 1
        if support < minimum_anchor_matches:
            if diagnostics is not None:
                diagnostics.below_anchor_support += 1
            continue
        member = reads[member_index].sequence
        if reverse:
            member = reverse_complement(member)
        start = max(0, offset)
        stop = min(len(target), offset + len(member))
        if stop - start < minimum_overlap:
            if diagnostics is not None:
                diagnostics.short_overlap += 1
            continue
        member_start = start - offset
        target_overlap = target[start:stop]
        member_overlap = member[member_start : member_start + stop - start]
        if _identity(target_overlap, member_overlap) < minimum_identity:
            if diagnostics is not None:
                diagnostics.dna_rejected += 1
            continue
        if _identity(_ry(target_overlap), _ry(member_overlap)) < minimum_ry_identity:
            if diagnostics is not None:
                diagnostics.ry_rejected += 1
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
            if diagnostics is not None:
                diagnostics.alignments += 1
        elif diagnostics is not None:
            diagnostics.molecule_ambiguous += 1
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
    minimum_internal_support: int | None = None,
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
        internal_support = (
            minimum_extension_support
            if minimum_internal_support is None
            else minimum_internal_support
        )
        if allow_internal_consensus and base != target[position] and (
            support >= internal_support
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
    damage_aware_ranking: bool = False,
    minimum_confidence_margin: float = 0.0,
    damage_mismatch_penalty: float = 0.25,
    ranking_damage_end_window: int = 5,
    diagnostics: _RankingDiagnostics | None = None,
) -> tuple[list[_Alignment], set[int], int]:
    """Select at most one clearly best proper extension on each side."""
    sides: dict[str, list[tuple[tuple[float, int, int, int], _Alignment]]] = {
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
        confidence = None
        if damage_aware_ranking:
            confidence = _overlap_confidence(
                target,
                candidate,
                damage_mismatch_penalty=damage_mismatch_penalty,
                damage_end_window=ranking_damage_end_window,
            )
        if left_extension:
            score = (
                (confidence, overlap, candidate.support, left_extension)
                if confidence is not None
                else (float(overlap), candidate.support, left_extension, 0)
            )
            sides["left"].append((score, candidate))
        if right_extension:
            score = (
                (confidence, overlap, candidate.support, right_extension)
                if confidence is not None
                else (float(overlap), candidate.support, right_extension, 0)
            )
            sides["right"].append((score, candidate))

    selected: dict[int, _Alignment] = {}
    ambiguous = 0
    for ranked in sides.values():
        if damage_aware_ranking and ranked and diagnostics is not None:
            diagnostics.ranked_sides += 1
            legacy_best = max(
                ranked,
                key=lambda item: (item[0][1], item[0][2], item[0][3]),
            )[1]
        else:
            legacy_best = None
        ranked.sort(key=lambda item: item[0], reverse=True)
        if not ranked:
            continue
        best_score, best = ranked[0]
        if (
            diagnostics is not None
            and legacy_best is not None
            and best.read_index != legacy_best.read_index
        ):
            diagnostics.changed_winners += 1
        if len(ranked) > 1:
            second_score = ranked[1][0]
            required_margin = (
                minimum_confidence_margin
                if damage_aware_ranking
                else minimum_overlap_margin
            )
            if best_score[0] - second_score[0] < required_margin:
                ambiguous += 1
                if damage_aware_ranking and diagnostics is not None:
                    diagnostics.ambiguous_extensions += 1
                continue
        selected[best.read_index] = best
    return list(selected.values()), contained, ambiguous


def _overlap_confidence(
    target: str,
    candidate: _Alignment,
    *,
    damage_mismatch_penalty: float,
    damage_end_window: int,
) -> float:
    """Return a conservative posterior match-rate score for one overlap."""
    start = max(0, candidate.offset)
    stop = min(len(target), candidate.offset + len(candidate.sequence))
    effective_errors = 0.0
    for target_position in range(start, stop):
        candidate_position = target_position - candidate.offset
        target_base = target[target_position]
        candidate_base = candidate.sequence[candidate_position]
        if target_base == candidate_base:
            continue
        terminal = (
            target_position < damage_end_window
            or target_position >= len(target) - damage_end_window
            or candidate_position < damage_end_window
            or candidate_position >= len(candidate.sequence) - damage_end_window
        )
        if terminal and frozenset((target_base, candidate_base)) in _DAMAGE_PAIRS:
            effective_errors += damage_mismatch_penalty
        else:
            effective_errors += 1.0
    overlap = stop - start
    effective_matches = overlap - effective_errors
    alpha = effective_matches + 1.0
    beta = effective_errors + 1.0
    total = alpha + beta
    mean = alpha / total
    variance = alpha * beta / (total * total * (total + 1.0))
    return mean - sqrt(variance)


def _filter_extension_boundary_support(
    target: str,
    candidates: list[_Alignment],
    selected: list[_Alignment],
    minimum_support: int,
) -> list[_Alignment]:
    """Keep extensions whose first appended base has independent support."""
    if minimum_support <= 1:
        return selected
    supported: list[_Alignment] = []
    for choice in selected:
        position = -1 if choice.offset < 0 else len(target)
        choice_base = choice.sequence[position - choice.offset]
        support = sum(
            candidate.offset <= position < candidate.offset + len(candidate.sequence)
            and candidate.sequence[position - candidate.offset] == choice_base
            for candidate in candidates
        )
        if support >= minimum_support:
            supported.append(choice)
    return supported


def _filter_raw_supported_extensions(
    target: str,
    candidates: list[_Alignment],
    selected: list[_Alignment],
    molecule_ids: list[int],
    derived: list[bool],
    minimum_support: int,
) -> tuple[list[_Alignment], int]:
    """Require appended bases to be corroborated by pristine source molecules."""
    supported: list[_Alignment] = []
    rejected = 0
    for choice in selected:
        position = -1 if choice.offset < 0 else len(target)
        choice_base = choice.sequence[position - choice.offset]
        support = {
            molecule_ids[candidate.read_index]
            for candidate in candidates
            if not derived[candidate.read_index]
            and candidate.offset <= position < candidate.offset + len(candidate.sequence)
            and candidate.sequence[position - candidate.offset] == choice_base
        }
        if len(support) >= minimum_support:
            supported.append(choice)
        else:
            rejected += 1
    return supported, rejected


def _raw_supports_target_mismatches(
    target: str,
    candidate: _Alignment,
    raw_candidates: list[_Alignment],
    raw_reads: list[Read],
    molecule_ids: list[int],
    minimum_support: int,
    minimum_margin: int,
    minimum_base_quality: int,
) -> bool:
    """Confirm one derived mismatch using quality-filtered raw molecules."""
    start = max(0, candidate.offset)
    stop = min(len(target), candidate.offset + len(candidate.sequence))
    mismatches = [
        position
        for position in range(start, stop)
        if target[position] != candidate.sequence[position - candidate.offset]
    ]
    if len(mismatches) != 1:
        return False
    position = mismatches[0]
    support_by_base: dict[str, set[int]] = defaultdict(set)
    for raw in raw_candidates:
        if raw.offset <= position < raw.offset + len(raw.sequence):
            raw_position = position - raw.offset
            record = raw_reads[raw.read_index]
            if record.qualities is not None:
                quality_position = (
                    raw_position
                    if raw.sequence == record.sequence
                    else len(record.sequence) - raw_position - 1
                )
                if record.qualities[quality_position] < minimum_base_quality:
                    continue
            support_by_base[raw.sequence[position - raw.offset]].add(
                molecule_ids[raw.read_index]
            )
    target_support = len(support_by_base[target[position]])
    competing_support = max(
        (
            len(support)
            for base, support in support_by_base.items()
            if base != target[position]
        ),
        default=0,
    )
    return (
        target_support >= minimum_support
        and target_support - competing_support >= minimum_margin
    )


def _projection_regresses(
    previous_n50: int,
    previous_bases: int,
    current_lengths: list[int],
) -> bool:
    """Reject a round that loses sequence without improving continuity."""
    return (
        _n50(current_lengths) <= previous_n50
        and sum(current_lengths) < previous_bases
    )


def _filter_reciprocal_best_extensions(
    target: str,
    selected: list[_Alignment],
    cluster_indices: set[int],
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
    minimum_overlap_margin: int,
    diagnostics: _ReciprocalDiagnostics,
) -> list[_Alignment]:
    """Keep extensions whose opposite end uniquely points to this cluster."""
    def side_supported(
        reciprocal_candidates: list[_Alignment],
        side: str,
    ) -> bool:
        ranked: list[tuple[int, _Alignment, str]] = []
        for candidate in reciprocal_candidates:
            start = max(0, candidate.offset)
            stop = min(len(choice.sequence), candidate.offset + len(candidate.sequence))
            overlap = stop - start
            if side == "left" and candidate.offset < 0:
                extension = candidate.sequence[:-candidate.offset]
            elif (
                side == "right"
                and candidate.offset + len(candidate.sequence) > len(choice.sequence)
            ):
                extension = candidate.sequence[len(choice.sequence) - candidate.offset :]
            else:
                continue
            ranked.append((overlap, candidate, extension))
        if not ranked:
            return False
        best_overlap = max(item[0] for item in ranked)
        contenders = [
            item
            for item in ranked
            if best_overlap - item[0] < minimum_overlap_margin
        ]
        cluster_extensions = [
            extension
            for _, candidate, extension in contenders
            if candidate.read_index in cluster_indices
        ]
        for cluster_extension in cluster_extensions:
            compatible = True
            for _, _, extension in contenders:
                shared = min(len(cluster_extension), len(extension))
                if side == "left":
                    first = cluster_extension[-shared:]
                    second = extension[-shared:]
                else:
                    first = cluster_extension[:shared]
                    second = extension[:shared]
                if (
                    _identity(first, second) < minimum_identity
                    or _identity(_ry(first), _ry(second)) < minimum_ry_identity
                ):
                    compatible = False
                    break
            if compatible:
                return True
        return False

    supported: list[_Alignment] = []
    for choice in selected:
        diagnostics.checks += 1
        reciprocal_candidates = _candidate_alignments(
            choice.sequence,
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
        needs_left = choice.offset + len(choice.sequence) > len(target)
        needs_right = choice.offset < 0
        has_left = side_supported(reciprocal_candidates, "left")
        has_right = side_supported(reciprocal_candidates, "right")
        if (not needs_left or has_left) and (not needs_right or has_right):
            supported.append(choice)
        else:
            diagnostics.rejections += 1
    return supported


def _audit_cross_cluster_recruitment(
    contigs: list[Read],
    cluster_indices: list[set[int]],
    reads: list[Read],
    read_clusters: list[int | None],
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
    minimum_overlap_margin: int,
    minimum_consensus_support: int,
    dominance_ratio: float,
    damage_aware_ranking: bool,
    minimum_confidence_margin: float,
    damage_mismatch_penalty: float,
    ranking_damage_end_window: int,
) -> _RecruitmentDiagnostics:
    """Measure safe cross-cluster end recruitment without changing contigs."""
    diagnostics = _RecruitmentDiagnostics(
        deferred_reads=sum(cluster is None for cluster in read_clusters),
        cross_cluster_reads=(
            sum(cluster is not None for cluster in read_clusters)
            if len(contigs) > 1
            else 0
        ),
    )
    deferred_candidate_reads: set[int] = set()
    cross_cluster_candidate_reads: set[int] = set()
    deferred_extension_reads: set[int] = set()
    cross_cluster_extension_reads: set[int] = set()
    placements: dict[int, list[tuple[tuple[float, int, int, int], int, _Alignment]]] = (
        defaultdict(list)
    )
    for contig_index, contig in enumerate(contigs):
        own_cluster = cluster_indices[contig_index]
        candidates = _candidate_alignments(
            contig.sequence,
            reads,
            molecule_ids,
            anchors,
            own_cluster,
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
        for candidate in candidates:
            candidate_cluster = read_clusters[candidate.read_index]
            candidate_read_set = (
                deferred_candidate_reads
                if candidate_cluster is None
                else cross_cluster_candidate_reads
            )
            candidate_read_set.add(candidate.read_index)
            start = max(0, candidate.offset)
            stop = min(
                len(contig.sequence),
                candidate.offset + len(candidate.sequence),
            )
            overlap = stop - start
            extension = max(
                0,
                -candidate.offset,
                candidate.offset + len(candidate.sequence) - len(contig.sequence),
            )
            if not extension:
                continue
            extension_read_set = (
                deferred_extension_reads
                if candidate_cluster is None
                else cross_cluster_extension_reads
            )
            extension_read_set.add(candidate.read_index)
            primary_score = (
                _overlap_confidence(
                    contig.sequence,
                    candidate,
                    damage_mismatch_penalty=damage_mismatch_penalty,
                    damage_end_window=ranking_damage_end_window,
                )
                if damage_aware_ranking
                else float(overlap)
            )
            placements[candidate.read_index].append(
                (
                    (primary_score, overlap, candidate.support, extension),
                    contig_index,
                    candidate,
                )
            )

    diagnostics.deferred_candidate_reads = len(deferred_candidate_reads)
    diagnostics.cross_cluster_candidate_reads = len(cross_cluster_candidate_reads)
    diagnostics.deferred_extension_candidates = len(deferred_extension_reads)
    diagnostics.cross_cluster_extension_candidates = len(cross_cluster_extension_reads)
    selected_by_contig: dict[int, list[_Alignment]] = defaultdict(list)
    for read_placements in placements.values():
        best_by_contig: dict[int, tuple[tuple[float, int, int, int], int, _Alignment]] = {}
        for placement in read_placements:
            current = best_by_contig.get(placement[1])
            if current is None or placement[0] > current[0]:
                best_by_contig[placement[1]] = placement
        ranked = sorted(best_by_contig.values(), reverse=True, key=lambda item: item[0])
        if len(ranked) > 1:
            difference = ranked[0][0][0] - ranked[1][0][0]
            ambiguous = (
                difference <= minimum_confidence_margin
                if damage_aware_ranking
                else difference < minimum_overlap_margin
            )
            if ambiguous:
                read_index = ranked[0][2].read_index
                if read_clusters[read_index] is None:
                    diagnostics.deferred_ambiguous_assignments += 1
                else:
                    diagnostics.cross_cluster_ambiguous_assignments += 1
                continue
        _, contig_index, candidate = ranked[0]
        selected_by_contig[contig_index].append(candidate)
        if read_clusters[candidate.read_index] is None:
            diagnostics.deferred_unique_assignments += 1
        else:
            diagnostics.cross_cluster_unique_assignments += 1

    reciprocal_diagnostics = _ReciprocalDiagnostics()
    reciprocal_by_contig: dict[int, list[_Alignment]] = {}
    for contig_index, selected in selected_by_contig.items():
        reciprocal_by_contig[contig_index] = _filter_reciprocal_best_extensions(
            contigs[contig_index].sequence,
            selected,
            cluster_indices[contig_index],
            reads,
            molecule_ids,
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
            minimum_overlap_margin=minimum_overlap_margin,
            diagnostics=reciprocal_diagnostics,
        )
    for recruited in reciprocal_by_contig.values():
        for candidate in recruited:
            if read_clusters[candidate.read_index] is None:
                diagnostics.deferred_reciprocal_assignments += 1
            else:
                diagnostics.cross_cluster_reciprocal_assignments += 1

    def record_supported_extension(
        contig: Read,
        recruited: list[_Alignment],
    ) -> int:
        if not recruited:
            return 0
        result = _consensus(
            contig.sequence,
            recruited,
            minimum_extension_support=minimum_consensus_support,
            dominance_ratio=dominance_ratio,
            damage_end_window=0,
            allow_internal_consensus=False,
        )
        return len(result.sequence) - len(contig.sequence)

    for contig_index, recruited in reciprocal_by_contig.items():
        contig = contigs[contig_index]
        added = record_supported_extension(contig, recruited)
        deferred_added = record_supported_extension(
            contig,
            [
                candidate
                for candidate in recruited
                if read_clusters[candidate.read_index] is None
            ],
        )
        cross_cluster_added = record_supported_extension(
            contig,
            [
                candidate
                for candidate in recruited
                if read_clusters[candidate.read_index] is not None
            ],
        )
        if added:
            diagnostics.supported_contigs += 1
            diagnostics.supported_bases += added
        if deferred_added:
            diagnostics.deferred_supported_contigs += 1
            diagnostics.deferred_supported_bases += deferred_added
        if cross_cluster_added:
            diagnostics.cross_cluster_supported_contigs += 1
            diagnostics.cross_cluster_supported_bases += cross_cluster_added
    return diagnostics


def _spanning_molecule_support(
    left: str,
    right: str,
    overlap: int,
    read_indices: set[int],
    reads: list[Read],
    molecule_ids: list[int],
    *,
    minimum_identity: float,
    minimum_ry_identity: float,
) -> int:
    """Count molecules spanning unique sequence on both sides of an overlap."""
    merged = left + right[overlap:]
    overlap_start = len(left) - overlap
    overlap_end = len(left)
    supporting: set[int] = set()
    for read_index in read_indices:
        sequence = reads[read_index].sequence
        for oriented in (sequence, reverse_complement(sequence)):
            first_offset = max(0, overlap_end - len(oriented) + 1)
            last_offset = min(overlap_start - 1, len(merged) - len(oriented))
            if first_offset > last_offset:
                continue
            for offset in range(first_offset, last_offset + 1):
                target = merged[offset : offset + len(oriented)]
                if (
                    _identity(target, oriented) >= minimum_identity
                    and _identity(_ry(target), _ry(oriented)) >= minimum_ry_identity
                ):
                    supporting.add(molecule_ids[read_index])
                    break
            if molecule_ids[read_index] in supporting:
                break
    return len(supporting)


def _audit_read_supported_contig_links(
    contigs: list[Read],
    cluster_indices: list[set[int]],
    reads: list[Read],
    molecule_ids: list[int],
    *,
    anchor_k: int,
    anchors_per_read: int,
    maximum_anchor_occurrences: int,
    minimum_anchor_matches: int,
    minimum_overlap: int,
    minimum_identity: float,
    minimum_ry_identity: float,
    minimum_overlap_margin: int,
    minimum_read_support: int,
) -> _ContigLinkDiagnostics:
    """Project only reciprocal contig links crossed by independent molecules."""
    diagnostics = _ContigLinkDiagnostics()
    if len(contigs) < 2:
        return diagnostics
    position_bits = max(
        1, max(len(contig.sequence) for contig in contigs).bit_length()
    )
    target_window = max(len(contig.sequence) for contig in contigs)
    anchors = _anchor_index(
        contigs,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    contig_ids = list(range(len(contigs)))
    reciprocal_diagnostics = _ReciprocalDiagnostics()
    evidenced: dict[tuple[int, str, int, str], _SupportedContigLink] = {}
    for source_index, source in enumerate(contigs):
        candidates = _candidate_alignments(
            source.sequence,
            contigs,
            contig_ids,
            anchors,
            {source_index},
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
        proper = [
            candidate
            for candidate in candidates
            if candidate.offset < 0
            or candidate.offset + len(candidate.sequence) > len(source.sequence)
        ]
        diagnostics.candidate_overlaps += len(proper)
        selected, _, ambiguous = _unique_best_extensions(
            source.sequence,
            proper,
            minimum_overlap_margin=minimum_overlap_margin,
        )
        diagnostics.unique_overlaps += len(selected)
        diagnostics.ambiguous_ends += ambiguous
        reciprocal = _filter_reciprocal_best_extensions(
            source.sequence,
            selected,
            {source_index},
            contigs,
            contig_ids,
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
            minimum_overlap_margin=minimum_overlap_margin,
            diagnostics=reciprocal_diagnostics,
        )
        diagnostics.reciprocal_overlaps += len(reciprocal)
        for choice in reciprocal:
            target_index = choice.read_index
            left_extension = max(0, -choice.offset)
            right_extension = max(
                0,
                choice.offset + len(choice.sequence) - len(source.sequence),
            )
            if bool(left_extension) == bool(right_extension):
                continue
            reversed_target = choice.sequence != contigs[target_index].sequence
            if right_extension:
                left_sequence = source.sequence
                right_sequence = choice.sequence
                overlap = len(source.sequence) - choice.offset
                source_end = "right"
                target_end = "right" if reversed_target else "left"
            else:
                left_sequence = choice.sequence
                right_sequence = source.sequence
                overlap = choice.offset + len(choice.sequence)
                source_end = "left"
                target_end = "left" if reversed_target else "right"
            evidence_reads = cluster_indices[source_index] | cluster_indices[target_index]
            support = _spanning_molecule_support(
                left_sequence,
                right_sequence,
                overlap,
                evidence_reads,
                reads,
                molecule_ids,
                minimum_identity=minimum_identity,
                minimum_ry_identity=minimum_ry_identity,
            )
            if source_index < target_index:
                key = (source_index, source_end, target_index, target_end)
                link = _SupportedContigLink(
                    source_index,
                    target_index,
                    source_end,
                    target_end,
                    overlap,
                    support,
                )
            else:
                key = (target_index, target_end, source_index, source_end)
                link = _SupportedContigLink(
                    target_index,
                    source_index,
                    target_end,
                    source_end,
                    overlap,
                    support,
                )
            current = evidenced.get(key)
            if current is None or (link.support, link.overlap) > (
                current.support,
                current.overlap,
            ):
                evidenced[key] = link

    support_values = [link.support for link in evidenced.values()]
    diagnostics.support_at_least_1 = sum(value >= 1 for value in support_values)
    diagnostics.support_at_least_2 = sum(value >= 2 for value in support_values)
    diagnostics.support_at_least_3 = sum(value >= 3 for value in support_values)
    diagnostics.support_at_least_5 = sum(value >= 5 for value in support_values)
    diagnostics.max_read_support = max(support_values, default=0)
    supported = {
        key: link
        for key, link in evidenced.items()
        if link.support >= minimum_read_support
    }
    diagnostics.read_supported_overlaps = len(supported)
    end_counts: dict[tuple[int, str], int] = defaultdict(int)
    for link in supported.values():
        end_counts[(link.first, link.first_end)] += 1
        end_counts[(link.second, link.second_end)] += 1
    diagnostics.ambiguous_ends += sum(count > 1 for count in end_counts.values())
    linear_links = [
        link
        for link in supported.values()
        if end_counts[(link.first, link.first_end)] == 1
        and end_counts[(link.second, link.second_end)] == 1
    ]
    adjacency: dict[int, list[tuple[int, _SupportedContigLink]]] = defaultdict(list)
    for link in linear_links:
        adjacency[link.first].append((link.second, link))
        adjacency[link.second].append((link.first, link))
    visited: set[int] = set()
    for start in adjacency:
        if start in visited:
            continue
        stack = [start]
        component: set[int] = set()
        component_links: dict[tuple[int, str, int, str], _SupportedContigLink] = {}
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            visited.add(node)
            for neighbor, link in adjacency[node]:
                component_links[
                    (link.first, link.first_end, link.second, link.second_end)
                ] = link
                stack.append(neighbor)
        if len(component_links) == len(component):
            diagnostics.cyclic_components += 1
            continue
        diagnostics.linear_chains += 1
        diagnostics.projected_joins += len(component_links)
        merged_bases = sum(len(contigs[index].sequence) for index in component)
        merged_bases -= sum(link.overlap for link in component_links.values())
        diagnostics.projected_merged_bases += merged_bases
        diagnostics.projected_longest_contig = max(
            diagnostics.projected_longest_contig,
            merged_bases,
        )
    return diagnostics


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
    reciprocal_best_extension: bool = False,
    reciprocal_diagnostics: _ReciprocalDiagnostics | None = None,
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
        cluster_indices = {center_index}
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
            if reciprocal_best_extension and selected:
                if reciprocal_diagnostics is None:
                    reciprocal_diagnostics = _ReciprocalDiagnostics()
                selected = _filter_reciprocal_best_extensions(
                    target,
                    selected,
                    cluster_indices,
                    contigs,
                    molecule_ids,
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
                    minimum_overlap_margin=minimum_overlap_margin,
                    diagnostics=reciprocal_diagnostics,
                )
            for read_index in contained:
                assigned[read_index] = True
            if not selected:
                break
            for candidate in selected:
                assigned[candidate.read_index] = True
                cluster_indices.add(candidate.read_index)
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


def _n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    threshold = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            return length
    return 0


def audit_two_tier_redundancy(
    primary_contigs: list[Read],
    rescue_contigs: list[Read],
    *,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_identity: float = 0.97,
    minimum_ry_identity: float = 0.99,
    minimum_coverage: float = 0.99,
) -> TwoTierRedundancyDiagnostics:
    """Project novel rescue contigs without changing the primary assembly."""
    if not rescue_contigs:
        lengths = [len(contig.sequence) for contig in primary_contigs]
        return TwoTierRedundancyDiagnostics(
            replacement_projected_bases=sum(lengths),
            replacement_projected_n50=_n50(lengths),
            replacement_projected_longest_contig=max(lengths, default=0),
            projected_contigs=len(primary_contigs),
            projected_bases=sum(lengths),
            projected_n50=_n50(lengths),
            projected_longest_contig=max(lengths, default=0),
        )
    if not 0.0 <= minimum_identity <= 1.0:
        raise ValueError("minimum_identity must be between zero and one")
    if not 0.0 <= minimum_ry_identity <= 1.0:
        raise ValueError("minimum_ry_identity must be between zero and one")
    if not 0.0 < minimum_coverage <= 1.0:
        raise ValueError("minimum_coverage must be greater than zero and at most one")

    ordered_rescue = sorted(
        rescue_contigs,
        key=lambda contig: (-len(contig.sequence), contig.name, contig.sequence),
    )
    representatives = list(primary_contigs) + ordered_rescue
    primary_count = len(primary_contigs)
    position_bits = max(
        1, max(len(contig.sequence) for contig in representatives).bit_length()
    )
    anchors = _anchor_index(
        representatives,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    unavailable = [True] * len(representatives)
    for index in range(primary_count):
        unavailable[index] = False
    molecule_ids = list(range(len(representatives)))
    contained_by_primary = redundant_with_rescue = extends_primary = 0
    ambiguous_primary_extensions = 0
    replacements: dict[int, Read] = {}
    novel: list[Read] = []

    for rescue_offset, rescue in enumerate(ordered_rescue):
        rescue_index = primary_count + rescue_offset
        candidates = _candidate_alignments(
            rescue.sequence,
            representatives,
            molecule_ids,
            anchors,
            unavailable,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=anchor_k,
            minimum_identity=minimum_identity,
            minimum_ry_identity=minimum_ry_identity,
            position_bits=position_bits,
            target_window=max(len(rescue.sequence), anchor_k),
        )
        query_contained_by_primary = False
        query_contained_by_rescue = False
        extended_primary_indices: set[int] = set()
        for candidate in candidates:
            start = max(0, candidate.offset)
            stop = min(
                len(rescue.sequence),
                candidate.offset + len(candidate.sequence),
            )
            overlap = max(0, stop - start)
            query_coverage = overlap / len(rescue.sequence)
            representative_coverage = overlap / len(candidate.sequence)
            if query_coverage >= minimum_coverage:
                if candidate.read_index < primary_count:
                    query_contained_by_primary = True
                else:
                    query_contained_by_rescue = True
            elif (
                candidate.read_index < primary_count
                and representative_coverage >= minimum_coverage
            ):
                extended_primary_indices.add(candidate.read_index)

        if query_contained_by_primary:
            contained_by_primary += 1
            continue
        if query_contained_by_rescue:
            redundant_with_rescue += 1
            continue
        if extended_primary_indices:
            extends_primary += 1
            if len(extended_primary_indices) != 1:
                ambiguous_primary_extensions += 1
                continue
            primary_index = next(iter(extended_primary_indices))
            previous = replacements.get(primary_index)
            if previous is None or len(rescue.sequence) > len(previous.sequence):
                replacements[primary_index] = rescue
            continue
        novel.append(rescue)
        unavailable[rescue_index] = False

    projected_lengths = [
        len(contig.sequence) for contig in primary_contigs + novel
    ]
    replacement_lengths = [
        len(replacements.get(index, contig).sequence)
        for index, contig in enumerate(primary_contigs)
    ]
    return TwoTierRedundancyDiagnostics(
        rescue_contigs=len(rescue_contigs),
        rescue_bases=sum(len(contig.sequence) for contig in rescue_contigs),
        contained_by_primary=contained_by_primary,
        redundant_with_rescue=redundant_with_rescue,
        extends_primary=extends_primary,
        ambiguous_primary_extensions=ambiguous_primary_extensions,
        replacement_contigs=len(replacements),
        replacement_added_bases=sum(
            len(replacement.sequence) - len(primary_contigs[index].sequence)
            for index, replacement in replacements.items()
        ),
        replacement_projected_bases=sum(replacement_lengths),
        replacement_projected_n50=_n50(replacement_lengths),
        replacement_projected_longest_contig=max(
            replacement_lengths, default=0
        ),
        novel_contigs=len(novel),
        novel_bases=sum(len(contig.sequence) for contig in novel),
        projected_contigs=len(projected_lengths),
        projected_bases=sum(projected_lengths),
        projected_n50=_n50(projected_lengths),
        projected_longest_contig=max(projected_lengths, default=0),
    )


def _filter_redundant_projection(
    contigs: list[Read],
    *,
    anchor_k: int,
    anchors_per_read: int,
    maximum_anchor_occurrences: int,
    minimum_anchor_matches: int,
    minimum_identity: float = 0.97,
    minimum_ry_identity: float = 0.99,
    minimum_coverage: float = 0.99,
) -> tuple[list[Read], int]:
    """Keep longest representatives under final-output containment rules."""
    if not contigs:
        return [], 0
    ordered = sorted(
        contigs,
        key=lambda contig: (-len(contig.sequence), contig.name, contig.sequence),
    )
    exact_representatives: list[Read] = []
    exact_keys: set[str] = set()
    for contig in ordered:
        reverse = reverse_complement(contig.sequence)
        key = min(contig.sequence, reverse)
        if key in exact_keys:
            continue
        exact_keys.add(key)
        exact_representatives.append(contig)
    exact_redundant = len(ordered) - len(exact_representatives)
    ordered = exact_representatives
    position_bits = max(
        1, max(len(contig.sequence) for contig in ordered).bit_length()
    )
    anchors = _anchor_index(
        ordered,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    unavailable = [True] * len(ordered)
    molecule_ids = list(range(len(ordered)))
    retained: list[Read] = []
    for index, query in enumerate(ordered):
        candidates = _candidate_alignments(
            query.sequence,
            ordered,
            molecule_ids,
            anchors,
            unavailable,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=anchor_k,
            minimum_identity=minimum_identity,
            minimum_ry_identity=minimum_ry_identity,
            position_bits=position_bits,
            target_window=max(len(query.sequence), anchor_k),
        )
        contained = any(
            (
                min(
                    len(query.sequence),
                    candidate.offset + len(candidate.sequence),
                )
                - max(0, candidate.offset)
            )
            / len(query.sequence)
            >= minimum_coverage
            for candidate in candidates
        )
        if contained:
            continue
        retained.append(query)
        unavailable[index] = False
    return retained, exact_redundant + len(ordered) - len(retained)


def audit_master_overlap_graph(
    contigs: list[Read],
    *,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
) -> tuple[list[Read], MasterOverlapGraphDiagnostics]:
    """Compatibility wrapper for the graph projection implementation."""
    from anvaya.overlap_graph import audit_master_overlap_graph

    return audit_master_overlap_graph(
        contigs,
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_overlap=minimum_overlap,
    )


def audit_raw_confirmed_master_overlap_graph(
    contigs: list[Read],
    raw_reads: list[Read],
    *,
    molecule_ids: list[int] | None = None,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
    minimum_identity: float = 0.98,
    maximum_mismatches: int = 1,
    minimum_primary_allele_support: int = 3,
    minimum_alternate_allele_support: int = 2,
    minimum_support_margin: int = 2,
    minimum_base_quality: int = 20,
    damage_end_window: int = 5,
) -> tuple[list[Read], RawConfirmedMasterGraphDiagnostics]:
    """Compatibility wrapper for the graph projection implementation."""
    from anvaya.overlap_graph import audit_raw_confirmed_master_overlap_graph

    return audit_raw_confirmed_master_overlap_graph(
        contigs,
        raw_reads,
        molecule_ids=molecule_ids,
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_overlap=minimum_overlap,
        minimum_identity=minimum_identity,
        maximum_mismatches=maximum_mismatches,
        minimum_primary_allele_support=minimum_primary_allele_support,
        minimum_alternate_allele_support=minimum_alternate_allele_support,
        minimum_support_margin=minimum_support_margin,
        minimum_base_quality=minimum_base_quality,
        damage_end_window=damage_end_window,
    )


def audit_strain_safe_containment(
    contigs: list[Read],
    raw_reads: list[Read],
    *,
    molecule_ids: list[int] | None = None,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_identity: float = 0.99,
    minimum_coverage: float = 0.99,
    minimum_primary_allele_support: int = 3,
    minimum_alternate_allele_support: int = 2,
    minimum_support_margin: int = 2,
    minimum_base_quality: int = 20,
    damage_end_window: int = 5,
) -> tuple[list[Read], StrainSafeContainmentDiagnostics]:
    """Compatibility wrapper for the graph projection implementation."""
    from anvaya.overlap_graph import audit_strain_safe_containment

    return audit_strain_safe_containment(
        contigs,
        raw_reads,
        molecule_ids=molecule_ids,
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_identity=minimum_identity,
        minimum_coverage=minimum_coverage,
        minimum_primary_allele_support=minimum_primary_allele_support,
        minimum_alternate_allele_support=minimum_alternate_allele_support,
        minimum_support_margin=minimum_support_margin,
        minimum_base_quality=minimum_base_quality,
        damage_end_window=damage_end_window,
    )


def write_iterative_reclustering_report(
    diagnostics: IterativeReclusteringDiagnostics,
    path: str | Path,
) -> None:
    """Compatibility wrapper for the reclustering report writer."""
    from anvaya.overlap_reclustering import write_iterative_reclustering_report

    write_iterative_reclustering_report(diagnostics, path)


def audit_iterative_reclustering(
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
    dominance_ratio: float = 4.0,
    maximum_iterations: int = 4,
    minimum_overlap_margin: int = 3,
    minimum_ranked_extension_support: int = 1,
    reciprocal_best_extension: bool = True,
    damage_aware_ranking: bool = True,
    minimum_confidence_margin: float = 0.0,
    damage_mismatch_penalty: float = 0.25,
    ranking_damage_end_window: int = 5,
    extension_consensus: bool = True,
    minimum_output_length: int = 0,
    derived_minimum_identity: float = 0.99,
    derived_minimum_raw_support: int = 2,
    mismatch_minimum_raw_support: int = 3,
    mismatch_minimum_raw_margin: int = 2,
    mismatch_minimum_base_quality: int = 20,
) -> tuple[list[Read], IterativeReclusteringDiagnostics]:
    """Compatibility wrapper for the iterative reclustering implementation."""
    from anvaya.overlap_reclustering import audit_iterative_reclustering

    return audit_iterative_reclustering(
        reads,
        molecule_ids=molecule_ids,
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_overlap=minimum_overlap,
        minimum_identity=minimum_identity,
        minimum_ry_identity=minimum_ry_identity,
        minimum_cluster_size=minimum_cluster_size,
        minimum_consensus_support=minimum_consensus_support,
        dominance_ratio=dominance_ratio,
        maximum_iterations=maximum_iterations,
        minimum_overlap_margin=minimum_overlap_margin,
        minimum_ranked_extension_support=minimum_ranked_extension_support,
        reciprocal_best_extension=reciprocal_best_extension,
        damage_aware_ranking=damage_aware_ranking,
        minimum_confidence_margin=minimum_confidence_margin,
        damage_mismatch_penalty=damage_mismatch_penalty,
        ranking_damage_end_window=ranking_damage_end_window,
        extension_consensus=extension_consensus,
        minimum_output_length=minimum_output_length,
        derived_minimum_identity=derived_minimum_identity,
        derived_minimum_raw_support=derived_minimum_raw_support,
        mismatch_minimum_raw_support=mismatch_minimum_raw_support,
        mismatch_minimum_raw_margin=mismatch_minimum_raw_margin,
        mismatch_minimum_base_quality=mismatch_minimum_base_quality,
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
    ranked_extension: bool = False,
    extension_consensus: bool = False,
    minimum_ranked_extension_support: int = 1,
    reciprocal_best_extension: bool = False,
    damage_aware_ranking: bool = False,
    minimum_confidence_margin: float = 0.0,
    damage_mismatch_penalty: float = 0.25,
    ranking_damage_end_window: int = 5,
    cross_cluster_recruitment_audit: bool = False,
    read_supported_contig_link_audit: bool = False,
    minimum_contig_link_read_support: int = 2,
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
    if extension_consensus and not ranked_extension:
        raise ValueError("extension_consensus requires ranked_extension")
    if damage_aware_ranking and not ranked_extension:
        raise ValueError("damage_aware_ranking requires ranked_extension")
    if minimum_ranked_extension_support < 1:
        raise ValueError("minimum_ranked_extension_support must be at least 1")
    if minimum_confidence_margin < 0.0:
        raise ValueError("minimum_confidence_margin must not be negative")
    if not 0.0 <= damage_mismatch_penalty <= 1.0:
        raise ValueError("damage_mismatch_penalty must be between zero and one")
    if ranking_damage_end_window < 0:
        raise ValueError("ranking_damage_end_window must not be negative")
    if cross_cluster_recruitment_audit and maximum_contig_iterations:
        raise ValueError(
            "cross_cluster_recruitment_audit requires maximum_contig_iterations=0"
        )
    if read_supported_contig_link_audit and maximum_contig_iterations:
        raise ValueError(
            "read_supported_contig_link_audit requires maximum_contig_iterations=0"
        )
    if minimum_contig_link_read_support < 1:
        raise ValueError("minimum_contig_link_read_support must be at least 1")
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
    read_clusters: list[int | None] = [None] * len(reads)
    contigs: list[Read] = []
    contig_cluster_indices: list[set[int]] = []
    clusters = clustered_reads = extension_rounds = 0
    extended_contigs = added_bases = corrected_bases = 0
    ambiguous_extensions = 0
    correction_candidates = correction_terminal_only = 0
    correction_insufficient_support = correction_ambiguous = 0
    candidate_diagnostics = _CandidateDiagnostics()
    reciprocal_diagnostics = _ReciprocalDiagnostics()
    ranking_diagnostics = _RankingDiagnostics()
    targets_without_candidates = clusters_below_minimum_size = 0
    consensus_without_extension = 0
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
                diagnostics=candidate_diagnostics,
            )
            if not candidates:
                targets_without_candidates += 1
                break
            if not members and len(candidates) + 1 < minimum_cluster_size:
                clusters_below_minimum_size += 1
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
            if ranked_extension:
                internal = _consensus(
                    target,
                    tentative_members,
                    minimum_extension_support=len(tentative_members) + 2,
                    dominance_ratio=dominance_ratio,
                    damage_end_window=0,
                    minimum_internal_support=minimum_consensus_support,
                )
                selected, _, ambiguous = _unique_best_extensions(
                    internal.sequence,
                    candidates,
                    minimum_overlap_margin=minimum_overlap_margin,
                    damage_aware_ranking=damage_aware_ranking,
                    minimum_confidence_margin=minimum_confidence_margin,
                    damage_mismatch_penalty=damage_mismatch_penalty,
                    ranking_damage_end_window=ranking_damage_end_window,
                    diagnostics=ranking_diagnostics,
                )
                selected = _filter_extension_boundary_support(
                    internal.sequence,
                    candidates,
                    selected,
                    minimum_ranked_extension_support,
                )
                if reciprocal_best_extension and selected:
                    selected = _filter_reciprocal_best_extensions(
                        internal.sequence,
                        selected,
                        cluster_indices,
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
                        minimum_overlap_margin=minimum_overlap_margin,
                        diagnostics=reciprocal_diagnostics,
                    )
                ambiguous_extensions += ambiguous
                template = _consensus(
                    internal.sequence,
                    selected,
                    minimum_extension_support=1,
                    dominance_ratio=dominance_ratio,
                    damage_end_window=0,
                    allow_internal_consensus=False,
                )
                result = template
                if extension_consensus and selected:
                    shifted_members = tentative_members
                    if template.left_extension:
                        shifted_members = [
                            replace(
                                member,
                                offset=member.offset + template.left_extension,
                            )
                            for member in tentative_members
                        ]
                    recalled = _consensus(
                        template.sequence,
                        shifted_members,
                        minimum_extension_support=len(shifted_members) + 2,
                        minimum_internal_support=minimum_consensus_support,
                        dominance_ratio=dominance_ratio,
                        damage_end_window=0,
                    )
                    result = replace(
                        recalled,
                        left_extension=template.left_extension,
                        right_extension=template.right_extension,
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
                consensus_without_extension += 1
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
        contig_cluster_indices.append(set(cluster_indices))
        for read_index in cluster_indices:
            read_clusters[read_index] = len(contigs) - 1

    recruitment_diagnostics = _RecruitmentDiagnostics()
    if cross_cluster_recruitment_audit:
        recruitment_diagnostics = _audit_cross_cluster_recruitment(
            contigs,
            contig_cluster_indices,
            reads,
            read_clusters,
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
            minimum_overlap_margin=minimum_overlap_margin,
            minimum_consensus_support=minimum_consensus_support,
            dominance_ratio=dominance_ratio,
            damage_aware_ranking=damage_aware_ranking,
            minimum_confidence_margin=minimum_confidence_margin,
            damage_mismatch_penalty=damage_mismatch_penalty,
            ranking_damage_end_window=ranking_damage_end_window,
        )

    contig_link_diagnostics = _ContigLinkDiagnostics()
    if read_supported_contig_link_audit:
        contig_link_diagnostics = _audit_read_supported_contig_links(
            contigs,
            contig_cluster_indices,
            reads,
            molecules,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=minimum_overlap,
            minimum_identity=minimum_identity,
            minimum_ry_identity=minimum_ry_identity,
            minimum_overlap_margin=minimum_overlap_margin,
            minimum_read_support=minimum_contig_link_read_support,
        )

    contig_iterations = contig_merges = 0
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
            reciprocal_best_extension=reciprocal_best_extension,
            reciprocal_diagnostics=reciprocal_diagnostics,
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
        candidate_offsets=candidate_diagnostics.offsets,
        candidate_below_anchor_support=candidate_diagnostics.below_anchor_support,
        candidate_short_overlap=candidate_diagnostics.short_overlap,
        candidate_dna_rejected=candidate_diagnostics.dna_rejected,
        candidate_ry_rejected=candidate_diagnostics.ry_rejected,
        candidate_molecule_ambiguous=candidate_diagnostics.molecule_ambiguous,
        candidate_alignments=candidate_diagnostics.alignments,
        unavailable_anchor_hits=candidate_diagnostics.unavailable_anchor_hits,
        targets_without_candidates=targets_without_candidates,
        clusters_below_minimum_size=clusters_below_minimum_size,
        consensus_without_extension=consensus_without_extension,
        reciprocal_extension_checks=reciprocal_diagnostics.checks,
        reciprocal_extension_rejections=reciprocal_diagnostics.rejections,
        confidence_ranked_sides=ranking_diagnostics.ranked_sides,
        confidence_changed_winners=ranking_diagnostics.changed_winners,
        confidence_ambiguous_extensions=ranking_diagnostics.ambiguous_extensions,
        deferred_reads=recruitment_diagnostics.deferred_reads,
        cross_cluster_reads=recruitment_diagnostics.cross_cluster_reads,
        deferred_candidate_reads=recruitment_diagnostics.deferred_candidate_reads,
        cross_cluster_candidate_reads=(
            recruitment_diagnostics.cross_cluster_candidate_reads
        ),
        deferred_extension_candidates=(
            recruitment_diagnostics.deferred_extension_candidates
        ),
        cross_cluster_extension_candidates=(
            recruitment_diagnostics.cross_cluster_extension_candidates
        ),
        deferred_ambiguous_assignments=(
            recruitment_diagnostics.deferred_ambiguous_assignments
        ),
        cross_cluster_ambiguous_assignments=(
            recruitment_diagnostics.cross_cluster_ambiguous_assignments
        ),
        deferred_unique_assignments=(
            recruitment_diagnostics.deferred_unique_assignments
        ),
        cross_cluster_unique_assignments=(
            recruitment_diagnostics.cross_cluster_unique_assignments
        ),
        deferred_reciprocal_assignments=(
            recruitment_diagnostics.deferred_reciprocal_assignments
        ),
        cross_cluster_reciprocal_assignments=(
            recruitment_diagnostics.cross_cluster_reciprocal_assignments
        ),
        recruitment_supported_contigs=recruitment_diagnostics.supported_contigs,
        recruitment_supported_bases=recruitment_diagnostics.supported_bases,
        deferred_supported_contigs=(
            recruitment_diagnostics.deferred_supported_contigs
        ),
        deferred_supported_bases=recruitment_diagnostics.deferred_supported_bases,
        cross_cluster_supported_contigs=(
            recruitment_diagnostics.cross_cluster_supported_contigs
        ),
        cross_cluster_supported_bases=(
            recruitment_diagnostics.cross_cluster_supported_bases
        ),
        contig_link_candidate_overlaps=contig_link_diagnostics.candidate_overlaps,
        contig_link_unique_overlaps=contig_link_diagnostics.unique_overlaps,
        contig_link_reciprocal_overlaps=contig_link_diagnostics.reciprocal_overlaps,
        contig_link_read_supported_overlaps=(
            contig_link_diagnostics.read_supported_overlaps
        ),
        contig_link_support_at_least_1=contig_link_diagnostics.support_at_least_1,
        contig_link_support_at_least_2=contig_link_diagnostics.support_at_least_2,
        contig_link_support_at_least_3=contig_link_diagnostics.support_at_least_3,
        contig_link_support_at_least_5=contig_link_diagnostics.support_at_least_5,
        contig_link_max_read_support=contig_link_diagnostics.max_read_support,
        contig_link_ambiguous_ends=contig_link_diagnostics.ambiguous_ends,
        contig_link_cyclic_components=contig_link_diagnostics.cyclic_components,
        contig_link_linear_chains=contig_link_diagnostics.linear_chains,
        contig_link_projected_joins=contig_link_diagnostics.projected_joins,
        contig_link_projected_merged_bases=(
            contig_link_diagnostics.projected_merged_bases
        ),
        contig_link_projected_longest_contig=(
            contig_link_diagnostics.projected_longest_contig
        ),
    )
