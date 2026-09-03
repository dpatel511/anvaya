"""Fixed-pool iterative overlap reclustering and reporting."""

import csv
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

from anvaya.damage_consensus import _anchor_index, _identity
from anvaya.overlap_assembly import (
    IterativeReclusteringDiagnostics,
    IterativeReclusteringRound,
    _Alignment,
    _ReciprocalDiagnostics,
    _candidate_alignments,
    _consensus,
    _filter_extension_boundary_support,
    _filter_raw_supported_extensions,
    _filter_reciprocal_best_extensions,
    _filter_redundant_projection,
    _n50,
    _overlap_confidence,
    _projection_regresses,
    _raw_supports_target_mismatches,
    _unique_best_extensions,
)
from anvaya.reads import Read


def write_iterative_reclustering_report(
    diagnostics: IterativeReclusteringDiagnostics,
    path: str | Path,
) -> None:
    """Write one deterministic row per global reclustering round."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        field.name
        for field in IterativeReclusteringRound.__dataclass_fields__.values()
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for round_diagnostics in diagnostics.rounds:
            writer.writerow(
                {field: getattr(round_diagnostics, field) for field in fields}
            )


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
    """Project fixed-pool iterative assembly without changing accepted output."""
    if not reads:
        return [], IterativeReclusteringDiagnostics(converged=True)
    if maximum_iterations < 1:
        raise ValueError("maximum_iterations must be at least 1")
    if minimum_cluster_size < 2:
        raise ValueError("minimum_cluster_size must be at least 2")
    if minimum_consensus_support < 2:
        raise ValueError("minimum_consensus_support must be at least 2")
    if minimum_ranked_extension_support < 1:
        raise ValueError("minimum_ranked_extension_support must be at least 1")
    if minimum_output_length < 0:
        raise ValueError("minimum_output_length must not be negative")
    if not 0.0 <= derived_minimum_identity <= 1.0:
        raise ValueError("derived_minimum_identity must be between 0 and 1")
    if derived_minimum_raw_support < 1:
        raise ValueError("derived_minimum_raw_support must be at least 1")
    if mismatch_minimum_raw_support < 1:
        raise ValueError("mismatch_minimum_raw_support must be at least 1")
    if mismatch_minimum_raw_margin < 1:
        raise ValueError("mismatch_minimum_raw_margin must be at least 1")
    if mismatch_minimum_base_quality < 0:
        raise ValueError("mismatch_minimum_base_quality must not be negative")
    molecules = list(range(len(reads))) if molecule_ids is None else molecule_ids
    if len(molecules) != len(reads):
        raise ValueError("molecule IDs must align one-to-one with reads")

    raw_reads = [Read(read.name, read.sequence, read.qualities) for read in reads]
    current = list(raw_reads)
    source_lengths = [len(read.sequence) for read in reads]
    assembled = [False] * len(reads)
    derived = [False] * len(reads)
    round_reports: list[IterativeReclusteringRound] = []
    total_extended = total_added = total_candidates = 0
    total_reused = total_conflicts = 0
    total_reciprocal_checks = total_reciprocal_rejections = 0
    total_derived_rejections = total_raw_support_rejections = 0
    total_raw_confirmed_mismatches = 0
    converged = False
    rollback_rejected = False
    projected: list[Read] = []
    raw_position_bits = max(
        1, max(len(read.sequence) for read in raw_reads).bit_length()
    )
    raw_target_window = max(len(read.sequence) for read in raw_reads)
    raw_anchors = _anchor_index(
        raw_reads,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    raw_provenance = [False] * len(raw_reads)

    for iteration in range(1, maximum_iterations + 1):
        current_before = list(current)
        assembled_before = list(assembled)
        derived_before = list(derived)
        position_bits = max(
            1, max(len(read.sequence) for read in current).bit_length()
        )
        target_window = max(len(read.sequence) for read in current)
        anchors = _anchor_index(
            current,
            anchor_k,
            0,
            anchors_per_read,
            maximum_anchor_occurrences,
        )
        proposals: dict[
            int,
            tuple[str, list[_Alignment], list[_Alignment]],
        ] = {}
        reciprocal = _ReciprocalDiagnostics()
        derived_overlap_rejections = raw_support_rejections = 0
        raw_confirmed_mismatch_overlaps = 0
        order = sorted(
            range(len(current)),
            key=lambda index: (-len(current[index].sequence), index),
        )
        for center_index in order:
            target = current[center_index].sequence
            candidates = _candidate_alignments(
                target,
                current,
                molecules,
                anchors,
                {center_index},
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
            if iteration > 1:
                raw_candidates: list[_Alignment] = []
                if any(
                    derived[center_index] or derived[candidate.read_index]
                    for candidate in candidates
                ):
                    raw_candidates = _candidate_alignments(
                        target,
                        raw_reads,
                        molecules,
                        raw_anchors,
                        set(),
                        anchor_k=anchor_k,
                        anchors_per_read=anchors_per_read,
                        maximum_anchor_occurrences=maximum_anchor_occurrences,
                        minimum_anchor_matches=minimum_anchor_matches,
                        minimum_overlap=minimum_overlap,
                        minimum_identity=minimum_identity,
                        minimum_ry_identity=minimum_ry_identity,
                        position_bits=raw_position_bits,
                        target_window=raw_target_window,
                    )
                retained: list[_Alignment] = []
                for candidate in candidates:
                    if not (derived[center_index] or derived[candidate.read_index]):
                        retained.append(candidate)
                        continue
                    start = max(0, candidate.offset)
                    stop = min(
                        len(target), candidate.offset + len(candidate.sequence)
                    )
                    candidate_start = start - candidate.offset
                    candidate_stop = stop - candidate.offset
                    if _identity(
                        target[start:stop],
                        candidate.sequence[candidate_start:candidate_stop],
                    ) >= derived_minimum_identity:
                        retained.append(candidate)
                    elif _raw_supports_target_mismatches(
                        target,
                        candidate,
                        raw_candidates,
                        raw_reads,
                        molecules,
                        mismatch_minimum_raw_support,
                        mismatch_minimum_raw_margin,
                        mismatch_minimum_base_quality,
                    ):
                        retained.append(candidate)
                        raw_confirmed_mismatch_overlaps += 1
                    else:
                        derived_overlap_rejections += 1
                candidates = retained
            if len(candidates) + 1 < minimum_cluster_size:
                continue
            internal = _consensus(
                target,
                candidates,
                minimum_extension_support=len(candidates) + 2,
                minimum_internal_support=minimum_consensus_support,
                dominance_ratio=dominance_ratio,
                damage_end_window=0,
            )
            selected, _, _ = _unique_best_extensions(
                internal.sequence,
                candidates,
                minimum_overlap_margin=minimum_overlap_margin,
                damage_aware_ranking=damage_aware_ranking and iteration == 1,
                minimum_confidence_margin=minimum_confidence_margin,
                damage_mismatch_penalty=damage_mismatch_penalty,
                ranking_damage_end_window=ranking_damage_end_window,
            )
            selected = _filter_extension_boundary_support(
                internal.sequence,
                candidates,
                selected,
                minimum_ranked_extension_support,
            )
            if iteration > 1 and selected:
                raw_candidates = _candidate_alignments(
                    internal.sequence,
                    raw_reads,
                    molecules,
                    raw_anchors,
                    set(),
                    anchor_k=anchor_k,
                    anchors_per_read=anchors_per_read,
                    maximum_anchor_occurrences=maximum_anchor_occurrences,
                    minimum_anchor_matches=minimum_anchor_matches,
                    minimum_overlap=minimum_overlap,
                    minimum_identity=minimum_identity,
                    minimum_ry_identity=minimum_ry_identity,
                    position_bits=raw_position_bits,
                    target_window=raw_target_window,
                )
                selected, rejected = _filter_raw_supported_extensions(
                    internal.sequence,
                    raw_candidates,
                    selected,
                    molecules,
                    raw_provenance,
                    derived_minimum_raw_support,
                )
                raw_support_rejections += rejected
            if reciprocal_best_extension and selected:
                selected = _filter_reciprocal_best_extensions(
                    internal.sequence,
                    selected,
                    {
                        center_index,
                        *(candidate.read_index for candidate in candidates),
                    },
                    current,
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
                    diagnostics=reciprocal,
                )
            if selected:
                proposals[center_index] = (internal.sequence, candidates, selected)

        placements: dict[
            int,
            list[tuple[tuple[float, int, int, int], int, _Alignment]],
        ] = defaultdict(list)
        for center_index, (target, _, selected) in proposals.items():
            for candidate in selected:
                start = max(0, candidate.offset)
                stop = min(len(target), candidate.offset + len(candidate.sequence))
                overlap = stop - start
                extension = max(
                    0,
                    -candidate.offset,
                    candidate.offset + len(candidate.sequence) - len(target),
                )
                confidence = (
                    _overlap_confidence(
                        target,
                        candidate,
                        damage_mismatch_penalty=damage_mismatch_penalty,
                        damage_end_window=ranking_damage_end_window,
                    )
                    if damage_aware_ranking and iteration == 1
                    else float(overlap)
                )
                placements[molecules[candidate.read_index]].append(
                    (
                        (confidence, overlap, candidate.support, extension),
                        center_index,
                        candidate,
                    )
                )

        accepted_by_center: dict[int, list[_Alignment]] = defaultdict(list)
        reused_candidates = conflicting_candidates = 0
        for alternatives in placements.values():
            alternatives.sort(key=lambda item: (item[0], -item[1]), reverse=True)
            if len(alternatives) == 1:
                accepted_by_center[alternatives[0][1]].append(alternatives[0][2])
                continue
            reused_candidates += 1
            first, second = alternatives[:2]
            difference = first[0][0] - second[0][0]
            ambiguous = (
                difference <= minimum_confidence_margin
                if damage_aware_ranking and iteration == 1
                else difference < minimum_overlap_margin
            )
            if ambiguous:
                conflicting_candidates += 1
                continue
            accepted_by_center[first[1]].append(first[2])

        next_pool = list(current)
        extended_sequences = added_bases = 0
        for center_index, selected in accepted_by_center.items():
            target, candidates, _ = proposals[center_index]
            template = _consensus(
                target,
                selected,
                minimum_extension_support=1,
                dominance_ratio=dominance_ratio,
                damage_end_window=0,
                allow_internal_consensus=False,
            )
            result = template
            if extension_consensus:
                shifted = candidates
                if template.left_extension:
                    shifted = [
                        replace(candidate, offset=candidate.offset + template.left_extension)
                        for candidate in candidates
                    ]
                recalled = _consensus(
                    template.sequence,
                    shifted,
                    minimum_extension_support=len(shifted) + 2,
                    minimum_internal_support=minimum_consensus_support,
                    dominance_ratio=dominance_ratio,
                    damage_end_window=0,
                )
                result = replace(
                    recalled,
                    left_extension=template.left_extension,
                    right_extension=template.right_extension,
                )
            added = len(result.sequence) - len(current[center_index].sequence)
            if added <= 0:
                continue
            next_pool[center_index] = Read(current[center_index].name, result.sequence)
            assembled[center_index] = True
            derived[center_index] = True
            extended_sequences += 1
            added_bases += added

        current = next_pool
        projected_candidates = [
            read
            for index, read in enumerate(current)
            if assembled[index]
            and len(read.sequence) > source_lengths[index]
            and len(read.sequence) >= minimum_output_length
        ]
        round_projection, redundant = _filter_redundant_projection(
            projected_candidates,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
        )
        lengths = [len(contig.sequence) for contig in round_projection]
        if round_reports and _projection_regresses(
            round_reports[-1].projected_n50,
            round_reports[-1].projected_bases,
            lengths,
        ):
            current = current_before
            assembled = assembled_before
            derived = derived_before
            converged = True
            rollback_rejected = True
            break
        projected = round_projection
        round_reports.append(IterativeReclusteringRound(
            iteration=iteration,
            extended_sequences=extended_sequences,
            added_bases=added_bases,
            candidate_extensions=sum(len(item[2]) for item in proposals.values()),
            reused_candidates=reused_candidates,
            conflicting_candidates=conflicting_candidates,
            reciprocal_checks=reciprocal.checks,
            reciprocal_rejections=reciprocal.rejections,
            derived_overlap_rejections=derived_overlap_rejections,
            raw_confirmed_mismatch_overlaps=raw_confirmed_mismatch_overlaps,
            raw_support_rejections=raw_support_rejections,
            assembled_sequences=len(projected_candidates),
            redundant_sequences=redundant,
            projected_contigs=len(projected),
            projected_bases=sum(lengths),
            projected_n50=_n50(lengths),
            projected_longest_contig=max(lengths, default=0),
        ))
        total_extended += extended_sequences
        total_added += added_bases
        total_candidates += round_reports[-1].candidate_extensions
        total_reused += reused_candidates
        total_conflicts += conflicting_candidates
        total_reciprocal_checks += reciprocal.checks
        total_reciprocal_rejections += reciprocal.rejections
        total_derived_rejections += derived_overlap_rejections
        total_raw_confirmed_mismatches += raw_confirmed_mismatch_overlaps
        total_raw_support_rejections += raw_support_rejections
        if added_bases == 0:
            converged = True
            break

    final = round_reports[-1]
    return projected, IterativeReclusteringDiagnostics(
        iterations=len(round_reports),
        converged=converged,
        rollback_rejected=rollback_rejected,
        extended_sequences=total_extended,
        added_bases=total_added,
        candidate_extensions=total_candidates,
        reused_candidates=total_reused,
        conflicting_candidates=total_conflicts,
        reciprocal_checks=total_reciprocal_checks,
        reciprocal_rejections=total_reciprocal_rejections,
        derived_overlap_rejections=total_derived_rejections,
        raw_confirmed_mismatch_overlaps=total_raw_confirmed_mismatches,
        raw_support_rejections=total_raw_support_rejections,
        assembled_sequences=final.assembled_sequences,
        redundant_sequences=final.redundant_sequences,
        projected_contigs=final.projected_contigs,
        projected_bases=final.projected_bases,
        projected_n50=final.projected_n50,
        projected_longest_contig=final.projected_longest_contig,
        rounds=tuple(round_reports),
    )
