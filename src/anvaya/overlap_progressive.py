"""Sequence lifecycle primitives for progressive overlap assembly."""

from dataclasses import dataclass, replace
from enum import Enum

from anvaya.damage_consensus import _anchor_index
from anvaya.overlap_assembly import (
    _Alignment,
    _RankingDiagnostics,
    _ReciprocalDiagnostics,
    _candidate_alignments,
    _consensus,
    _filter_extension_boundary_support,
    _filter_reciprocal_best_extensions,
    _unique_best_extensions,
)
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


_DAMAGE_PAIRS = {frozenset(("C", "T")), frozenset(("G", "A"))}


class SequenceState(str, Enum):
    """The role of a representative in progressive overlap assembly."""

    RAW = "raw"
    CORRECTED_CENTER = "corrected_center"
    EXTENDED_CONTIG = "extended_contig"
    CONSUMED = "consumed"


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    """A mutable-by-replacement representative backed by immutable raw evidence."""

    index: int
    molecule_id: int
    raw: Read
    current: Read
    state: SequenceState = SequenceState.RAW
    generation: int = 0
    contributing_molecules: frozenset[int] = frozenset()

    def corrected(self, read: Read) -> "SequenceRecord":
        """Return a corrected center without altering its raw observation."""
        if self.state is SequenceState.CONSUMED:
            raise ValueError("a consumed sequence cannot be corrected")
        return replace(
            self,
            current=read,
            state=SequenceState.CORRECTED_CENTER,
            generation=self.generation + 1,
        )

    def extended(
        self,
        read: Read,
        contributing_molecules: frozenset[int],
    ) -> "SequenceRecord":
        """Return an extended contig with its independent molecule provenance."""
        if self.state is SequenceState.CONSUMED:
            raise ValueError("a consumed sequence cannot be extended")
        return replace(
            self,
            current=read,
            state=SequenceState.EXTENDED_CONTIG,
            generation=self.generation + 1,
            contributing_molecules=(
                self.contributing_molecules | contributing_molecules
            ),
        )

    def consumed(self) -> "SequenceRecord":
        """Retire a representative while retaining its immutable evidence."""
        return replace(self, state=SequenceState.CONSUMED)


@dataclass(frozen=True, slots=True)
class ProgressiveSequencePool:
    """Raw observations and current representatives kept in separate views."""

    records: tuple[SequenceRecord, ...]

    @classmethod
    def from_reads(
        cls,
        reads: list[Read],
        molecule_ids: list[int] | None = None,
    ) -> "ProgressiveSequencePool":
        molecules = list(range(len(reads))) if molecule_ids is None else molecule_ids
        if len(molecules) != len(reads):
            raise ValueError("molecule IDs must align one-to-one with reads")
        return cls(
            tuple(
                SequenceRecord(
                    index=index,
                    molecule_id=molecules[index],
                    raw=read,
                    current=read,
                    contributing_molecules=frozenset({molecules[index]}),
                )
                for index, read in enumerate(reads)
            )
        )

    @property
    def raw_evidence(self) -> tuple[Read, ...]:
        """Return every original observation, including consumed representatives."""
        return tuple(record.raw for record in self.records)

    @property
    def active_raw(self) -> tuple[SequenceRecord, ...]:
        """Return unmodified raw representatives eligible for raw clustering."""
        return tuple(
            record for record in self.records if record.state is SequenceState.RAW
        )

    @property
    def active_derived(self) -> tuple[SequenceRecord, ...]:
        """Return corrected or extended representatives for contig merging."""
        return tuple(
            record
            for record in self.records
            if record.state
            in {SequenceState.CORRECTED_CENTER, SequenceState.EXTENDED_CONTIG}
        )

    def replace_record(self, record: SequenceRecord) -> "ProgressiveSequencePool":
        """Replace one record while preserving stable input ordering."""
        if not 0 <= record.index < len(self.records):
            raise IndexError("sequence record index is outside the pool")
        if self.records[record.index].index != record.index:
            raise ValueError("sequence record index does not match the pool")
        records = list(self.records)
        records[record.index] = record
        return ProgressiveSequencePool(tuple(records))


@dataclass(frozen=True, slots=True)
class ProgressiveRawCluster:
    """A deterministic raw-fragment cluster around its longest center."""

    center_index: int
    member_indices: tuple[int, ...]
    alignments: tuple[_Alignment, ...]


@dataclass(frozen=True, slots=True)
class ProgressiveClusteringDiagnostics:
    """Admission counts from one raw-only clustering pass."""

    candidate_centers: int = 0
    centers_with_candidates: int = 0
    clusters: int = 0
    clusters_below_minimum_size: int = 0
    clustered_reads: int = 0
    retained_unclustered_reads: int = 0


@dataclass(frozen=True, slots=True)
class ProgressiveRawExtensionDiagnostics:
    """Results from the raw-only correction and extension phase."""

    clusters: int = 0
    corrected_centers: int = 0
    extended_centers: int = 0
    extension_rounds: int = 0
    recruited_reads: int = 0
    added_bases: int = 0
    consumed_reads: int = 0
    ambiguous_extensions: int = 0
    reciprocal_checks: int = 0
    reciprocal_rejections: int = 0


@dataclass(frozen=True, slots=True)
class ProgressiveIterationDiagnostics:
    """Results from repeatedly extending derived centers with unused raw reads."""

    iterations: int = 0
    converged: bool = False
    candidate_alignments: int = 0
    extended_centers: int = 0
    added_bases: int = 0
    consumed_reads: int = 0
    ambiguous_extensions: int = 0
    insufficient_consensus_support: int = 0
    evidence_priority_sweeps: int = 0
    evidence_priority_reordered_centers: int = 0
    evidence_priority_claim_conflicts: int = 0
    rejected_left_sides: int = 0
    rejected_right_sides: int = 0
    rejected_support_1: int = 0
    rejected_support_2: int = 0
    rejected_support_3: int = 0
    rejected_support_4: int = 0
    rejected_support_5_plus: int = 0
    rejected_agreeing_sides: int = 0
    rejected_conflicting_sides: int = 0
    rejected_extension_1_5: int = 0
    rejected_extension_6_10: int = 0
    rejected_extension_11_20: int = 0
    rejected_extension_21_plus: int = 0
    rejected_high_quality_boundary_observations: int = 0
    rejected_low_quality_boundary_observations: int = 0
    rejected_missing_quality_boundary_observations: int = 0
    rejected_damage_compatible_mismatches: int = 0
    rejected_ordinary_mismatches: int = 0
    strict_boundary_accepted_sides: int = 0
    strict_boundary_conflict_rejections: int = 0
    strict_boundary_quality_rejections: int = 0
    projected_contigs: int = 0
    projected_bases: int = 0
    projected_n50: int = 0
    projected_longest_contig: int = 0


@dataclass(slots=True)
class _RejectedExtensionDiagnostics:
    left_sides: int = 0
    right_sides: int = 0
    support_1: int = 0
    support_2: int = 0
    support_3: int = 0
    support_4: int = 0
    support_5_plus: int = 0
    agreeing_sides: int = 0
    conflicting_sides: int = 0
    extension_1_5: int = 0
    extension_6_10: int = 0
    extension_11_20: int = 0
    extension_21_plus: int = 0
    high_quality_boundary_observations: int = 0
    low_quality_boundary_observations: int = 0
    missing_quality_boundary_observations: int = 0
    damage_compatible_mismatches: int = 0
    ordinary_mismatches: int = 0


def _oriented_quality(read: Read, sequence: str, position: int) -> int | None:
    """Return the quality corresponding to an oriented alignment base."""
    if read.qualities is None:
        return None
    if sequence == read.sequence:
        return read.qualities[position]
    if sequence == reverse_complement(read.sequence):
        return read.qualities[len(read.sequence) - position - 1]
    return None


def _n50(lengths: list[int]) -> int:
    total = sum(lengths)
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative * 2 >= total:
            return length
    return 0


def _filter_strict_boundary_evidence(
    target: str,
    candidates: list[_Alignment],
    selected: list[_Alignment],
    reads: list[Read],
    *,
    minimum_support: int,
    minimum_base_quality: int,
) -> tuple[list[_Alignment], int, int, int]:
    """Keep unanimous extensions backed by enough high-quality molecules."""
    accepted: list[_Alignment] = []
    accepted_sides = conflict_rejections = quality_rejections = 0
    for choice in selected:
        valid = True
        choice_sides = []
        if choice.offset < 0:
            choice_sides.append(-1)
        if choice.offset + len(choice.sequence) > len(target):
            choice_sides.append(len(target))
        for boundary in choice_sides:
            choice_base = choice.sequence[boundary - choice.offset]
            supporting = competing = 0
            for candidate in candidates:
                if not (
                    candidate.offset
                    <= boundary
                    < candidate.offset + len(candidate.sequence)
                ):
                    continue
                position = boundary - candidate.offset
                base = candidate.sequence[position]
                if base != choice_base:
                    competing += 1
                    continue
                quality = _oriented_quality(
                    reads[candidate.read_index], candidate.sequence, position
                )
                if quality is not None and quality >= minimum_base_quality:
                    supporting += 1
            if competing:
                conflict_rejections += 1
                valid = False
            elif supporting < minimum_support:
                quality_rejections += 1
                valid = False
        if valid:
            accepted.append(choice)
            accepted_sides += len(choice_sides)
    return accepted, accepted_sides, conflict_rejections, quality_rejections


def _audit_rejected_extension_sides(
    target: str,
    candidates: list[_Alignment],
    reads: list[Read],
    diagnostics: _RejectedExtensionDiagnostics,
    *,
    minimum_base_quality: int,
    damage_end_window: int,
) -> None:
    """Summarize evidence on candidate sides that failed extension."""
    for side in ("left", "right"):
        boundary = -1 if side == "left" else len(target)
        side_candidates = [
            candidate
            for candidate in candidates
            if (
                candidate.offset < 0
                if side == "left"
                else candidate.offset + len(candidate.sequence) > len(target)
            )
        ]
        if not side_candidates:
            continue
        if side == "left":
            diagnostics.left_sides += 1
        else:
            diagnostics.right_sides += 1

        support_by_base = {base: 0 for base in "ACGT"}
        extensions_by_base: dict[str, list[int]] = {
            base: [] for base in "ACGT"
        }
        for candidate in side_candidates:
            candidate_position = boundary - candidate.offset
            base = candidate.sequence[candidate_position]
            if base not in support_by_base:
                continue
            support_by_base[base] += 1
            extension = (
                -candidate.offset
                if side == "left"
                else candidate.offset + len(candidate.sequence) - len(target)
            )
            extensions_by_base[base].append(extension)
            quality = _oriented_quality(
                reads[candidate.read_index],
                candidate.sequence,
                candidate_position,
            )
            if quality is None:
                diagnostics.missing_quality_boundary_observations += 1
            elif quality >= minimum_base_quality:
                diagnostics.high_quality_boundary_observations += 1
            else:
                diagnostics.low_quality_boundary_observations += 1

            overlap_start = max(0, candidate.offset)
            overlap_stop = min(
                len(target), candidate.offset + len(candidate.sequence)
            )
            for target_position in range(overlap_start, overlap_stop):
                aligned_position = target_position - candidate.offset
                target_base = target[target_position]
                candidate_base = candidate.sequence[aligned_position]
                if target_base == candidate_base:
                    continue
                terminal = (
                    target_position < damage_end_window
                    or target_position >= len(target) - damage_end_window
                    or aligned_position < damage_end_window
                    or aligned_position
                    >= len(candidate.sequence) - damage_end_window
                )
                if (
                    terminal
                    and frozenset((target_base, candidate_base))
                    in _DAMAGE_PAIRS
                ):
                    diagnostics.damage_compatible_mismatches += 1
                else:
                    diagnostics.ordinary_mismatches += 1

        ranked = sorted(support_by_base.items(), key=lambda item: item[1], reverse=True)
        winning_base, support = ranked[0]
        runner_up = ranked[1][1]
        if support == 1:
            diagnostics.support_1 += 1
        elif support == 2:
            diagnostics.support_2 += 1
        elif support == 3:
            diagnostics.support_3 += 1
        elif support == 4:
            diagnostics.support_4 += 1
        elif support >= 5:
            diagnostics.support_5_plus += 1
        if runner_up:
            diagnostics.conflicting_sides += 1
        else:
            diagnostics.agreeing_sides += 1

        extension = max(extensions_by_base[winning_base], default=0)
        if extension <= 5:
            diagnostics.extension_1_5 += 1
        elif extension <= 10:
            diagnostics.extension_6_10 += 1
        elif extension <= 20:
            diagnostics.extension_11_20 += 1
        else:
            diagnostics.extension_21_plus += 1


def _extension_support(target: str, candidates: list[_Alignment]) -> int:
    """Return independent candidate support for the stronger contig end."""
    left = {
        candidate.read_index for candidate in candidates if candidate.offset < 0
    }
    right = {
        candidate.read_index
        for candidate in candidates
        if candidate.offset + len(candidate.sequence) > len(target)
    }
    return max(len(left), len(right))


def _evidence_priority_key(
    record: SequenceRecord,
    candidates: list[_Alignment],
) -> tuple[int, int, int, int]:
    """Rank extension targets by raw support, provenance, then length."""
    return (
        _extension_support(record.current.sequence, candidates),
        len(record.contributing_molecules),
        len(record.current.sequence),
        -record.index,
    )


def discover_progressive_raw_clusters(
    pool: ProgressiveSequencePool,
    *,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
    minimum_identity: float = 0.90,
    minimum_ry_identity: float = 0.99,
    minimum_cluster_size: int = 5,
) -> tuple[tuple[ProgressiveRawCluster, ...], ProgressiveClusteringDiagnostics]:
    """Cluster raw representatives without consuming undersized centers.

    Candidate discovery follows the longest-first center policy used by Linclust
    and CarpeDeam. Only members of an admitted cluster become unavailable;
    undersized centers and their candidates remain eligible for later centers.
    """
    if minimum_cluster_size < 2:
        raise ValueError("minimum_cluster_size must be at least 2")
    if anchors_per_read < minimum_anchor_matches:
        raise ValueError("anchors_per_read must cover minimum_anchor_matches")

    active = pool.active_raw
    if not active:
        return (), ProgressiveClusteringDiagnostics()
    reads = [record.current for record in active]
    molecules = [record.molecule_id for record in active]
    global_indices = [record.index for record in active]
    position_bits = max(1, max(len(read.sequence) for read in reads).bit_length())
    target_window = max(len(read.sequence) for read in reads)
    anchors = _anchor_index(
        reads,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    order = sorted(
        range(len(reads)),
        key=lambda index: (-len(reads[index].sequence), global_indices[index]),
    )
    claimed: set[int] = set()
    clusters: list[ProgressiveRawCluster] = []
    centers_with_candidates = clusters_below_minimum_size = 0

    for center_index in order:
        if center_index in claimed:
            continue
        candidates = _candidate_alignments(
            reads[center_index].sequence,
            reads,
            molecules,
            anchors,
            claimed | {center_index},
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
        if candidates:
            centers_with_candidates += 1
        member_indices = {center_index}
        member_indices.update(candidate.read_index for candidate in candidates)
        independent_molecules = {
            molecules[index] for index in member_indices
        }
        if len(independent_molecules) < minimum_cluster_size:
            clusters_below_minimum_size += 1
            continue
        claimed.update(member_indices)
        clusters.append(
            ProgressiveRawCluster(
                center_index=global_indices[center_index],
                member_indices=tuple(
                    sorted(global_indices[index] for index in member_indices)
                ),
                alignments=tuple(
                    replace(
                        candidate,
                        read_index=global_indices[candidate.read_index],
                    )
                    for candidate in candidates
                ),
            )
        )

    return (
        tuple(clusters),
        ProgressiveClusteringDiagnostics(
            candidate_centers=len(order),
            centers_with_candidates=centers_with_candidates,
            clusters=len(clusters),
            clusters_below_minimum_size=clusters_below_minimum_size,
            clustered_reads=len(claimed),
            retained_unclustered_reads=len(active) - len(claimed),
        ),
    )


def extend_progressive_raw_clusters(
    pool: ProgressiveSequencePool,
    clusters: tuple[ProgressiveRawCluster, ...],
    *,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
    minimum_identity: float = 0.90,
    minimum_ry_identity: float = 0.99,
    minimum_consensus_support: int = 5,
    minimum_correction_support: int = 3,
    dominance_ratio: float = 4.0,
    correction_dominance_ratio: float = 2.0,
    minimum_overlap_margin: int = 3,
    minimum_extension_support: int = 1,
    minimum_confidence_margin: float = 0.0,
    damage_mismatch_penalty: float = 0.25,
    damage_end_window: int = 5,
    reciprocal_best_extension: bool = True,
    extension_consensus: bool = True,
    maximum_rounds: int = 3,
) -> tuple[ProgressiveSequencePool, ProgressiveRawExtensionDiagnostics]:
    """Iteratively extend admitted centers using immutable raw evidence."""
    if len(pool.active_raw) != len(pool.records):
        raise ValueError("raw extension requires a pristine raw sequence pool")
    if minimum_consensus_support < 2:
        raise ValueError("minimum_consensus_support must be at least 2")
    if minimum_correction_support < 2:
        raise ValueError("minimum_correction_support must be at least 2")
    if maximum_rounds < 1:
        raise ValueError("maximum_rounds must be at least 1")

    reads = list(pool.raw_evidence)
    molecules = [record.molecule_id for record in pool.records]
    position_bits = max(1, max(len(read.sequence) for read in reads).bit_length())
    target_window = max(len(read.sequence) for read in reads)
    anchors = _anchor_index(
        reads,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    updated = pool
    corrected_centers = extended_centers = extension_rounds = 0
    recruited_reads = added_bases = consumed_reads = 0
    ambiguous_extensions = 0
    reciprocal = _ReciprocalDiagnostics()
    ranking = _RankingDiagnostics()

    claimed = {
        member_index
        for cluster in clusters
        for member_index in cluster.member_indices
    }

    for cluster in clusters:
        center = updated.records[cluster.center_index]
        target = center.current.sequence
        seed = target
        members = list(cluster.alignments)
        cluster_indices = set(cluster.member_indices)

        for round_index in range(maximum_rounds):
            if round_index == 0:
                candidates = list(cluster.alignments)
            else:
                candidates = _candidate_alignments(
                    target,
                    reads,
                    molecules,
                    anchors,
                    claimed,
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

            tentative_members = members if round_index == 0 else members + candidates
            internal = _consensus(
                target,
                tentative_members,
                minimum_extension_support=len(tentative_members) + 2,
                minimum_internal_support=minimum_consensus_support,
                minimum_correction_support=minimum_correction_support,
                correction_dominance_ratio=correction_dominance_ratio,
                dominance_ratio=dominance_ratio,
                damage_end_window=0,
            )
            selected, _, ambiguous = _unique_best_extensions(
                internal.sequence,
                candidates,
                minimum_overlap_margin=minimum_overlap_margin,
                damage_aware_ranking=True,
                minimum_confidence_margin=minimum_confidence_margin,
                damage_mismatch_penalty=damage_mismatch_penalty,
                ranking_damage_end_window=damage_end_window,
                diagnostics=ranking,
            )
            selected = _filter_extension_boundary_support(
                internal.sequence,
                candidates,
                selected,
                minimum_extension_support,
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
                    diagnostics=reciprocal,
                )
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

            extension_rounds += 1
            ambiguous_extensions += ambiguous
            changed = result.sequence != target
            if round_index > 0 and changed:
                recruited = {candidate.read_index for candidate in candidates}
                claimed.update(recruited)
                cluster_indices.update(recruited)
                recruited_reads += len(recruited)
                members = tentative_members
            if not changed:
                break
            if result.left_extension:
                members = [
                    replace(
                        member,
                        offset=member.offset + result.left_extension,
                    )
                    for member in members
                ]
            target = result.sequence

        current = Read(
            f"progressive_contig_{cluster.center_index + 1}",
            target,
        )
        contributing = frozenset(
            molecules[index] for index in cluster_indices
        )
        if len(target) > len(seed):
            updated_center = center.extended(current, contributing)
            extended_centers += 1
            added_bases += len(target) - len(seed)
        else:
            updated_center = center.corrected(current)
            corrected_centers += target != seed
        updated = updated.replace_record(updated_center)
        for member_index in cluster_indices:
            if member_index == cluster.center_index:
                continue
            updated = updated.replace_record(updated.records[member_index].consumed())
            consumed_reads += 1

    return (
        updated,
        ProgressiveRawExtensionDiagnostics(
            clusters=len(clusters),
            corrected_centers=corrected_centers,
            extended_centers=extended_centers,
            extension_rounds=extension_rounds,
            recruited_reads=recruited_reads,
            added_bases=added_bases,
            consumed_reads=consumed_reads,
            ambiguous_extensions=ambiguous_extensions,
            reciprocal_checks=reciprocal.checks,
            reciprocal_rejections=reciprocal.rejections,
        ),
    )


def iterate_progressive_raw_extension(
    pool: ProgressiveSequencePool,
    *,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
    minimum_identity: float = 0.90,
    minimum_ry_identity: float = 0.99,
    minimum_consensus_support: int = 5,
    dominance_ratio: float = 4.0,
    maximum_iterations: int = 3,
    minimum_overlap_margin: int = 3,
    minimum_confidence_margin: float = 0.0,
    damage_mismatch_penalty: float = 0.25,
    damage_end_window: int = 5,
    evidence_priority: bool = False,
    audit_rejected_extensions: bool = False,
    audit_minimum_base_quality: int = 20,
    require_strict_boundary_evidence: bool = False,
) -> tuple[ProgressiveSequencePool, ProgressiveIterationDiagnostics]:
    """Extend persistent derived centers using only unused raw fragments."""
    if maximum_iterations < 1:
        raise ValueError("maximum_iterations must be at least 1")
    if minimum_consensus_support < 2:
        raise ValueError("minimum_consensus_support must be at least 2")

    updated = pool
    iterations = candidate_alignments = extended_centers = added_bases = 0
    consumed_reads = ambiguous_extensions = insufficient_support = 0
    priority_sweeps = priority_reordered = priority_claim_conflicts = 0
    rejected = _RejectedExtensionDiagnostics()
    strict_accepted = strict_conflicts = strict_quality = 0
    converged = False

    for iteration in range(1, maximum_iterations + 1):
        raw_records = updated.active_raw
        derived_records = updated.active_derived
        if not raw_records or not derived_records:
            converged = True
            break
        raw_reads = [record.raw for record in raw_records]
        raw_molecules = [record.molecule_id for record in raw_records]
        raw_global_indices = [record.index for record in raw_records]
        anchors = _anchor_index(
            raw_reads,
            anchor_k,
            0,
            anchors_per_read,
            maximum_anchor_occurrences,
        )
        position_bits = max(
            1, max(len(read.sequence) for read in raw_reads).bit_length()
        )
        target_window = max(len(read.sequence) for read in raw_reads)
        claimed_local: set[int] = set()
        extended_this_round = 0
        longest_first = sorted(
            derived_records,
            key=lambda record: (-len(record.current.sequence), record.index),
        )
        cached_candidates: dict[int, list[_Alignment]] = {}
        if evidence_priority:
            for center in derived_records:
                cached_candidates[center.index] = _candidate_alignments(
                    center.current.sequence,
                    raw_reads,
                    raw_molecules,
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
            order = sorted(
                derived_records,
                key=lambda record: _evidence_priority_key(
                    record, cached_candidates[record.index]
                ),
                reverse=True,
            )
            priority_sweeps += 1
            priority_reordered += sum(
                first.index != second.index
                for first, second in zip(longest_first, order)
            )
        else:
            order = longest_first

        for center in order:
            target = center.current.sequence
            if evidence_priority:
                cached = cached_candidates[center.index]
                candidates = [
                    candidate
                    for candidate in cached
                    if candidate.read_index not in claimed_local
                ]
                priority_claim_conflicts += len(cached) - len(candidates)
            else:
                candidates = _candidate_alignments(
                    target,
                    raw_reads,
                    raw_molecules,
                    anchors,
                    claimed_local,
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
            candidate_alignments += len(candidates)
            if not candidates:
                continue
            internal = _consensus(
                target,
                candidates,
                minimum_extension_support=len(candidates) + 2,
                minimum_internal_support=minimum_consensus_support,
                dominance_ratio=dominance_ratio,
                damage_end_window=0,
            )
            selected, _, ambiguous = _unique_best_extensions(
                internal.sequence,
                candidates,
                minimum_overlap_margin=minimum_overlap_margin,
                damage_aware_ranking=True,
                minimum_confidence_margin=minimum_confidence_margin,
                damage_mismatch_penalty=damage_mismatch_penalty,
                ranking_damage_end_window=damage_end_window,
            )
            ambiguous_extensions += ambiguous
            selected = _filter_extension_boundary_support(
                internal.sequence,
                candidates,
                selected,
                minimum_consensus_support,
            )
            if require_strict_boundary_evidence and selected:
                (
                    selected,
                    accepted_sides,
                    conflict_rejections,
                    quality_rejections,
                ) = _filter_strict_boundary_evidence(
                    internal.sequence,
                    candidates,
                    selected,
                    raw_reads,
                    minimum_support=minimum_consensus_support,
                    minimum_base_quality=audit_minimum_base_quality,
                )
                strict_accepted += accepted_sides
                strict_conflicts += conflict_rejections
                strict_quality += quality_rejections
            if not selected:
                if audit_rejected_extensions:
                    _audit_rejected_extension_sides(
                        internal.sequence,
                        candidates,
                        raw_reads,
                        rejected,
                        minimum_base_quality=audit_minimum_base_quality,
                        damage_end_window=damage_end_window,
                    )
                insufficient_support += 1
                continue

            supporting: dict[int, _Alignment] = {}
            for choice in selected:
                boundary = -1 if choice.offset < 0 else len(internal.sequence)
                base = choice.sequence[boundary - choice.offset]
                for candidate in candidates:
                    if not (
                        candidate.offset
                        <= boundary
                        < candidate.offset + len(candidate.sequence)
                    ):
                        continue
                    if candidate.sequence[boundary - candidate.offset] != base:
                        continue
                    supporting[candidate.read_index] = candidate
            result = _consensus(
                internal.sequence,
                list(supporting.values()),
                minimum_extension_support=minimum_consensus_support,
                dominance_ratio=dominance_ratio,
                damage_end_window=0,
                allow_internal_consensus=False,
            )
            if len(result.sequence) <= len(target):
                if audit_rejected_extensions:
                    _audit_rejected_extension_sides(
                        internal.sequence,
                        candidates,
                        raw_reads,
                        rejected,
                        minimum_base_quality=audit_minimum_base_quality,
                        damage_end_window=damage_end_window,
                    )
                insufficient_support += 1
                continue

            global_supporters = {
                raw_global_indices[index] for index in supporting
            }
            current = Read(center.current.name, result.sequence)
            updated = updated.replace_record(
                center.extended(
                    current,
                    frozenset(
                        updated.records[index].molecule_id
                        for index in global_supporters
                    ),
                )
            )
            for local_index in supporting:
                claimed_local.add(local_index)
            for global_index in global_supporters:
                updated = updated.replace_record(
                    updated.records[global_index].consumed()
                )
            extension = len(result.sequence) - len(target)
            extended_centers += 1
            extended_this_round += 1
            added_bases += extension
            consumed_reads += len(global_supporters)

        iterations = iteration
        if not extended_this_round:
            converged = True
            break

    projected_lengths = [
        len(record.current.sequence) for record in updated.active_derived
    ]
    return (
        updated,
        ProgressiveIterationDiagnostics(
            iterations=iterations,
            converged=converged,
            candidate_alignments=candidate_alignments,
            extended_centers=extended_centers,
            added_bases=added_bases,
            consumed_reads=consumed_reads,
            ambiguous_extensions=ambiguous_extensions,
            insufficient_consensus_support=insufficient_support,
            evidence_priority_sweeps=priority_sweeps,
            evidence_priority_reordered_centers=priority_reordered,
            evidence_priority_claim_conflicts=priority_claim_conflicts,
            rejected_left_sides=rejected.left_sides,
            rejected_right_sides=rejected.right_sides,
            rejected_support_1=rejected.support_1,
            rejected_support_2=rejected.support_2,
            rejected_support_3=rejected.support_3,
            rejected_support_4=rejected.support_4,
            rejected_support_5_plus=rejected.support_5_plus,
            rejected_agreeing_sides=rejected.agreeing_sides,
            rejected_conflicting_sides=rejected.conflicting_sides,
            rejected_extension_1_5=rejected.extension_1_5,
            rejected_extension_6_10=rejected.extension_6_10,
            rejected_extension_11_20=rejected.extension_11_20,
            rejected_extension_21_plus=rejected.extension_21_plus,
            rejected_high_quality_boundary_observations=(
                rejected.high_quality_boundary_observations
            ),
            rejected_low_quality_boundary_observations=(
                rejected.low_quality_boundary_observations
            ),
            rejected_missing_quality_boundary_observations=(
                rejected.missing_quality_boundary_observations
            ),
            rejected_damage_compatible_mismatches=(
                rejected.damage_compatible_mismatches
            ),
            rejected_ordinary_mismatches=rejected.ordinary_mismatches,
            strict_boundary_accepted_sides=strict_accepted,
            strict_boundary_conflict_rejections=strict_conflicts,
            strict_boundary_quality_rejections=strict_quality,
            projected_contigs=len(projected_lengths),
            projected_bases=sum(projected_lengths),
            projected_n50=_n50(projected_lengths),
            projected_longest_contig=max(projected_lengths, default=0),
        ),
    )
