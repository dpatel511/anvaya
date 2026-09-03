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
        if len(member_indices) < minimum_cluster_size:
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
) -> tuple[ProgressiveSequencePool, ProgressiveIterationDiagnostics]:
    """Extend persistent derived centers using only unused raw fragments."""
    if maximum_iterations < 1:
        raise ValueError("maximum_iterations must be at least 1")
    if minimum_consensus_support < 2:
        raise ValueError("minimum_consensus_support must be at least 2")

    updated = pool
    iterations = candidate_alignments = extended_centers = added_bases = 0
    consumed_reads = ambiguous_extensions = insufficient_support = 0
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
        order = sorted(
            derived_records,
            key=lambda record: (-len(record.current.sequence), record.index),
        )

        for center in order:
            target = center.current.sequence
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
            if not selected:
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
        ),
    )
