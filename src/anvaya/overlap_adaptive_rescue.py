"""Low-support rescue constrained by trusted contigs and raw evidence."""

from collections import defaultdict
from dataclasses import dataclass

from anvaya.damage_consensus import _anchor_index
from anvaya.overlap_assembly import (
    _Alignment,
    _MasterOverlapEdge,
    _candidate_alignments,
    _n50,
    project_two_tier_rescue,
)
from anvaya.overlap_graph import _confirm_raw_mismatches, _project_master_edges
from anvaya.overlap_progressive import (
    ProgressiveRawCluster,
    ProgressiveSequencePool,
    _oriented_quality,
    discover_progressive_raw_clusters,
    extend_progressive_raw_clusters,
)
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


@dataclass(frozen=True, slots=True)
class AdaptiveRescueDiagnostics:
    """Admission and projection counts for the low-support rescue tier."""

    primary_contigs: int = 0
    eligible_raw_reads: int = 0
    rescue_clusters: int = 0
    rescue_clustered_reads: int = 0
    rescue_contigs: int = 0
    candidate_primary_links: int = 0
    candidate_rescue_links: int = 0
    exact_primary_links: int = 0
    exact_rescue_links: int = 0
    near_exact_candidates: int = 0
    near_exact_raw_confirmed: int = 0
    near_exact_accepted: int = 0
    mismatch_positions: int = 0
    insufficient_raw_support: int = 0
    strain_conflicts: int = 0
    ambiguous_near_molecules: int = 0
    near_exact_unique_support_rejections: int = 0
    exact_preferred: int = 0
    anchored_components: int = 0
    rejected_unanchored_components: int = 0
    promoted_rescue_contigs: int = 0
    rejected_unattached_contigs: int = 0
    ambiguous_ends: int = 0
    reciprocal_edges: int = 0
    corrected_overlap_bases: int = 0
    projected_contigs: int = 0
    projected_bases: int = 0
    projected_n50: int = 0
    projected_longest_contig: int = 0


@dataclass(frozen=True, slots=True)
class SelectiveRescueDiagnostics:
    """Counts for independently supported, nonredundant rescue contigs."""

    primary_contigs: int = 0
    eligible_raw_reads: int = 0
    rescue_clusters: int = 0
    rescue_clustered_reads: int = 0
    rescue_contigs: int = 0
    insufficient_molecule_support: int = 0
    contained_by_primary: int = 0
    redundant_with_rescue: int = 0
    primary_extensions: int = 0
    ambiguous_primary_extensions: int = 0
    novel_contigs: int = 0
    novel_bases: int = 0
    projected_contigs: int = 0
    projected_bases: int = 0
    projected_n50: int = 0
    projected_longest_contig: int = 0


@dataclass(frozen=True, slots=True)
class SupportTwoRescueDiagnostics:
    """Counts for the high-confidence two-molecule rescue tier."""

    eligible_raw_reads: int = 0
    candidate_clusters: int = 0
    admitted_clusters: int = 0
    ordinary_mismatch_rejections: int = 0
    internal_damage_rejections: int = 0
    boundary_quality_rejections: int = 0
    rescue_contigs: int = 0
    contained_by_primary: int = 0
    redundant_with_rescue: int = 0
    novel_contigs: int = 0
    novel_bases: int = 0
    projected_contigs: int = 0
    projected_bases: int = 0
    projected_n50: int = 0
    projected_longest_contig: int = 0


def _oriented_sequence(contigs: list[Read], node: tuple[int, bool]) -> str:
    sequence = contigs[node[0]].sequence
    return reverse_complement(sequence) if node[1] else sequence


def _add_bidirected_edge(
    edges: dict[tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge],
    contigs: list[Read],
    first: tuple[int, bool],
    second: tuple[int, bool],
    overlap: int,
    corrections: tuple[tuple[int, str], ...] = (),
) -> None:
    edge = _MasterOverlapEdge(
        first,
        second,
        len(_oriented_sequence(contigs, first)) - overlap,
        overlap,
        corrections,
    )
    edges[(first, second)] = edge
    reverse_first = (second[0], not second[1])
    reverse_second = (first[0], not first[1])
    second_sequence = _oriented_sequence(contigs, second)
    reverse_corrections = tuple(
        (
            len(second_sequence) - 1 - (
                position - (len(_oriented_sequence(contigs, first)) - overlap)
            ),
            reverse_complement(base),
        )
        for position, base in corrections
    )
    edges[(reverse_first, reverse_second)] = _MasterOverlapEdge(
        reverse_first,
        reverse_second,
        len(_oriented_sequence(contigs, reverse_first)) - overlap,
        overlap,
        reverse_corrections,
    )


def _reverse_node(node: tuple[int, bool]) -> tuple[int, bool]:
    return node[0], not node[1]


def _physical_key(
    first: tuple[int, bool],
    second: tuple[int, bool],
) -> tuple[tuple[int, bool], tuple[int, bool]]:
    reverse = (_reverse_node(second), _reverse_node(first))
    return min((first, second), reverse)


def _retain_primary_anchored_edges(
    edges: dict[tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge],
    primary_count: int,
) -> tuple[
    dict[tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge],
    int,
    int,
]:
    """Keep physical overlap components containing a trusted primary contig."""
    neighbors: dict[int, set[int]] = defaultdict(set)
    for edge in edges.values():
        source = edge.source[0]
        target = edge.target[0]
        neighbors[source].add(target)
        neighbors[target].add(source)

    anchored_nodes: set[int] = set()
    visited: set[int] = set()
    anchored_components = rejected_components = 0
    for start in sorted(neighbors):
        if start in visited:
            continue
        component: set[int] = set()
        pending = [start]
        while pending:
            node = pending.pop()
            if node in component:
                continue
            component.add(node)
            pending.extend(neighbors[node] - component)
        visited.update(component)
        if any(node < primary_count for node in component):
            anchored_nodes.update(component)
            anchored_components += 1
        else:
            rejected_components += 1

    return (
        {
            key: edge
            for key, edge in edges.items()
            if edge.source[0] in anchored_nodes
            and edge.target[0] in anchored_nodes
        },
        anchored_components,
        rejected_components,
    )


def _assemble_support_rescue_contigs(
    pool: ProgressiveSequencePool,
    *,
    name_prefix: str,
    enforce_molecule_support: bool,
    minimum_rescue_support: int,
    anchor_k: int,
    anchors_per_read: int,
    maximum_anchor_occurrences: int,
    minimum_anchor_matches: int,
    minimum_overlap: int,
) -> tuple[list[Read], int, int, int]:
    """Assemble unused raw molecules and retain provenance-backed contigs."""
    raw_records = list(pool.active_raw)
    if not raw_records:
        return [], 0, 0, 0
    rescue_pool = ProgressiveSequencePool.from_reads(
        [record.raw for record in raw_records],
        [record.molecule_id for record in raw_records],
    )
    clusters, clustering = discover_progressive_raw_clusters(
        rescue_pool,
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_overlap=minimum_overlap,
        minimum_cluster_size=minimum_rescue_support,
    )
    rescue_pool, _ = extend_progressive_raw_clusters(
        rescue_pool,
        clusters,
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_overlap=minimum_overlap,
        minimum_consensus_support=minimum_rescue_support,
        minimum_correction_support=minimum_rescue_support,
        minimum_extension_support=minimum_rescue_support - 1,
        reciprocal_best_extension=True,
    )
    records = list(rescue_pool.active_derived)
    cluster_support = {
        cluster.center_index: len(
            {
                rescue_pool.records[index].molecule_id
                for index in cluster.member_indices
            }
        )
        for cluster in clusters
    }
    accepted = [
        record
        for record in records
        if not enforce_molecule_support
        or cluster_support.get(record.index, 0) >= minimum_rescue_support
    ]
    return (
        [
            Read(f"{name_prefix}_{index + 1}", record.current.sequence)
            for index, record in enumerate(accepted)
        ],
        clustering.clusters,
        clustering.clustered_reads,
        len(records) - len(accepted),
    )


def project_selective_support_rescue(
    pool: ProgressiveSequencePool,
    *,
    minimum_rescue_support: int = 3,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
) -> tuple[list[Read], SelectiveRescueDiagnostics]:
    """Project novel support-three contigs without collapsing strain variants."""
    if minimum_rescue_support < 3:
        raise ValueError("selective rescue support must be at least 3")
    primary = [record.current for record in pool.active_derived]
    rescue, clusters, clustered_reads, insufficient = (
        _assemble_support_rescue_contigs(
            pool,
            name_prefix="selective_rescue",
            enforce_molecule_support=True,
            minimum_rescue_support=minimum_rescue_support,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=minimum_overlap,
        )
    )
    projected, redundancy = project_two_tier_rescue(
        primary,
        rescue,
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_identity=1.0,
        minimum_ry_identity=1.0,
        minimum_coverage=1.0,
    )
    return projected, SelectiveRescueDiagnostics(
        primary_contigs=len(primary),
        eligible_raw_reads=len(pool.active_raw),
        rescue_clusters=clusters,
        rescue_clustered_reads=clustered_reads,
        rescue_contigs=len(rescue),
        insufficient_molecule_support=insufficient,
        contained_by_primary=redundancy.contained_by_primary,
        redundant_with_rescue=redundancy.redundant_with_rescue,
        primary_extensions=redundancy.extends_primary,
        ambiguous_primary_extensions=(
            redundancy.ambiguous_primary_extensions
        ),
        novel_contigs=redundancy.novel_contigs,
        novel_bases=redundancy.novel_bases,
        projected_contigs=redundancy.projected_contigs,
        projected_bases=redundancy.projected_bases,
        projected_n50=redundancy.projected_n50,
        projected_longest_contig=redundancy.projected_longest_contig,
    )


def _support_two_cluster_rejection(
    pool: ProgressiveSequencePool,
    cluster: ProgressiveRawCluster,
    *,
    minimum_base_quality: int,
    damage_end_window: int,
) -> str | None:
    center = pool.records[cluster.center_index].current.sequence
    for alignment in cluster.alignments:
        candidate = alignment.sequence
        start = max(0, alignment.offset)
        stop = min(len(center), alignment.offset + len(candidate))
        candidate_start = start - alignment.offset
        for position in range(start, stop):
            candidate_position = candidate_start + position - start
            left = center[position]
            right = candidate[candidate_position]
            if left == right:
                continue
            if (left in "AG") != (right in "AG"):
                return "ordinary_mismatch"
            terminal = (
                position < damage_end_window
                or len(center) - position <= damage_end_window
                or candidate_position < damage_end_window
                or len(candidate) - candidate_position <= damage_end_window
            )
            if not terminal:
                return "internal_damage"

        extension_positions = range(0, candidate_start)
        if stop - alignment.offset < len(candidate):
            extension_positions = tuple(extension_positions) + tuple(
                range(stop - alignment.offset, len(candidate))
            )
        read = pool.records[alignment.read_index].current
        for position in extension_positions:
            quality = _oriented_quality(read, candidate, position)
            if quality is None or quality < minimum_base_quality:
                return "boundary_quality"
    return None


def project_high_confidence_support_two_rescue(
    pool: ProgressiveSequencePool,
    primary: list[Read],
    *,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 20,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
    minimum_base_quality: int = 20,
    damage_end_window: int = 5,
) -> tuple[list[Read], SupportTwoRescueDiagnostics]:
    """Project a strictly filtered two-molecule rescue tier."""
    if minimum_anchor_matches < 2:
        raise ValueError("support-two rescue requires at least two anchor matches")
    if minimum_base_quality < 0:
        raise ValueError("support-two rescue base quality must not be negative")
    raw_records = list(pool.active_raw)
    if not raw_records:
        return list(primary), SupportTwoRescueDiagnostics(
            projected_contigs=len(primary),
            projected_bases=sum(len(contig.sequence) for contig in primary),
            projected_n50=_n50([len(contig.sequence) for contig in primary]),
            projected_longest_contig=max(
                (len(contig.sequence) for contig in primary), default=0
            ),
        )

    rescue_pool = ProgressiveSequencePool.from_reads(
        [record.raw for record in raw_records],
        [record.molecule_id for record in raw_records],
    )
    clusters, _ = discover_progressive_raw_clusters(
        rescue_pool,
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_overlap=minimum_overlap,
        minimum_cluster_size=2,
    )
    admitted: list[ProgressiveRawCluster] = []
    ordinary = internal_damage = boundary_quality = 0
    for cluster in clusters:
        rejection = _support_two_cluster_rejection(
            rescue_pool,
            cluster,
            minimum_base_quality=minimum_base_quality,
            damage_end_window=damage_end_window,
        )
        if rejection == "ordinary_mismatch":
            ordinary += 1
        elif rejection == "internal_damage":
            internal_damage += 1
        elif rejection == "boundary_quality":
            boundary_quality += 1
        else:
            admitted.append(cluster)

    rescue_pool, _ = extend_progressive_raw_clusters(
        rescue_pool,
        tuple(admitted),
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_overlap=minimum_overlap,
        minimum_consensus_support=2,
        minimum_correction_support=2,
        minimum_extension_support=1,
        reciprocal_best_extension=True,
    )
    rescue = [
        Read(f"support_two_rescue_{index + 1}", record.current.sequence)
        for index, record in enumerate(rescue_pool.active_derived)
    ]
    projected, redundancy = project_two_tier_rescue(
        primary,
        rescue,
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_identity=1.0,
        minimum_ry_identity=1.0,
        minimum_coverage=1.0,
    )
    return projected, SupportTwoRescueDiagnostics(
        eligible_raw_reads=len(raw_records),
        candidate_clusters=len(clusters),
        admitted_clusters=len(admitted),
        ordinary_mismatch_rejections=ordinary,
        internal_damage_rejections=internal_damage,
        boundary_quality_rejections=boundary_quality,
        rescue_contigs=len(rescue),
        contained_by_primary=redundancy.contained_by_primary,
        redundant_with_rescue=redundancy.redundant_with_rescue,
        novel_contigs=redundancy.novel_contigs,
        novel_bases=redundancy.novel_bases,
        projected_contigs=redundancy.projected_contigs,
        projected_bases=redundancy.projected_bases,
        projected_n50=redundancy.projected_n50,
        projected_longest_contig=redundancy.projected_longest_contig,
    )
def audit_adaptive_support_rescue(
    pool: ProgressiveSequencePool,
    *,
    minimum_rescue_support: int = 3,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
    allow_rescue_chains: bool = False,
    allow_near_exact: bool = False,
    near_exact_minimum_identity: float = 0.98,
    near_exact_maximum_mismatches: int = 1,
    minimum_primary_allele_support: int = 3,
    minimum_alternate_allele_support: int = 2,
    minimum_support_margin: int = 2,
    minimum_base_quality: int = 20,
    damage_end_window: int = 5,
) -> tuple[list[Read], AdaptiveRescueDiagnostics]:
    """Project anchored support-three rescue through conservative dovetails."""
    if minimum_rescue_support < 3:
        raise ValueError("adaptive rescue support must be at least 3")
    if not 0.0 <= near_exact_minimum_identity <= 1.0:
        raise ValueError("adaptive rescue identity must be between 0 and 1")
    if near_exact_maximum_mismatches < 1:
        raise ValueError("adaptive rescue maximum mismatches must be at least 1")
    if minimum_primary_allele_support < 1:
        raise ValueError("adaptive rescue primary support must be at least 1")
    if minimum_alternate_allele_support < 1:
        raise ValueError("adaptive rescue alternate support must be at least 1")
    if minimum_support_margin < 1:
        raise ValueError("adaptive rescue support margin must be at least 1")
    if minimum_base_quality < 0:
        raise ValueError("adaptive rescue base quality must not be negative")
    if damage_end_window < 0:
        raise ValueError("adaptive rescue damage end window must not be negative")
    primary = [record.current for record in pool.active_derived]
    raw_records = list(pool.active_raw)
    if not primary or not raw_records:
        lengths = [len(contig.sequence) for contig in primary]
        return primary, AdaptiveRescueDiagnostics(
            primary_contigs=len(primary),
            eligible_raw_reads=len(raw_records),
            projected_contigs=len(primary),
            projected_bases=sum(lengths),
            projected_n50=_n50(lengths),
            projected_longest_contig=max(lengths, default=0),
        )

    rescue, rescue_clusters, rescue_clustered_reads, _ = (
        _assemble_support_rescue_contigs(
            pool,
            name_prefix="adaptive_rescue",
            enforce_molecule_support=False,
            minimum_rescue_support=minimum_rescue_support,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=minimum_overlap,
        )
    )
    combined = primary + rescue
    if not rescue:
        lengths = [len(contig.sequence) for contig in primary]
        return primary, AdaptiveRescueDiagnostics(
            primary_contigs=len(primary),
            eligible_raw_reads=len(raw_records),
            rescue_clusters=rescue_clusters,
            rescue_clustered_reads=rescue_clustered_reads,
            projected_contigs=len(primary),
            projected_bases=sum(lengths),
            projected_n50=_n50(lengths),
            projected_longest_contig=max(lengths, default=0),
        )

    targets = combined if allow_rescue_chains else primary
    anchors = _anchor_index(
        targets,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    target_ids = list(range(len(targets)))
    position_bits = max(1, max(len(read.sequence) for read in targets).bit_length())
    target_window = max(len(read.sequence) for read in targets)
    exact_physical: dict[
        tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge
    ] = {}
    near_edges: dict[
        tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge
    ] = {}
    near_molecules: dict[
        tuple[tuple[int, bool], tuple[int, bool]], set[int]
    ] = {}
    raw_reads = list(pool.raw_evidence)
    molecule_ids = [record.molecule_id for record in pool.records]
    raw_anchors = _anchor_index(
        raw_reads,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    raw_position_bits = max(
        1, max((len(read.sequence) for read in raw_reads), default=1).bit_length()
    )
    candidate_primary_keys: set[
        tuple[tuple[int, bool], tuple[int, bool]]
    ] = set()
    candidate_rescue_keys: set[
        tuple[tuple[int, bool], tuple[int, bool]]
    ] = set()
    seen_candidates: set[
        tuple[tuple[int, bool], tuple[int, bool]]
    ] = set()
    exact_primary_keys: set[
        tuple[tuple[int, bool], tuple[int, bool]]
    ] = set()
    exact_rescue_keys: set[
        tuple[tuple[int, bool], tuple[int, bool]]
    ] = set()
    near_candidates = mismatch_positions = 0
    insufficient_support = strain_conflicts = 0
    for rescue_offset, item in enumerate(rescue):
        rescue_index = len(primary) + rescue_offset
        candidates = _candidate_alignments(
            item.sequence,
            targets,
            target_ids,
            anchors,
            {rescue_index},
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=minimum_overlap,
            minimum_identity=(
                near_exact_minimum_identity if allow_near_exact else 1.0
            ),
            minimum_ry_identity=(
                near_exact_minimum_identity if allow_near_exact else 1.0
            ),
            position_bits=position_bits,
            target_window=target_window,
        )
        raw_candidates: list[_Alignment] | None = None
        for candidate in candidates:
            left_extension = max(0, -candidate.offset)
            right_extension = max(
                0,
                candidate.offset + len(candidate.sequence) - len(item.sequence),
            )
            if bool(left_extension) == bool(right_extension):
                continue
            target_node = (
                candidate.read_index,
                candidate.sequence != targets[candidate.read_index].sequence,
            )
            rescue_node = (rescue_index, False)
            if right_extension:
                first, second = rescue_node, target_node
                overlap = len(item.sequence) - candidate.offset
            else:
                first, second = target_node, rescue_node
                overlap = candidate.offset + len(candidate.sequence)
            key = _physical_key(first, second)
            if key in seen_candidates:
                continue
            seen_candidates.add(key)
            target_is_primary = candidate.read_index < len(primary)
            if target_is_primary:
                candidate_primary_keys.add(key)
            else:
                candidate_rescue_keys.add(key)
            start = max(0, candidate.offset)
            stop = min(
                len(item.sequence),
                candidate.offset + len(candidate.sequence),
            )
            mismatches = [
                position
                for position in range(start, stop)
                if item.sequence[position]
                != candidate.sequence[position - candidate.offset]
            ]
            if not mismatches:
                edge = _MasterOverlapEdge(
                    first,
                    second,
                    len(_oriented_sequence(combined, first)) - overlap,
                    overlap,
                )
                current = exact_physical.get(key)
                if current is None or edge.overlap > current.overlap:
                    exact_physical[key] = edge
                if target_is_primary:
                    exact_primary_keys.add(key)
                else:
                    exact_rescue_keys.add(key)
                continue

            near_candidates += 1
            mismatch_positions += len(mismatches)
            if len(mismatches) > near_exact_maximum_mismatches:
                insufficient_support += 1
                continue
            if raw_candidates is None:
                raw_candidates = _candidate_alignments(
                    item.sequence,
                    raw_reads,
                    molecule_ids,
                    raw_anchors,
                    set(),
                    anchor_k=anchor_k,
                    anchors_per_read=anchors_per_read,
                    maximum_anchor_occurrences=maximum_anchor_occurrences,
                    minimum_anchor_matches=minimum_anchor_matches,
                    minimum_overlap=anchor_k,
                    minimum_identity=0.90,
                    minimum_ry_identity=0.90,
                    position_bits=raw_position_bits,
                    target_window=len(item.sequence),
                )
            confirmation = _confirm_raw_mismatches(
                item.sequence,
                candidate,
                mismatches,
                raw_candidates,
                raw_reads,
                molecule_ids,
                minimum_primary_allele_support=minimum_primary_allele_support,
                minimum_alternate_allele_support=(
                    minimum_alternate_allele_support
                ),
                minimum_support_margin=minimum_support_margin,
                minimum_base_quality=minimum_base_quality,
                damage_end_window=damage_end_window,
            )
            if confirmation.strain_conflict:
                strain_conflicts += 1
                continue
            if confirmation.insufficient_support:
                insufficient_support += 1
                continue
            corrections = tuple(
                (
                    position
                    if right_extension
                    else position - candidate.offset,
                    winner,
                )
                for position, winner in confirmation.winners
            )
            edge = _MasterOverlapEdge(
                first,
                second,
                len(_oriented_sequence(combined, first)) - overlap,
                overlap,
                corrections,
            )
            current = near_edges.get(key)
            if current is None or edge.overlap > current.overlap:
                near_edges[key] = edge
                near_molecules[key] = set(confirmation.supporting_molecules)
            elif edge.overlap == current.overlap:
                near_molecules[key].update(confirmation.supporting_molecules)

    exact_edges: dict[
        tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge
    ] = {}
    for edge in exact_physical.values():
        _add_bidirected_edge(
            exact_edges,
            combined,
            edge.source,
            edge.target,
            edge.overlap,
        )

    near_edges_by_molecule: dict[
        int, set[tuple[tuple[int, bool], tuple[int, bool]]]
    ] = defaultdict(set)
    for key, molecules in near_molecules.items():
        for molecule in molecules:
            near_edges_by_molecule[molecule].add(key)
    unique_near = {
        key: near_edges[key]
        for key, molecules in near_molecules.items()
        if sum(len(near_edges_by_molecule[molecule]) == 1 for molecule in molecules)
        >= minimum_primary_allele_support
    }
    near_bidirected: dict[
        tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge
    ] = {}
    for edge in unique_near.values():
        _add_bidirected_edge(
            near_bidirected,
            combined,
            edge.source,
            edge.target,
            edge.overlap,
            edge.corrections,
        )

    exact_outgoing = {edge.source for edge in exact_edges.values()}
    exact_incoming = {edge.target for edge in exact_edges.values()}
    preferred_near = {
        key: edge
        for key, edge in near_bidirected.items()
        if edge.source not in exact_outgoing and edge.target not in exact_incoming
    }
    exact_preferred = (len(near_bidirected) - len(preferred_near)) // 2
    edges = dict(exact_edges)
    edges.update(preferred_near)
    edges, anchored_components, rejected_components = (
        _retain_primary_anchored_edges(edges, len(primary))
    )
    accepted_near_keys = {
        _physical_key(edge.source, edge.target)
        for edge in edges.values()
        if _physical_key(edge.source, edge.target) in unique_near
    }

    projection, graph = _project_master_edges(
        combined,
        edges,
        name_prefix="adaptive_support_contig",
    )
    unattached = [
        contig for contig in projection if contig.name.startswith("adaptive_rescue_")
    ]
    projected = [
        contig for contig in projection
        if not contig.name.startswith("adaptive_rescue_")
    ]
    promoted = len(rescue) - len(unattached)
    lengths = [len(contig.sequence) for contig in projected]
    return projected, AdaptiveRescueDiagnostics(
        primary_contigs=len(primary),
        eligible_raw_reads=len(raw_records),
        rescue_clusters=rescue_clusters,
        rescue_clustered_reads=rescue_clustered_reads,
        rescue_contigs=len(rescue),
        candidate_primary_links=len(candidate_primary_keys),
        candidate_rescue_links=len(candidate_rescue_keys),
        exact_primary_links=len(exact_primary_keys),
        exact_rescue_links=len(exact_rescue_keys),
        near_exact_candidates=near_candidates,
        near_exact_raw_confirmed=len(near_edges),
        near_exact_accepted=len(accepted_near_keys),
        mismatch_positions=mismatch_positions,
        insufficient_raw_support=insufficient_support,
        strain_conflicts=strain_conflicts,
        ambiguous_near_molecules=sum(
            len(keys) > 1 for keys in near_edges_by_molecule.values()
        ),
        near_exact_unique_support_rejections=(len(near_edges) - len(unique_near)),
        exact_preferred=exact_preferred,
        anchored_components=anchored_components,
        rejected_unanchored_components=rejected_components,
        promoted_rescue_contigs=promoted,
        rejected_unattached_contigs=len(unattached),
        ambiguous_ends=graph.ambiguous_ends,
        reciprocal_edges=graph.reciprocal_edges,
        corrected_overlap_bases=graph.corrected_overlap_bases,
        projected_contigs=len(projected),
        projected_bases=sum(lengths),
        projected_n50=_n50(lengths),
        projected_longest_contig=max(lengths, default=0),
    )
