"""Raw-evidence-gated links between exhausted progressive contigs."""

from collections import defaultdict
from dataclasses import dataclass

from anvaya.damage_consensus import _anchor_index
from anvaya.overlap_assembly import (
    _Alignment,
    _MasterOverlapEdge,
    _candidate_alignments,
    _n50,
)
from anvaya.overlap_graph import _confirm_raw_mismatches, _project_master_edges
from anvaya.overlap_progressive import ProgressiveSequencePool
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


_Node = tuple[int, bool]
_EdgeKey = tuple[_Node, _Node]


@dataclass(frozen=True, slots=True)
class ProgressiveLinkDiagnostics:
    """Evidence and projection counts for progressive-contig links."""

    input_contigs: int = 0
    candidate_dovetails: int = 0
    exact_dovetails: int = 0
    near_exact_candidates: int = 0
    near_exact_raw_confirmed: int = 0
    near_exact_accepted: int = 0
    mismatch_positions: int = 0
    insufficient_raw_support: int = 0
    strain_conflicts: int = 0
    ambiguous_near_molecules: int = 0
    near_exact_unique_support_rejections: int = 0
    exact_preferred: int = 0
    bridge_candidate_molecules: int = 0
    ambiguous_bridge_molecules: int = 0
    support_at_least_1: int = 0
    support_at_least_2: int = 0
    support_at_least_3: int = 0
    support_at_least_5: int = 0
    max_read_support: int = 0
    supported_dovetails: int = 0
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


def _reverse_node(node: _Node) -> _Node:
    return node[0], not node[1]


def _physical_key(source: _Node, target: _Node) -> _EdgeKey:
    reverse = (_reverse_node(target), _reverse_node(source))
    return min((source, target), reverse)


def _oriented_sequence(contigs: list[Read], node: _Node) -> str:
    sequence = contigs[node[0]].sequence
    return reverse_complement(sequence) if node[1] else sequence


def _add_bidirected_edge(
    destination: dict[_EdgeKey, _MasterOverlapEdge],
    contigs: list[Read],
    edge: _MasterOverlapEdge,
) -> None:
    destination[(edge.source, edge.target)] = edge
    reverse_source = _reverse_node(edge.target)
    reverse_target = _reverse_node(edge.source)
    second_sequence = _oriented_sequence(contigs, edge.target)
    reverse_corrections = tuple(
        (
            len(second_sequence) - 1 - (position - edge.shift),
            reverse_complement(base),
        )
        for position, base in edge.corrections
    )
    destination[(reverse_source, reverse_target)] = _MasterOverlapEdge(
        reverse_source,
        reverse_target,
        len(_oriented_sequence(contigs, reverse_source)) - edge.overlap,
        edge.overlap,
        reverse_corrections,
    )


def _spanning_molecules(
    edge: _MasterOverlapEdge,
    contigs: list[Read],
    raw_reads: list[Read],
    molecule_ids: list[int],
    raw_anchors: dict[int, list[int]],
    *,
    anchor_k: int,
    anchors_per_read: int,
    maximum_anchor_occurrences: int,
    minimum_anchor_matches: int,
    minimum_overlap: int,
    minimum_identity: float,
    minimum_ry_identity: float,
    position_bits: int,
    maximum_read_length: int,
) -> set[int]:
    """Return raw molecules spanning unique sequence on both sides of an edge."""
    left = _oriented_sequence(contigs, edge.source)
    right = _oriented_sequence(contigs, edge.target)
    merged = left + right[edge.overlap :]
    overlap_start = len(left) - edge.overlap
    overlap_end = len(left)
    flank = maximum_read_length - 1
    window_start = max(0, overlap_start - flank)
    window_end = min(len(merged), overlap_end + flank)
    junction = merged[window_start:window_end]
    local_start = overlap_start - window_start
    local_end = overlap_end - window_start
    candidates = _candidate_alignments(
        junction,
        raw_reads,
        molecule_ids,
        raw_anchors,
        set(),
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_overlap=minimum_overlap,
        minimum_identity=minimum_identity,
        minimum_ry_identity=minimum_ry_identity,
        position_bits=position_bits,
        target_window=len(junction),
    )
    return {
        molecule_ids[candidate.read_index]
        for candidate in candidates
        if candidate.offset < local_start
        and candidate.offset + len(candidate.sequence) > local_end
    }


def audit_raw_supported_progressive_links(
    pool: ProgressiveSequencePool,
    *,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
    minimum_identity: float = 0.90,
    minimum_ry_identity: float = 0.99,
    minimum_read_support: int = 2,
    allow_near_exact: bool = False,
    near_exact_minimum_identity: float = 0.98,
    near_exact_maximum_mismatches: int = 1,
    minimum_primary_allele_support: int = 3,
    minimum_alternate_allele_support: int = 2,
    minimum_support_margin: int = 2,
    minimum_base_quality: int = 20,
    damage_end_window: int = 5,
) -> tuple[list[Read], ProgressiveLinkDiagnostics]:
    """Project evidence-gated exact and optionally near-exact contig links.

    The progressive pool remains unchanged. A raw molecule counts only when it
    crosses unique sequence on both sides of exactly one physical dovetail.
    Near-exact links instead require an unambiguous Q20 raw allele at every
    mismatch. Exact links take precedence at their oriented contig ends.
    Branches, cycles, and non-reciprocal graph edges remain separate contigs.
    """
    if minimum_overlap < anchor_k:
        raise ValueError("progressive link overlap must be at least anchor k")
    if minimum_read_support < 1:
        raise ValueError("minimum read support must be at least 1")
    if not 0.0 <= near_exact_minimum_identity <= 1.0:
        raise ValueError("near-exact identity must be between 0 and 1")
    if near_exact_maximum_mismatches < 1:
        raise ValueError("near-exact maximum mismatches must be at least 1")
    if minimum_primary_allele_support < 1:
        raise ValueError("primary allele support must be at least 1")
    if minimum_alternate_allele_support < 1:
        raise ValueError("alternate allele support must be at least 1")
    if minimum_support_margin < 1:
        raise ValueError("support margin must be at least 1")
    if minimum_base_quality < 0:
        raise ValueError("minimum base quality must not be negative")
    if damage_end_window < 0:
        raise ValueError("damage end window must not be negative")

    records = pool.active_derived
    contigs = [record.current for record in records]
    if len(contigs) < 2:
        lengths = [len(contig.sequence) for contig in contigs]
        return contigs, ProgressiveLinkDiagnostics(
            input_contigs=len(contigs),
            projected_contigs=len(contigs),
            projected_bases=sum(lengths),
            projected_n50=_n50(lengths),
            projected_longest_contig=max(lengths, default=0),
        )

    contig_position_bits = max(
        1, max(len(contig.sequence) for contig in contigs).bit_length()
    )
    contig_window = max(len(contig.sequence) for contig in contigs)
    contig_anchors = _anchor_index(
        contigs,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    contig_ids = list(range(len(contigs)))
    physical_edges: dict[_EdgeKey, _MasterOverlapEdge] = {}
    near_edges: dict[_EdgeKey, _MasterOverlapEdge] = {}
    near_molecules: dict[_EdgeKey, set[int]] = {}
    candidate_dovetails = 0
    near_candidates = mismatch_positions = 0
    insufficient_support = strain_conflicts = 0

    raw_reads = list(pool.raw_evidence)
    molecule_ids = [record.molecule_id for record in pool.records]
    raw_position_bits = max(
        1, max(len(read.sequence) for read in raw_reads).bit_length()
    )
    maximum_read_length = max(len(read.sequence) for read in raw_reads)
    raw_anchors = _anchor_index(
        raw_reads,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    raw_candidates_by_source: dict[int, list[_Alignment]] = {}

    for source_index, source in enumerate(contigs):
        candidates = _candidate_alignments(
            source.sequence,
            contigs,
            contig_ids,
            contig_anchors,
            {source_index},
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
            position_bits=contig_position_bits,
            target_window=contig_window,
        )
        for candidate in candidates:
            left_extension = max(0, -candidate.offset)
            right_extension = max(
                0,
                candidate.offset + len(candidate.sequence) - len(source.sequence),
            )
            if bool(left_extension) == bool(right_extension):
                continue
            candidate_dovetails += 1
            target_node = (
                candidate.read_index,
                candidate.sequence != contigs[candidate.read_index].sequence,
            )
            source_node = (source_index, False)
            if right_extension:
                first, second = source_node, target_node
                overlap = len(source.sequence) - candidate.offset
            else:
                first, second = target_node, source_node
                overlap = candidate.offset + len(candidate.sequence)
            edge = _MasterOverlapEdge(
                first,
                second,
                len(_oriented_sequence(contigs, first)) - overlap,
                overlap,
            )
            key = _physical_key(first, second)
            start = max(0, candidate.offset)
            stop = min(
                len(source.sequence),
                candidate.offset + len(candidate.sequence),
            )
            mismatches = [
                position
                for position in range(start, stop)
                if source.sequence[position]
                != candidate.sequence[position - candidate.offset]
            ]
            if not mismatches:
                current = physical_edges.get(key)
                if current is None or edge.overlap > current.overlap:
                    physical_edges[key] = edge
                continue

            near_candidates += 1
            mismatch_positions += len(mismatches)
            if len(mismatches) > near_exact_maximum_mismatches:
                insufficient_support += 1
                continue
            raw_candidates = raw_candidates_by_source.get(source_index)
            if raw_candidates is None:
                raw_candidates = _candidate_alignments(
                    source.sequence,
                    raw_reads,
                    molecule_ids,
                    raw_anchors,
                    set(),
                    anchor_k=anchor_k,
                    anchors_per_read=anchors_per_read,
                    maximum_anchor_occurrences=maximum_anchor_occurrences,
                    minimum_anchor_matches=minimum_anchor_matches,
                    minimum_overlap=anchor_k,
                    minimum_identity=minimum_identity,
                    minimum_ry_identity=minimum_identity,
                    position_bits=raw_position_bits,
                    target_window=len(source.sequence),
                )
                raw_candidates_by_source[source_index] = raw_candidates
            confirmation = _confirm_raw_mismatches(
                source.sequence,
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
            confirmed = _MasterOverlapEdge(
                edge.source,
                edge.target,
                edge.shift,
                edge.overlap,
                corrections,
            )
            current = near_edges.get(key)
            if current is None or confirmed.overlap > current.overlap:
                near_edges[key] = confirmed
                near_molecules[key] = set(confirmation.supporting_molecules)
            elif confirmed.overlap == current.overlap:
                near_molecules[key].update(confirmation.supporting_molecules)

    molecules_by_edge: dict[_EdgeKey, set[int]] = {}
    edges_by_molecule: dict[int, set[_EdgeKey]] = defaultdict(set)
    for key, edge in physical_edges.items():
        molecules = _spanning_molecules(
            edge,
            contigs,
            raw_reads,
            molecule_ids,
            raw_anchors,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=minimum_overlap,
            minimum_identity=minimum_identity,
            minimum_ry_identity=minimum_ry_identity,
            position_bits=raw_position_bits,
            maximum_read_length=maximum_read_length,
        )
        molecules_by_edge[key] = molecules
        for molecule in molecules:
            edges_by_molecule[molecule].add(key)

    unique_support = {
        key: sum(len(edges_by_molecule[molecule]) == 1 for molecule in molecules)
        for key, molecules in molecules_by_edge.items()
    }
    supported = {
        key: physical_edges[key]
        for key, support in unique_support.items()
        if support >= minimum_read_support
    }
    exact_bidirected: dict[_EdgeKey, _MasterOverlapEdge] = {}
    for edge in supported.values():
        _add_bidirected_edge(exact_bidirected, contigs, edge)

    near_edges_by_molecule: dict[int, set[_EdgeKey]] = defaultdict(set)
    for key, molecules in near_molecules.items():
        for molecule in molecules:
            near_edges_by_molecule[molecule].add(key)
    unique_near_support = {
        key: {
            molecule
            for molecule in molecules
            if len(near_edges_by_molecule[molecule]) == 1
        }
        for key, molecules in near_molecules.items()
    }
    unique_near_edges = {
        key: near_edges[key]
        for key, molecules in unique_near_support.items()
        if len(molecules) >= minimum_primary_allele_support
    }
    near_bidirected: dict[_EdgeKey, _MasterOverlapEdge] = {}
    for edge in unique_near_edges.values():
        _add_bidirected_edge(near_bidirected, contigs, edge)

    exact_outgoing = {edge.source for edge in exact_bidirected.values()}
    exact_incoming = {edge.target for edge in exact_bidirected.values()}
    preferred_near = {
        key: edge
        for key, edge in near_bidirected.items()
        if edge.source not in exact_outgoing and edge.target not in exact_incoming
    }
    exact_preferred = (len(near_bidirected) - len(preferred_near)) // 2
    combined = dict(exact_bidirected)
    combined.update(preferred_near)

    projection, graph = _project_master_edges(
        contigs,
        combined,
        name_prefix="progressive_link_contig",
    )
    lengths = [len(contig.sequence) for contig in projection]
    support_values = list(unique_support.values())
    return projection, ProgressiveLinkDiagnostics(
        input_contigs=len(contigs),
        candidate_dovetails=candidate_dovetails,
        exact_dovetails=len(physical_edges),
        near_exact_candidates=near_candidates // 2,
        near_exact_raw_confirmed=len(near_edges),
        near_exact_accepted=len(preferred_near) // 2,
        mismatch_positions=mismatch_positions // 2,
        insufficient_raw_support=insufficient_support // 2,
        strain_conflicts=strain_conflicts // 2,
        ambiguous_near_molecules=sum(
            len(edges) > 1 for edges in near_edges_by_molecule.values()
        ),
        near_exact_unique_support_rejections=(
            len(near_edges) - len(unique_near_edges)
        ),
        exact_preferred=exact_preferred,
        bridge_candidate_molecules=len(edges_by_molecule),
        ambiguous_bridge_molecules=sum(
            len(edges) > 1 for edges in edges_by_molecule.values()
        ),
        support_at_least_1=sum(value >= 1 for value in support_values),
        support_at_least_2=sum(value >= 2 for value in support_values),
        support_at_least_3=sum(value >= 3 for value in support_values),
        support_at_least_5=sum(value >= 5 for value in support_values),
        max_read_support=max(support_values, default=0),
        supported_dovetails=len(supported),
        transitive_edges_removed=graph.transitive_edges_removed,
        ambiguous_ends=graph.ambiguous_ends,
        reciprocal_edges=graph.reciprocal_edges,
        linear_paths=graph.linear_paths,
        cyclic_components=graph.cyclic_components,
        merged_contigs=graph.merged_contigs,
        corrected_overlap_bases=graph.corrected_overlap_bases,
        added_bases=graph.added_bases,
        projected_contigs=len(projection),
        projected_bases=sum(lengths),
        projected_n50=_n50(lengths),
        projected_longest_contig=max(lengths, default=0),
    )
