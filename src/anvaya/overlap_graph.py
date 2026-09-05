"""Evidence-gated overlap graph construction and path projection."""

from collections import defaultdict
from dataclasses import dataclass

from anvaya.damage_consensus import _anchor_index, _identity
from anvaya.overlap_assembly import (
    MasterOverlapGraphDiagnostics,
    RawConfirmedMasterGraphDiagnostics,
    StrainSafeContainmentDiagnostics,
    _Alignment,
    _DAMAGE_PAIRS,
    _MasterGraphProjection,
    _MasterOverlapEdge,
    _candidate_alignments,
    _n50,
)
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


@dataclass(frozen=True, slots=True)
class _RawMismatchConfirmation:
    winners: tuple[tuple[int, str], ...] = ()
    supporting_molecules: frozenset[int] = frozenset()
    strain_conflict: bool = False
    insufficient_support: bool = False


def _confirm_raw_mismatches(
    source_sequence: str,
    candidate: _Alignment,
    mismatch_positions: list[int],
    raw_candidates: list[_Alignment],
    raw_reads: list[Read],
    molecule_ids: list[int],
    *,
    minimum_primary_allele_support: int,
    minimum_alternate_allele_support: int,
    minimum_support_margin: int,
    minimum_base_quality: int,
    damage_end_window: int,
) -> _RawMismatchConfirmation:
    """Resolve overlap mismatches from independent high-quality raw molecules."""
    winners: list[tuple[int, str]] = []
    supporting_molecules: set[int] | None = None
    for position in mismatch_positions:
        source_base = source_sequence[position]
        candidate_base = candidate.sequence[position - candidate.offset]
        alleles = frozenset((source_base, candidate_base))
        support: dict[str, set[int]] = defaultdict(set)
        for raw in raw_candidates:
            if not raw.offset <= position < raw.offset + len(raw.sequence):
                continue
            raw_position = position - raw.offset
            record = raw_reads[raw.read_index]
            quality_position = (
                raw_position
                if raw.sequence == record.sequence
                else len(record.sequence) - raw_position - 1
            )
            if record.qualities is None:
                continue
            if record.qualities[quality_position] < minimum_base_quality:
                continue
            if alleles in _DAMAGE_PAIRS and (
                quality_position < damage_end_window
                or quality_position >= len(record.sequence) - damage_end_window
            ):
                continue
            base = raw.sequence[raw_position]
            if base in alleles:
                support[base].add(molecule_ids[raw.read_index])
        source_support = support[source_base]
        candidate_support = support[candidate_base]
        if (
            len(source_support) >= minimum_alternate_allele_support
            and len(candidate_support) >= minimum_alternate_allele_support
        ):
            return _RawMismatchConfirmation(strain_conflict=True)
        if (
            len(source_support) >= minimum_primary_allele_support
            and len(source_support) - len(candidate_support) >= minimum_support_margin
        ):
            winner = source_base
            winner_support = source_support
        elif (
            len(candidate_support) >= minimum_primary_allele_support
            and len(candidate_support) - len(source_support) >= minimum_support_margin
        ):
            winner = candidate_base
            winner_support = candidate_support
        else:
            return _RawMismatchConfirmation(insufficient_support=True)
        winners.append((position, winner))
        if supporting_molecules is None:
            supporting_molecules = set(winner_support)
        else:
            supporting_molecules.intersection_update(winner_support)

    return _RawMismatchConfirmation(
        winners=tuple(winners),
        supporting_molecules=frozenset(supporting_molecules or ()),
    )


def audit_master_overlap_graph(
    contigs: list[Read],
    *,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_overlap: int = 30,
) -> tuple[list[Read], MasterOverlapGraphDiagnostics]:
    """Project exact, reciprocal, non-branching paths without changing input."""
    if not contigs:
        return [], MasterOverlapGraphDiagnostics()
    if minimum_overlap < anchor_k:
        raise ValueError("master graph overlap must be at least anchor k")

    def oriented_sequence(node: tuple[int, bool]) -> str:
        sequence = contigs[node[0]].sequence
        return reverse_complement(sequence) if node[1] else sequence

    def add_edge(
        edges: dict[tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge],
        source: tuple[int, bool],
        target: tuple[int, bool],
        shift: int,
        overlap: int,
    ) -> None:
        edge = _MasterOverlapEdge(source, target, shift, overlap)
        key = (source, target)
        current = edges.get(key)
        if current is None or edge.overlap > current.overlap:
            edges[key] = edge

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
    edges: dict[
        tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge
    ] = {}
    candidate_dovetails = 0
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
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            position_bits=position_bits,
            target_window=target_window,
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
            reversed_target = (
                candidate.sequence != contigs[candidate.read_index].sequence
            )
            source_node = (source_index, False)
            candidate_node = (candidate.read_index, reversed_target)
            if right_extension:
                first, second = source_node, candidate_node
                overlap = len(source.sequence) - candidate.offset
            else:
                first, second = candidate_node, source_node
                overlap = candidate.offset + len(candidate.sequence)
            shift = len(oriented_sequence(first)) - overlap
            add_edge(edges, first, second, shift, overlap)
            add_edge(
                edges,
                (second[0], not second[1]),
                (first[0], not first[1]),
                shift,
                overlap,
            )

    outgoing: dict[tuple[int, bool], list[_MasterOverlapEdge]] = defaultdict(list)
    for edge in edges.values():
        outgoing[edge.source].append(edge)
    transitive: set[tuple[tuple[int, bool], tuple[int, bool]]] = set()
    for edge in edges.values():
        for first in outgoing[edge.source]:
            if first.target == edge.target:
                continue
            for second in outgoing.get(first.target, ()):
                if (
                    second.target == edge.target
                    and first.shift + second.shift == edge.shift
                ):
                    transitive.add((edge.source, edge.target))
                    break
            if (edge.source, edge.target) in transitive:
                break
    reduced = {
        key: edge for key, edge in edges.items() if key not in transitive
    }

    outgoing = defaultdict(list)
    incoming: dict[tuple[int, bool], list[_MasterOverlapEdge]] = defaultdict(list)
    for edge in reduced.values():
        outgoing[edge.source].append(edge)
        incoming[edge.target].append(edge)
    ambiguous_nodes = {
        node
        for node in set(outgoing) | set(incoming)
        if len(outgoing[node]) > 1 or len(incoming[node]) > 1
    }
    reciprocal = {
        key: edge
        for key, edge in reduced.items()
        if len(outgoing[edge.source]) == 1 and len(incoming[edge.target]) == 1
    }
    linear_out = {edge.source: edge for edge in reciprocal.values()}
    linear_in = {edge.target: edge for edge in reciprocal.values()}

    projected_paths: dict[str, tuple[str, set[int]]] = {}
    visited_edges: set[tuple[tuple[int, bool], tuple[int, bool]]] = set()
    starts = sorted(node for node in linear_out if node not in linear_in)
    for start in starts:
        path = [start]
        node = start
        while node in linear_out:
            edge = linear_out[node]
            key = (edge.source, edge.target)
            if key in visited_edges or edge.target in path:
                break
            visited_edges.add(key)
            path.append(edge.target)
            node = edge.target
        if len(path) < 2:
            continue
        sequence = oriented_sequence(path[0])
        for node in path[1:]:
            edge = linear_in[node]
            sequence += oriented_sequence(node)[edge.overlap:]
        canonical = min(sequence, reverse_complement(sequence))
        projected_paths.setdefault(
            canonical,
            (sequence, {node[0] for node in path}),
        )

    cyclic_components: set[frozenset[int]] = set()
    for edge in reciprocal.values():
        key = (edge.source, edge.target)
        if key in visited_edges:
            continue
        component: set[int] = set()
        node = edge.source
        while node in linear_out:
            next_edge = linear_out[node]
            next_key = (next_edge.source, next_edge.target)
            if next_key in visited_edges:
                break
            visited_edges.add(next_key)
            component.add(node[0])
            node = next_edge.target
        if component:
            cyclic_components.add(frozenset(component))

    merged: list[Read] = []
    merged_indices: set[int] = set()
    extension_bases = 0
    for index, (_, (sequence, physical_nodes)) in enumerate(
        sorted(projected_paths.items()),
        start=1,
    ):
        merged.append(Read(f"master_overlap_contig_{index}", sequence))
        merged_indices.update(physical_nodes)
        extension_bases += len(sequence) - max(
            len(contigs[node].sequence) for node in physical_nodes
        )
    projection = merged + [
        contig for index, contig in enumerate(contigs) if index not in merged_indices
    ]
    lengths = [len(contig.sequence) for contig in projection]
    return projection, MasterOverlapGraphDiagnostics(
        input_contigs=len(contigs),
        candidate_dovetails=candidate_dovetails,
        exact_dovetails=len(edges) // 2,
        transitive_edges_removed=len(transitive) // 2,
        ambiguous_ends=len(ambiguous_nodes) // 2,
        reciprocal_edges=len(reciprocal) // 2,
        linear_paths=len(projected_paths),
        cyclic_components=len(cyclic_components),
        merged_contigs=len(merged_indices),
        added_bases=extension_bases,
        projected_contigs=len(projection),
        projected_bases=sum(lengths),
        projected_n50=_n50(lengths),
        projected_longest_contig=max(lengths, default=0),
    )

def _project_master_edges(
    contigs: list[Read],
    edges: dict[tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge],
    *,
    name_prefix: str,
) -> tuple[list[Read], _MasterGraphProjection]:
    """Reduce a bidirected overlap graph and spell reciprocal linear paths."""

    def oriented_sequence(node: tuple[int, bool]) -> str:
        sequence = contigs[node[0]].sequence
        return reverse_complement(sequence) if node[1] else sequence

    outgoing: dict[tuple[int, bool], list[_MasterOverlapEdge]] = defaultdict(list)
    for edge in edges.values():
        outgoing[edge.source].append(edge)
    transitive: set[tuple[tuple[int, bool], tuple[int, bool]]] = set()
    for edge in edges.values():
        for first in outgoing[edge.source]:
            if first.target == edge.target:
                continue
            for second in outgoing.get(first.target, ()):
                if (
                    second.target == edge.target
                    and first.shift + second.shift == edge.shift
                ):
                    transitive.add((edge.source, edge.target))
                    break
            if (edge.source, edge.target) in transitive:
                break
    reduced = {
        key: edge for key, edge in edges.items() if key not in transitive
    }

    outgoing = defaultdict(list)
    incoming: dict[tuple[int, bool], list[_MasterOverlapEdge]] = defaultdict(list)
    for edge in reduced.values():
        outgoing[edge.source].append(edge)
        incoming[edge.target].append(edge)
    ambiguous_nodes = {
        node
        for node in set(outgoing) | set(incoming)
        if len(outgoing[node]) > 1 or len(incoming[node]) > 1
    }
    reciprocal = {
        key: edge
        for key, edge in reduced.items()
        if len(outgoing[edge.source]) == 1 and len(incoming[edge.target]) == 1
    }
    linear_out = {edge.source: edge for edge in reciprocal.values()}
    linear_in = {edge.target: edge for edge in reciprocal.values()}

    projected_paths: dict[str, tuple[str, set[int], int]] = {}
    visited_edges: set[tuple[tuple[int, bool], tuple[int, bool]]] = set()
    starts = sorted(node for node in linear_out if node not in linear_in)
    for start in starts:
        path = [start]
        node = start
        while node in linear_out:
            edge = linear_out[node]
            key = (edge.source, edge.target)
            if key in visited_edges or edge.target in path:
                break
            visited_edges.add(key)
            path.append(edge.target)
            node = edge.target
        if len(path) < 2:
            continue
        sequence = list(oriented_sequence(path[0]))
        corrected = 0
        for node in path[1:]:
            edge = linear_in[node]
            source_start = len(sequence) - len(oriented_sequence(edge.source))
            for position, base in edge.corrections:
                assembled_position = source_start + position
                if sequence[assembled_position] != base:
                    sequence[assembled_position] = base
                    corrected += 1
            sequence.extend(oriented_sequence(node)[edge.overlap:])
        assembled = "".join(sequence)
        canonical = min(assembled, reverse_complement(assembled))
        projected_paths.setdefault(
            canonical,
            (assembled, {path_node[0] for path_node in path}, corrected),
        )

    cyclic_components: set[frozenset[int]] = set()
    for edge in reciprocal.values():
        key = (edge.source, edge.target)
        if key in visited_edges:
            continue
        component: set[int] = set()
        node = edge.source
        while node in linear_out:
            next_edge = linear_out[node]
            next_key = (next_edge.source, next_edge.target)
            if next_key in visited_edges:
                break
            visited_edges.add(next_key)
            component.add(node[0])
            node = next_edge.target
        if component:
            cyclic_components.add(frozenset(component))

    merged: list[Read] = []
    merged_indices: set[int] = set()
    corrected_bases = extension_bases = 0
    for index, (_, (sequence, physical_nodes, corrected)) in enumerate(
        sorted(projected_paths.items()),
        start=1,
    ):
        merged.append(Read(f"{name_prefix}_{index}", sequence))
        merged_indices.update(physical_nodes)
        corrected_bases += corrected
        extension_bases += len(sequence) - max(
            len(contigs[node].sequence) for node in physical_nodes
        )
    projection = merged + [
        contig for index, contig in enumerate(contigs) if index not in merged_indices
    ]
    return projection, _MasterGraphProjection(
        transitive_edges_removed=len(transitive) // 2,
        ambiguous_ends=len(ambiguous_nodes) // 2,
        reciprocal_edges=len(reciprocal) // 2,
        linear_paths=len(projected_paths),
        cyclic_components=len(cyclic_components),
        merged_contigs=len(merged_indices),
        corrected_overlap_bases=corrected_bases,
        added_bases=extension_bases,
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
    """Project exact and raw-confirmed near-exact branchless overlaps."""
    if not contigs:
        return [], RawConfirmedMasterGraphDiagnostics()
    if minimum_overlap < anchor_k:
        raise ValueError("raw-confirmed master overlap must be at least anchor k")
    if not 0.0 <= minimum_identity <= 1.0:
        raise ValueError("raw-confirmed master identity must be between 0 and 1")
    if maximum_mismatches < 1:
        raise ValueError("raw-confirmed maximum mismatches must be at least 1")
    if minimum_primary_allele_support < 1:
        raise ValueError("primary allele support must be at least 1")
    if minimum_alternate_allele_support < 1:
        raise ValueError("alternate allele support must be at least 1")
    if minimum_support_margin < 1:
        raise ValueError("raw-confirmed support margin must be at least 1")
    if minimum_base_quality < 0:
        raise ValueError("raw-confirmed base quality must not be negative")
    if damage_end_window < 0:
        raise ValueError("raw-confirmed damage end window must not be negative")

    molecules = list(range(len(raw_reads))) if molecule_ids is None else molecule_ids
    if len(molecules) != len(raw_reads):
        raise ValueError("molecule IDs must align one-to-one with raw reads")

    def oriented_sequence(node: tuple[int, bool]) -> str:
        sequence = contigs[node[0]].sequence
        return reverse_complement(sequence) if node[1] else sequence

    def add_edge(
        destination: dict[
            tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge
        ],
        first: tuple[int, bool],
        second: tuple[int, bool],
        overlap: int,
        corrections: tuple[tuple[int, str], ...] = (),
    ) -> None:
        shift = len(oriented_sequence(first)) - overlap
        edge = _MasterOverlapEdge(first, second, shift, overlap, corrections)
        key = (first, second)
        current = destination.get(key)
        if current is None or edge.overlap > current.overlap:
            destination[key] = edge

        second_sequence = oriented_sequence(second)
        reverse_corrections = tuple(
            (
                len(second_sequence) - 1 - (position - shift),
                reverse_complement(base),
            )
            for position, base in corrections
        )
        reverse_first = (second[0], not second[1])
        reverse_second = (first[0], not first[1])
        reverse_edge = _MasterOverlapEdge(
            reverse_first,
            reverse_second,
            len(oriented_sequence(reverse_first)) - overlap,
            overlap,
            reverse_corrections,
        )
        reverse_key = (reverse_first, reverse_second)
        current = destination.get(reverse_key)
        if current is None or reverse_edge.overlap > current.overlap:
            destination[reverse_key] = reverse_edge

    contig_bits = max(1, max(len(item.sequence) for item in contigs).bit_length())
    raw_bits = max(
        1, max((len(item.sequence) for item in raw_reads), default=1).bit_length()
    )
    target_window = max(len(item.sequence) for item in contigs)
    contig_anchors = _anchor_index(
        contigs,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    raw_anchors = (
        _anchor_index(
            raw_reads,
            anchor_k,
            0,
            anchors_per_read,
            maximum_anchor_occurrences,
        )
        if raw_reads
        else {}
    )
    contig_ids = list(range(len(contigs)))
    exact_edges: dict[
        tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge
    ] = {}
    near_edges: dict[
        tuple[tuple[int, bool], tuple[int, bool]], _MasterOverlapEdge
    ] = {}
    candidate_dovetails = near_candidates = mismatch_positions = 0
    insufficient_support = strain_conflicts = 0

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
            minimum_identity=minimum_identity,
            minimum_ry_identity=minimum_identity,
            position_bits=contig_bits,
            target_window=target_window,
        )
        raw_candidates: list[_Alignment] | None = None
        for candidate in candidates:
            left_extension = max(0, -candidate.offset)
            right_extension = max(
                0,
                candidate.offset + len(candidate.sequence) - len(source.sequence),
            )
            if bool(left_extension) == bool(right_extension):
                continue
            candidate_dovetails += 1
            start = max(0, candidate.offset)
            stop = min(
                len(source.sequence), candidate.offset + len(candidate.sequence)
            )
            mismatches = [
                position
                for position in range(start, stop)
                if source.sequence[position]
                != candidate.sequence[position - candidate.offset]
            ]
            reversed_target = (
                candidate.sequence != contigs[candidate.read_index].sequence
            )
            source_node = (source_index, False)
            candidate_node = (candidate.read_index, reversed_target)
            if right_extension:
                first, second = source_node, candidate_node
                overlap = len(source.sequence) - candidate.offset
            else:
                first, second = candidate_node, source_node
                overlap = candidate.offset + len(candidate.sequence)
            if not mismatches:
                add_edge(exact_edges, first, second, overlap)
                continue
            near_candidates += 1
            mismatch_positions += len(mismatches)
            if len(mismatches) > maximum_mismatches:
                insufficient_support += 1
                continue
            if raw_candidates is None:
                raw_candidates = _candidate_alignments(
                    source.sequence,
                    raw_reads,
                    molecules,
                    raw_anchors,
                    set(),
                    anchor_k=anchor_k,
                    anchors_per_read=anchors_per_read,
                    maximum_anchor_occurrences=maximum_anchor_occurrences,
                    minimum_anchor_matches=minimum_anchor_matches,
                    minimum_overlap=anchor_k,
                    minimum_identity=0.90,
                    minimum_ry_identity=0.90,
                    position_bits=raw_bits,
                    target_window=len(source.sequence),
                )
            confirmation = _confirm_raw_mismatches(
                source.sequence,
                candidate,
                mismatches,
                raw_candidates,
                raw_reads,
                molecules,
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
            add_edge(near_edges, first, second, overlap, corrections)

    exact_outgoing = {edge.source for edge in exact_edges.values()}
    exact_incoming = {edge.target for edge in exact_edges.values()}
    preferred_near_edges = {
        key: edge
        for key, edge in near_edges.items()
        if edge.source not in exact_outgoing and edge.target not in exact_incoming
    }
    exact_preferred = (len(near_edges) - len(preferred_near_edges)) // 2
    combined_edges = dict(exact_edges)
    combined_edges.update(preferred_near_edges)
    projection, graph = _project_master_edges(
        contigs,
        combined_edges,
        name_prefix="raw_confirmed_master_contig",
    )
    lengths = [len(contig.sequence) for contig in projection]
    return projection, RawConfirmedMasterGraphDiagnostics(
        input_contigs=len(contigs),
        candidate_dovetails=candidate_dovetails,
        exact_dovetails=len(exact_edges) // 2,
        near_exact_candidates=near_candidates // 2,
        near_exact_accepted=len(preferred_near_edges) // 2,
        mismatch_positions=mismatch_positions // 2,
        insufficient_raw_support=insufficient_support // 2,
        strain_conflicts=strain_conflicts // 2,
        exact_preferred=exact_preferred,
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
    """Remove only contained copies whose alternatives lack raw-read support."""
    if not contigs:
        return [], StrainSafeContainmentDiagnostics()
    if not 0.0 <= minimum_identity <= 1.0:
        raise ValueError("containment identity must be between 0 and 1")
    if not 0.0 <= minimum_coverage <= 1.0:
        raise ValueError("containment coverage must be between 0 and 1")
    if minimum_primary_allele_support < 1:
        raise ValueError("primary allele support must be at least 1")
    if minimum_alternate_allele_support < 1:
        raise ValueError("alternate allele support must be at least 1")
    if minimum_support_margin < 1:
        raise ValueError("containment support margin must be at least 1")
    if minimum_base_quality < 0:
        raise ValueError("containment base quality must not be negative")
    if damage_end_window < 0:
        raise ValueError("containment damage end window must not be negative")

    molecules = list(range(len(raw_reads))) if molecule_ids is None else molecule_ids
    if len(molecules) != len(raw_reads):
        raise ValueError("molecule IDs must align one-to-one with raw reads")
    contig_position_bits = max(
        1, max(len(contig.sequence) for contig in contigs).bit_length()
    )
    contig_anchors = _anchor_index(
        contigs,
        anchor_k,
        0,
        anchors_per_read,
        maximum_anchor_occurrences,
    )
    raw_position_bits = max(
        1, max((len(read.sequence) for read in raw_reads), default=1).bit_length()
    )
    raw_anchors = (
        _anchor_index(
            raw_reads,
            anchor_k,
            0,
            anchors_per_read,
            maximum_anchor_occurrences,
        )
        if raw_reads
        else {}
    )
    contig_ids = list(range(len(contigs)))
    order = sorted(
        range(len(contigs)),
        key=lambda index: (
            -len(contigs[index].sequence),
            contigs[index].name,
            contigs[index].sequence,
        ),
    )
    unavailable = [True] * len(contigs)
    removed: set[int] = set()
    candidate_containments = exact_duplicates = near_candidates = 0
    strain_protected = insufficient_evidence = unsupported_variants = 0

    def raw_support(
        raw_candidates: list[_Alignment],
        position: int,
        alleles: frozenset[str],
    ) -> dict[str, int]:
        support: dict[str, set[int]] = defaultdict(set)
        damage_like = alleles in _DAMAGE_PAIRS
        for raw in raw_candidates:
            if not raw.offset <= position < raw.offset + len(raw.sequence):
                continue
            raw_position = position - raw.offset
            record = raw_reads[raw.read_index]
            quality_position = (
                raw_position
                if raw.sequence == record.sequence
                else len(record.sequence) - raw_position - 1
            )
            if record.qualities is None:
                continue
            if record.qualities[quality_position] < minimum_base_quality:
                continue
            if damage_like and (
                quality_position < damage_end_window
                or quality_position >= len(record.sequence) - damage_end_window
            ):
                continue
            base = raw.sequence[raw_position]
            if base in alleles:
                support[base].add(molecules[raw.read_index])
        return {base: len(indices) for base, indices in support.items()}

    for query_index in order:
        query = contigs[query_index].sequence
        candidates = _candidate_alignments(
            query,
            contigs,
            contig_ids,
            contig_anchors,
            unavailable,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_overlap=max(anchor_k, int(len(query) * minimum_coverage)),
            minimum_identity=minimum_identity,
            minimum_ry_identity=minimum_identity,
            position_bits=contig_position_bits,
            target_window=len(query),
        )
        containments: list[tuple[float, int, _Alignment]] = []
        for candidate in candidates:
            start = max(0, candidate.offset)
            stop = min(len(query), candidate.offset + len(candidate.sequence))
            coverage = (stop - start) / len(query)
            if coverage < minimum_coverage:
                continue
            overlap = candidate.sequence[
                start - candidate.offset : stop - candidate.offset
            ]
            identity = _identity(query[start:stop], overlap)
            containments.append((identity, stop - start, candidate))
        candidate_containments += len(containments)
        discard = False
        for _, _, candidate in sorted(
            containments,
            key=lambda item: (item[0], item[1], len(item[2].sequence)),
            reverse=True,
        ):
            start = max(0, candidate.offset)
            stop = min(len(query), candidate.offset + len(candidate.sequence))
            if start != 0 or stop != len(query):
                near_candidates += 1
                insufficient_evidence += 1
                continue
            container = candidate.sequence[
                -candidate.offset : len(query) - candidate.offset
            ]
            mismatches = [
                position
                for position, (first, second) in enumerate(zip(query, container))
                if first != second
            ]
            if not mismatches:
                exact_duplicates += 1
                discard = True
                break
            near_candidates += 1
            raw_candidates = _candidate_alignments(
                query,
                raw_reads,
                molecules,
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
                target_window=len(query),
            )
            supported_alternate = False
            primary_confirmed = True
            for position in mismatches:
                query_base = query[position]
                primary_base = container[position]
                support = raw_support(
                    raw_candidates,
                    position,
                    frozenset((query_base, primary_base)),
                )
                alternate_support = support.get(query_base, 0)
                primary_support = support.get(primary_base, 0)
                if alternate_support >= minimum_alternate_allele_support:
                    supported_alternate = True
                    break
                if (
                    primary_support < minimum_primary_allele_support
                    or primary_support - alternate_support < minimum_support_margin
                ):
                    primary_confirmed = False
            if supported_alternate:
                strain_protected += 1
                continue
            if not primary_confirmed:
                insufficient_evidence += 1
                continue
            unsupported_variants += 1
            discard = True
            break
        if discard:
            removed.add(query_index)
        else:
            unavailable[query_index] = False

    projection = [
        contig for index, contig in enumerate(contigs) if index not in removed
    ]
    lengths = [len(contig.sequence) for contig in projection]
    return projection, StrainSafeContainmentDiagnostics(
        input_contigs=len(contigs),
        candidate_containments=candidate_containments,
        exact_duplicates=exact_duplicates,
        near_duplicate_candidates=near_candidates,
        strain_protected=strain_protected,
        insufficient_evidence=insufficient_evidence,
        unsupported_variants=unsupported_variants,
        removed_contigs=len(removed),
        removed_bases=sum(len(contigs[index].sequence) for index in removed),
        projected_contigs=len(projection),
        projected_bases=sum(lengths),
        projected_n50=_n50(lengths),
        projected_longest_contig=max(lengths, default=0),
    )
