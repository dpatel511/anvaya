"""Bounded exact and damage-aware read-overlap string graph."""
from collections import Counter, defaultdict
from dataclasses import asdict, replace
from time import perf_counter

from anvaya.overlap_assembly import (
    _candidate_alignments,
    _CandidateDiagnostics,
    _MasterOverlapEdge,
)
from anvaya.overlap_graph import _project_master_edges
from anvaya.overlap_index import _anchor_index
from anvaya.overlap_progressive import ProgressiveSequencePool
from anvaya.overlap_progressive_links import _add_bidirected_edge
from anvaya.phase_graph import phase_graph_components
from anvaya.raw_consensus import DamageProfile, RawPlacement
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


def _complement(base):
    return reverse_complement(base)


def _physical_edge_signature(edge, sequences):
    """Identify one bidirected edge, including its spelled corrections."""
    forward = (
        edge.source, edge.target, edge.shift, edge.overlap, edge.corrections,
    )
    target_length = len(sequences[edge.target[0]].sequence)
    reverse = (
        (edge.target[0], not edge.target[1]),
        (edge.source[0], not edge.source[1]),
        target_length - edge.overlap,
        edge.overlap,
        tuple(sorted(
            (target_length - 1 - (position - edge.shift), _complement(base))
            for position, base in edge.corrections
        )),
    )
    return min(forward, reverse)


def _exact_containment_candidates(sequences, maximum_seed_length=31):
    """Return every exact shorter-in-longer placement using indexed prefixes.

    Placement offsets are 0-based. Seeds discover candidates only; the complete
    oriented child is verified before a placement is returned. Grouping by seed
    length keeps reads shorter than ``maximum_seed_length`` complete as well.
    """
    seeds_by_length = defaultdict(lambda: defaultdict(list))
    for child_index, child in enumerate(sequences):
        orientations = [(False, child)]
        reverse_child = reverse_complement(child)
        if reverse_child != child:
            orientations.append((True, reverse_child))
        seed_length = min(maximum_seed_length, len(child))
        for reverse, oriented_child in orientations:
            seeds_by_length[seed_length][oriented_child[:seed_length]].append(
                (child_index, reverse, oriented_child)
            )

    candidates = defaultdict(set)
    for parent_index, parent in enumerate(sequences):
        for seed_length, seed_index in seeds_by_length.items():
            if seed_length >= len(parent):
                continue
            for offset in range(len(parent) - seed_length + 1):
                for child_index, reverse, oriented_child in seed_index.get(
                    parent[offset:offset + seed_length], ()
                ):
                    if len(oriented_child) >= len(parent):
                        continue
                    if parent.startswith(oriented_child, offset):
                        candidates[child_index].add(
                            (parent_index, reverse, offset)
                        )
    return {
        child_index: candidates[child_index]
        for child_index in sorted(candidates)
    }


def _damage_compatible_edge(
    reads, molecule_ids, source_sequence, target_sequence, offset, profile,
    raw_candidates,
):
    """Accept only mismatches explained by supplied terminal damage.

    All read and overlap positions are 0-based. ``offset`` places the target's
    first oriented base relative to the source's first oriented base.
    """
    start = max(0, offset)
    stop = min(len(source_sequence), offset + len(target_sequence))
    mismatches = [
        position for position in range(start, stop)
        if source_sequence[position] != target_sequence[position - offset]
    ]
    if not mismatches:
        return (), "exact"

    corrections = []
    for source_position in mismatches:
        observations = defaultdict(set)
        missing_quality = low_quality = False
        for raw in raw_candidates:
            if not raw.offset <= source_position < raw.offset + len(raw.sequence):
                continue
            oriented_position = source_position - raw.offset
            raw_index = raw.read_index
            effective_reverse = raw.reverse
            raw_position = (
                len(reads[raw_index].sequence) - 1 - oriented_position
                if effective_reverse else oriented_position
            )
            qualities = reads[raw_index].qualities
            if qualities is None:
                missing_quality = True
                continue
            if qualities[raw_position] < 20:
                low_quality = True
                continue
            raw_base = reads[raw_index].sequence[raw_position]
            observations[molecule_ids[raw_index]].add(
                (raw_base, raw_position, len(reads[raw_index].sequence), effective_reverse)
            )

        unambiguous = [next(iter(values)) for values in observations.values() if len(values) == 1]
        # The two edge endpoints propose the link; a third molecule supplies
        # independent allele evidence rather than letting a pair explain itself.
        if len(unambiguous) < 3:
            if not unambiguous and missing_quality:
                return (), "missing_quality"
            if not unambiguous and low_quality:
                return (), "low_quality"
            return (), "insufficient_molecules"
        observed_bases = {source_sequence[source_position], target_sequence[source_position - offset]}
        compatible_latents = []
        for latent in observed_bases:
            exact = damaged = 0
            for raw_base, raw_position, read_length, reverse in unambiguous:
                raw_latent = _complement(latent) if reverse else latent
                if raw_base == raw_latent:
                    exact += 1
                    continue
                damage_rate = 0.0
                if raw_latent == "C" and raw_base == "T" and raw_position < len(profile.five_prime_ct):
                    damage_rate = profile.five_prime_ct[raw_position]
                elif raw_latent == "G" and raw_base == "A":
                    distance_3p = read_length - 1 - raw_position
                    if distance_3p < len(profile.three_prime_ga):
                        damage_rate = profile.three_prime_ga[distance_3p]
                if damage_rate > 0:
                    damaged += 1
                    continue
                break
            else:
                if exact and damaged:
                    compatible_latents.append(latent)
        if len(compatible_latents) != 1:
            return (), "non_directional_or_internal"
        corrections.append((source_position, compatible_latents[0]))
    return tuple(corrections), "damage_compatible"


def assemble(
    reads, molecule_ids=None, damage_profile=None, maximum_reads=1000,
    stage_timings=None, reduced_graph_audit=None,
):
    if stage_timings is not None:
        stage_timings.clear()
    if not reads:
        raise ValueError("string graph requires at least one read")
    if maximum_reads < 1:
        raise ValueError("maximum reads must be at least 1")
    if len(reads) > maximum_reads:
        raise ValueError(f"bounded string graph requires <={maximum_reads} reads")
    if molecule_ids is None:
        molecule_ids = list(range(len(reads)))
    if len(molecule_ids) != len(reads):
        raise ValueError("molecule IDs must align one-to-one with reads")
    groups = {}
    for index, read in enumerate(reads):
        canonical = min(read.sequence, reverse_complement(read.sequence))
        groups.setdefault(canonical, []).append((index, read.sequence != canonical))
    sequences = list(groups)
    nodes = [Read(f"node_{i}", sequence) for i, sequence in enumerate(sequences)]
    started = perf_counter()
    anchors = _anchor_index(nodes, 15, 0, 8, 100)
    raw_anchors = _anchor_index(reads, 15, 0, 8, 100)
    if stage_timings is not None:
        stage_timings["indexing"] = perf_counter() - started
    maximum_length = max(len(read.sequence) for read in reads)
    position_bits = maximum_length.bit_length()
    edges = {}
    ambiguous_edges = set()
    containment_candidates = defaultdict(set)
    edge_counts = Counter()
    candidate_search = _CandidateDiagnostics()
    raw_candidate_search = _CandidateDiagnostics()
    started = perf_counter()
    for index, node in enumerate(nodes):
        raw_candidates = _candidate_alignments(
            node.sequence, reads, list(range(len(reads))), raw_anchors, set(),
            anchor_k=15, anchors_per_read=8, maximum_anchor_occurrences=100,
            minimum_anchor_matches=1, minimum_overlap=30, minimum_identity=.9,
            minimum_ry_identity=.99, position_bits=position_bits,
            target_window=maximum_length, diagnostics=raw_candidate_search,
        ) if damage_profile is not None else ()
        candidates = _candidate_alignments(node.sequence, nodes, list(range(len(nodes))), anchors, {index},
            anchor_k=15, anchors_per_read=8, maximum_anchor_occurrences=100,
            minimum_anchor_matches=1, minimum_overlap=30,
            minimum_identity=.9 if damage_profile is not None else 1.,
            minimum_ry_identity=.99 if damage_profile is not None else 1.,
            position_bits=position_bits, target_window=maximum_length,
            diagnostics=candidate_search)
        for candidate in candidates:
            candidate_node = (candidate.read_index, candidate.reverse)
            corrections, reason = _damage_compatible_edge(
                reads, molecule_ids, node.sequence, candidate.sequence,
                candidate.offset, damage_profile, raw_candidates,
            ) if damage_profile is not None else ((), "exact")
            edge_counts[reason] += 1
            if reason not in ("exact", "damage_compatible"):
                continue
            source_length = len(node.sequence)
            candidate_length = len(candidate.sequence)
            left_extension = max(0, -candidate.offset)
            right_extension = max(
                0, candidate.offset + candidate_length - source_length
            )
            if right_extension and not left_extension:
                edge = _MasterOverlapEdge((index, False), candidate_node,
                                          candidate.offset,
                                          source_length - candidate.offset,
                                          corrections)
                _add_bidirected_edge(edges, nodes, edge, ambiguous_edges)
            elif left_extension and not right_extension:
                reverse_corrections = tuple(
                    (position - candidate.offset, base) for position, base in corrections
                )
                edge = _MasterOverlapEdge(candidate_node, (index, False),
                                          -candidate.offset,
                                          candidate_length + candidate.offset,
                                          reverse_corrections)
                _add_bidirected_edge(edges, nodes, edge, ambiguous_edges)
    if stage_timings is not None:
        stage_timings["candidate_edges"] = perf_counter() - started
    # Enumerate exact 0-based placements independently of anchor tie-breaking.
    # A contained node is safe to retire only when it has one physical placement.
    started = perf_counter()
    containment_candidates = _exact_containment_candidates(sequences)
    containments = {}
    for child, candidates in containment_candidates.items():
        if len(candidates) == 1:
            containments[child] = next(iter(candidates))
    contained_nodes = set(containments)
    repeated_parent_containment_nodes = {
        child for child, candidates in containment_candidates.items()
        if any(count > 1 for count in Counter(
            candidate[0] for candidate in candidates
        ).values())
    }
    if stage_timings is not None:
        stage_timings["containment"] = perf_counter() - started

    # Audit phase evidence before containment retirement, transitive reduction,
    # reciprocal-edge selection or path spelling can discard alternatives.
    started = perf_counter()
    phase_edges = dict(edges)
    for child, (parent, reverse, offset) in containments.items():
        forward = _MasterOverlapEdge(
            (parent, False), (child, reverse), offset, len(sequences[child]), (),
        )
        reverse_offset = len(sequences[parent]) - offset - len(sequences[child])
        mirrored = _MasterOverlapEdge(
            (parent, True), (child, not reverse), reverse_offset,
            len(sequences[child]), (),
        )
        phase_edges[(forward.source, forward.target)] = forward
        phase_edges[(mirrored.source, mirrored.target)] = mirrored
    prelayout_phase = phase_graph_components(
        reads,
        molecule_ids,
        sequences,
        groups,
        phase_edges,
        damage_profile if damage_profile is not None else DamageProfile((), ()),
    )
    consistent_phase = [component for component in prelayout_phase
                        if component.phase is not None]
    if stage_timings is not None:
        stage_timings["phase_audit"] = perf_counter() - started

    def ambiguous_full_containment(edge):
        return any(
            node in repeated_parent_containment_nodes
            and edge.overlap >= len(sequences[node])
            for node in (edge.source[0], edge.target[0])
        )

    started = perf_counter()
    edges = {
        key: edge for key, edge in edges.items()
        if edge.source[0] not in contained_nodes
        and edge.target[0] not in contained_nodes
        and not ambiguous_full_containment(edge)
    }
    layouts = {}
    contigs, diagnostics = _project_master_edges(
        nodes, edges, name_prefix="path", path_layouts=layouts,
        verify_transitive_alleles=damage_profile is not None,
        excluded_contigs=contained_nodes,
        reduced_graph_audit=reduced_graph_audit,
    )
    children = defaultdict(list)
    for child, (parent, reverse, offset) in containments.items():
        children[parent].append((child, reverse, offset))

    def add_contained(layout, node, reverse, offset):
        layout.append((node, reverse, offset))
        for child, child_reverse, child_offset in children.get(node, ()):
            if reverse:
                placed_reverse = not child_reverse
                placed_offset = (
                    offset + len(sequences[node]) - child_offset - len(sequences[child])
                )
            else:
                placed_reverse = child_reverse
                placed_offset = offset + child_offset
            add_contained(layout, child, placed_reverse, placed_offset)

    for name, path_layout in tuple(layouts.items()):
        expanded = []
        for node, reverse, offset in path_layout:
            add_contained(expanded, node, reverse, offset)
        layouts[name] = tuple(expanded)
    if stage_timings is not None:
        stage_timings["projection"] = perf_counter() - started
    started = perf_counter()
    original_pool = ProgressiveSequencePool.from_reads(reads, molecule_ids)
    records = [record.consumed() for record in original_pool.records]
    # Representatives are storage slots only; all path evidence is explicit.
    for slot, contig in enumerate(contigs):
        placements = []
        for node, reverse, offset in layouts[contig.name]:
            for raw_index, raw_reverse in groups[sequences[node]]:
                raw_length = len(reads[raw_index].sequence)
                placements.append(
                    RawPlacement(raw_index, offset, reverse != raw_reverse, 0, raw_length)
                )
                if sequences[node] == reverse_complement(sequences[node]):
                    placements.append(
                        RawPlacement(raw_index, offset, reverse == raw_reverse, 0, raw_length)
                    )
        record = replace(original_pool.records[slot].corrected(contig),
                         raw_placements=tuple(placements),
                         contributing_molecules=frozenset(
                             original_pool.records[p.read_index].molecule_id
                             for p in placements
                         ))
        records[slot] = record
    pool = ProgressiveSequencePool(tuple(records))
    if stage_timings is not None:
        stage_timings["provenance"] = perf_counter() - started
    return pool, dict(
        asdict(diagnostics),
        nodes=len(nodes),
        directed_edges=len(edges),
        ambiguous_physical_edges=len(ambiguous_edges),
        contained_nodes=len(contained_nodes),
        ambiguous_containments=len(containment_candidates) - len(contained_nodes),
        candidate_search=asdict(candidate_search),
        raw_candidate_search=asdict(raw_candidate_search),
        candidate_classifications=dict(edge_counts),
        prelayout_phase_components=len(prelayout_phase),
        prelayout_phase_coordinate_consistent=len(consistent_phase),
        prelayout_phase_coordinate_unresolved=len(prelayout_phase) - len(consistent_phase),
        prelayout_phase_variable_sites=sum(
            len(component.phase.phase.variable_positions_0based)
            for component in consistent_phase
        ),
        prelayout_phase_resolved_links=sum(
            len(component.phase.phase.links) for component in consistent_phase
        ),
        prelayout_phase_ambiguous_links=sum(
            len(component.phase.phase.ambiguous_links) for component in consistent_phase
        ),
        prelayout_phase_unlinked_pairs=sum(
            len(component.phase.phase.unlinked_pairs_0based)
            for component in consistent_phase
        ),
    )
