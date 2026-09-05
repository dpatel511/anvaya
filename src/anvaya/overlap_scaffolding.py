"""Conservative paired-end scaffolding of progressive overlap contigs."""

from collections import defaultdict
from dataclasses import dataclass
from statistics import median

from anvaya.damage_consensus import _anchor_index, _identity
from anvaya.overlap_assembly import _candidate_alignments, _n50, _ry
from anvaya.reads import Read
from anvaya.sequences import canonical_sequence, reverse_complement


@dataclass(frozen=True, slots=True)
class _Placement:
    contig: int
    reverse: bool
    start: int
    end: int


@dataclass(frozen=True, slots=True)
class _ScaffoldEdge:
    first: tuple[int, str]
    second: tuple[int, str]
    support: int
    gap: int


@dataclass(frozen=True, slots=True)
class PairedScaffoldDiagnostics:
    """Counts describing a report-only paired progressive scaffold graph."""

    input_contigs: int = 0
    paired_molecules: int = 0
    uniquely_mapped_pairs: int = 0
    ambiguous_reads: int = 0
    same_contig_pairs: int = 0
    calibration_pairs: int = 0
    calibration_valid: bool = False
    fragment_mean: float = 0.0
    fragment_sd: float = 0.0
    cross_contig_pairs: int = 0
    negative_gap_pairs: int = 0
    candidate_links: int = 0
    supported_links: int = 0
    ambiguous_ends: int = 0
    reciprocal_links: int = 0
    cyclic_components: int = 0
    linear_paths: int = 0
    joined_contigs: int = 0
    inserted_gap_bases: int = 0
    projected_contigs: int = 0
    projected_bases: int = 0
    projected_n50: int = 0
    projected_longest_contig: int = 0


def _median_absolute_deviation(values: list[int], center: float) -> float:
    return float(median([abs(value - center) for value in values]))


def _unique_placement(
    read: Read,
    contigs: list[Read],
    contig_ids: list[int],
    anchors: dict[int, list[int]],
    *,
    anchor_k: int,
    anchors_per_read: int,
    maximum_anchor_occurrences: int,
    minimum_anchor_matches: int,
    minimum_identity: float,
    minimum_ry_identity: float,
    position_bits: int,
    target_window: int,
) -> tuple[_Placement | None, bool]:
    if len(read.sequence) < anchor_k:
        return None, False
    candidates = _candidate_alignments(
        read.sequence,
        contigs,
        contig_ids,
        anchors,
        set(),
        anchor_k=anchor_k,
        anchors_per_read=anchors_per_read,
        maximum_anchor_occurrences=maximum_anchor_occurrences,
        minimum_anchor_matches=minimum_anchor_matches,
        minimum_overlap=len(read.sequence),
        minimum_identity=minimum_identity,
        minimum_ry_identity=minimum_ry_identity,
        position_bits=position_bits,
        target_window=target_window,
    )
    ranked: list[tuple[tuple[float, float, int], _Placement]] = []
    for candidate in candidates:
        oriented_start = -candidate.offset
        if oriented_start < 0:
            continue
        oriented_end = oriented_start + len(read.sequence)
        if oriented_end > len(candidate.sequence):
            continue
        segment = candidate.sequence[oriented_start:oriented_end]
        dna_identity = _identity(read.sequence, segment)
        ry_identity = _identity(_ry(read.sequence), _ry(segment))
        original = contigs[candidate.read_index].sequence
        reversed_original = reverse_complement(original)
        if original == reversed_original:
            continue
        reverse = candidate.sequence == reversed_original
        start = (
            len(original) - oriented_end if reverse else oriented_start
        )
        placement = _Placement(
            candidate.read_index,
            reverse,
            start,
            start + len(read.sequence),
        )
        ranked.append(((dna_identity, ry_identity, candidate.support), placement))
    if not ranked:
        return None, False
    ranked.sort(key=lambda item: item[0], reverse=True)
    best_score = ranked[0][0]
    best = {item[1] for item in ranked if item[0] == best_score}
    if len(best) != 1:
        return None, True
    return next(iter(best)), False


def _calibration_span(first: _Placement, second: _Placement) -> int | None:
    if first.contig != second.contig or first.reverse == second.reverse:
        return None
    forward = second if first.reverse else first
    reverse = first if first.reverse else second
    if forward.start >= reverse.start:
        return None
    return reverse.end - forward.start


def _connected_end(placement: _Placement) -> str:
    return "left" if placement.reverse else "right"


def _distance_to_connected_end(placement: _Placement, contig_length: int) -> int:
    return placement.end if placement.reverse else contig_length - placement.start


def _canonical_link(
    first: tuple[int, str], second: tuple[int, str]
) -> tuple[tuple[int, str], tuple[int, str]]:
    return (first, second) if first < second else (second, first)


def _oriented_nodes(edge: _ScaffoldEdge) -> tuple[tuple[int, bool], tuple[int, bool]]:
    source = (edge.first[0], edge.first[1] == "left")
    target = (edge.second[0], edge.second[1] == "right")
    return source, target


def scaffold_progressive_contigs(
    contigs: list[Read],
    left_reads: list[Read],
    right_reads: list[Read],
    *,
    anchor_k: int = 15,
    anchors_per_read: int = 8,
    maximum_anchor_occurrences: int = 100,
    minimum_anchor_matches: int = 2,
    minimum_mapping_identity: float = 0.95,
    minimum_mapping_ry_identity: float = 0.99,
    minimum_link_support: int = 3,
    dominance_ratio: float = 3.0,
    fragment_mean: float | None = None,
    fragment_sd: float | None = None,
    minimum_calibration_pairs: int = 20,
) -> tuple[list[Read], PairedScaffoldDiagnostics]:
    """Project reciprocal paired links without changing input contigs."""
    if len(left_reads) != len(right_reads):
        raise ValueError("paired-end inputs must contain the same number of reads")
    if minimum_link_support < 1:
        raise ValueError("minimum_link_support must be at least 1")
    if dominance_ratio < 1.0:
        raise ValueError("dominance_ratio must be at least 1")
    if minimum_calibration_pairs < 1:
        raise ValueError("minimum_calibration_pairs must be at least 1")
    if (fragment_mean is None) != (fragment_sd is None):
        raise ValueError("fragment_mean and fragment_sd must be provided together")
    if fragment_mean is not None and fragment_mean <= 0:
        raise ValueError("fragment_mean must be positive")
    if fragment_sd is not None and fragment_sd < 0:
        raise ValueError("fragment_sd must not be negative")
    if not contigs or not left_reads:
        return list(contigs), PairedScaffoldDiagnostics(
            input_contigs=len(contigs),
            paired_molecules=len(left_reads),
            projected_contigs=len(contigs),
            projected_bases=sum(len(contig.sequence) for contig in contigs),
            projected_n50=_n50([len(contig.sequence) for contig in contigs]),
            projected_longest_contig=max(
                (len(contig.sequence) for contig in contigs), default=0
            ),
        )

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
    placements: list[tuple[_Placement, _Placement]] = []
    ambiguous_reads = 0
    for left, right in zip(left_reads, right_reads):
        first, first_ambiguous = _unique_placement(
            left,
            contigs,
            contig_ids,
            anchors,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_identity=minimum_mapping_identity,
            minimum_ry_identity=minimum_mapping_ry_identity,
            position_bits=position_bits,
            target_window=target_window,
        )
        second, second_ambiguous = _unique_placement(
            right,
            contigs,
            contig_ids,
            anchors,
            anchor_k=anchor_k,
            anchors_per_read=anchors_per_read,
            maximum_anchor_occurrences=maximum_anchor_occurrences,
            minimum_anchor_matches=minimum_anchor_matches,
            minimum_identity=minimum_mapping_identity,
            minimum_ry_identity=minimum_mapping_ry_identity,
            position_bits=position_bits,
            target_window=target_window,
        )
        ambiguous_reads += int(first_ambiguous) + int(second_ambiguous)
        if first is not None and second is not None:
            placements.append((first, second))

    calibration = [
        span
        for first, second in placements
        if (span := _calibration_span(first, second)) is not None
    ]
    if fragment_mean is None:
        if len(calibration) < minimum_calibration_pairs:
            lengths = [len(contig.sequence) for contig in contigs]
            return list(contigs), PairedScaffoldDiagnostics(
                input_contigs=len(contigs),
                paired_molecules=len(left_reads),
                uniquely_mapped_pairs=len(placements),
                ambiguous_reads=ambiguous_reads,
                same_contig_pairs=sum(
                    first.contig == second.contig for first, second in placements
                ),
                calibration_pairs=len(calibration),
                projected_contigs=len(contigs),
                projected_bases=sum(lengths),
                projected_n50=_n50(lengths),
                projected_longest_contig=max(lengths, default=0),
            )
        estimated_mean = float(median(calibration))
        typical_read_length = float(
            median(
                [
                    len(read.sequence)
                    for pair in zip(left_reads, right_reads)
                    for read in pair
                ]
            )
        )
        if estimated_mean < typical_read_length:
            lengths = [len(contig.sequence) for contig in contigs]
            return list(contigs), PairedScaffoldDiagnostics(
                input_contigs=len(contigs),
                paired_molecules=len(left_reads),
                uniquely_mapped_pairs=len(placements),
                ambiguous_reads=ambiguous_reads,
                same_contig_pairs=sum(
                    first.contig == second.contig for first, second in placements
                ),
                calibration_pairs=len(calibration),
                fragment_mean=estimated_mean,
                projected_contigs=len(contigs),
                projected_bases=sum(lengths),
                projected_n50=_n50(lengths),
                projected_longest_contig=max(lengths, default=0),
            )
        fragment_mean = estimated_mean
        fragment_sd = max(
            1.0,
            1.4826 * _median_absolute_deviation(calibration, fragment_mean),
        )
    assert fragment_sd is not None

    estimates: dict[tuple[tuple[int, str], tuple[int, str]], list[int]] = (
        defaultdict(list)
    )
    cross_contig_pairs = negative_gap_pairs = 0
    for first, second in placements:
        if first.contig == second.contig:
            continue
        cross_contig_pairs += 1
        first_end = (first.contig, _connected_end(first))
        second_end = (second.contig, _connected_end(second))
        gap = round(
            fragment_mean
            - _distance_to_connected_end(
                first, len(contigs[first.contig].sequence)
            )
            - _distance_to_connected_end(
                second, len(contigs[second.contig].sequence)
            )
        )
        if gap < 0:
            negative_gap_pairs += 1
            continue
        estimates[_canonical_link(first_end, second_end)].append(gap)

    edges: list[_ScaffoldEdge] = []
    for (first, second), gaps in estimates.items():
        center = float(median(gaps))
        tolerance = max(10.0, 3.0 * fragment_sd)
        consistent = [gap for gap in gaps if abs(gap - center) <= tolerance]
        if len(consistent) >= minimum_link_support:
            edges.append(
                _ScaffoldEdge(
                    first,
                    second,
                    len(consistent),
                    round(median(consistent)),
                )
            )

    by_end: dict[tuple[int, str], list[_ScaffoldEdge]] = defaultdict(list)
    for edge in edges:
        by_end[edge.first].append(edge)
        by_end[edge.second].append(edge)
    winners: dict[tuple[int, str], _ScaffoldEdge] = {}
    ambiguous_ends = 0
    for end, candidates in by_end.items():
        ranked = sorted(candidates, key=lambda edge: (-edge.support, edge.first, edge.second))
        winner = ranked[0]
        runner_up = ranked[1].support if len(ranked) > 1 else 0
        if runner_up and winner.support < dominance_ratio * runner_up:
            ambiguous_ends += 1
            continue
        winners[end] = winner
    accepted = [
        edge
        for edge in edges
        if winners.get(edge.first) == edge and winners.get(edge.second) == edge
    ]

    outgoing: dict[tuple[int, bool], tuple[tuple[int, bool], _ScaffoldEdge]] = {}
    incoming: dict[tuple[int, bool], tuple[tuple[int, bool], _ScaffoldEdge]] = {}
    for edge in accepted:
        source, target = _oriented_nodes(edge)
        outgoing[source] = (target, edge)
        incoming[target] = (source, edge)
        reverse_source = (target[0], not target[1])
        reverse_target = (source[0], not source[1])
        outgoing[reverse_source] = (reverse_target, edge)
        incoming[reverse_target] = (reverse_source, edge)

    visited_contigs: set[int] = set()
    projected: list[Read] = []
    linear_paths = inserted_gap_bases = 0
    starts = sorted(node for node in outgoing if node not in incoming)
    for start in starts:
        if start[0] in visited_contigs:
            continue
        path = [start]
        edges_in_path: list[_ScaffoldEdge] = []
        current = start
        while current in outgoing:
            target, edge = outgoing[current]
            if target[0] in {node[0] for node in path}:
                break
            edges_in_path.append(edge)
            path.append(target)
            current = target
        if len(path) < 2:
            continue
        sequence = (
            reverse_complement(contigs[path[0][0]].sequence)
            if path[0][1]
            else contigs[path[0][0]].sequence
        )
        for node, edge in zip(path[1:], edges_in_path):
            gap = max(1, edge.gap)
            sequence += "N" * gap
            sequence += (
                reverse_complement(contigs[node[0]].sequence)
                if node[1]
                else contigs[node[0]].sequence
            )
            inserted_gap_bases += gap
        canonical = canonical_sequence(sequence)
        projected.append(Read(f"paired_scaffold_{linear_paths + 1}", canonical))
        visited_contigs.update(node[0] for node in path)
        linear_paths += 1

    cyclic_components: set[frozenset[int]] = set()
    for source in outgoing:
        if source[0] in visited_contigs:
            continue
        component: set[int] = set()
        current = source
        while current in outgoing and current[0] not in component:
            component.add(current[0])
            current = outgoing[current][0]
        if current[0] in component:
            cyclic_components.add(frozenset(component))

    projected.extend(
        contig for index, contig in enumerate(contigs) if index not in visited_contigs
    )
    lengths = [len(contig.sequence) for contig in projected]
    return projected, PairedScaffoldDiagnostics(
        input_contigs=len(contigs),
        paired_molecules=len(left_reads),
        uniquely_mapped_pairs=len(placements),
        ambiguous_reads=ambiguous_reads,
        same_contig_pairs=sum(
            first.contig == second.contig for first, second in placements
        ),
        calibration_pairs=len(calibration),
        calibration_valid=True,
        fragment_mean=fragment_mean,
        fragment_sd=fragment_sd,
        cross_contig_pairs=cross_contig_pairs,
        negative_gap_pairs=negative_gap_pairs,
        candidate_links=len(estimates),
        supported_links=len(edges),
        ambiguous_ends=ambiguous_ends,
        reciprocal_links=len(accepted),
        cyclic_components=len(cyclic_components),
        linear_paths=linear_paths,
        joined_contigs=len(visited_contigs),
        inserted_gap_bases=inserted_gap_bases,
        projected_contigs=len(projected),
        projected_bases=sum(lengths),
        projected_n50=_n50(lengths),
        projected_longest_contig=max(lengths, default=0),
    )
