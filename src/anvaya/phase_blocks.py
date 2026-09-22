"""Truth-blind phase-link evidence from damage-aware molecule observations."""

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations

from anvaya.linked_allele_guard import variable_neighbours
from anvaya.raw_consensus import BASES, observation_likelihoods


@dataclass(frozen=True, slots=True)
class PhaseLinkEvidence:
    left_position_0based: int
    right_position_0based: int
    allele_pair_support: tuple[tuple[int, int, int], ...]
    allele_pair_molecules: tuple[tuple[int, int, tuple[int, ...]], ...]
    complete_allele_coverage: bool


@dataclass(frozen=True, slots=True)
class PhaseBlockEvidence:
    variable_positions_0based: tuple[int, ...]
    blocks: tuple[tuple[int, ...], ...]
    links: tuple[PhaseLinkEvidence, ...]
    ambiguous_links: tuple[PhaseLinkEvidence, ...]
    unlinked_pairs_0based: tuple[tuple[int, int], ...]
    marker_alleles: tuple[tuple[int, int, tuple[int, ...]], ...]
    excluded_alleles: tuple[tuple[int, int, tuple[int, ...]], ...]


@dataclass(frozen=True, slots=True)
class ComponentPlacementEvidence:
    read_index: int
    molecule_id: int
    offset: int
    reverse: bool
    read_start_0based: int
    read_stop_exclusive: int
    status: str


@dataclass(frozen=True, slots=True)
class ComponentPhaseDiagnostics:
    placements: int
    consistent_placements: int
    ambiguous_placements: int
    ambiguous_molecules: int
    usable_observations: int
    distinct_observed_molecules: int
    covered_sites: int
    uncovered_sites: int
    missing_quality_observations: int
    low_quality_observations: int
    ambiguous_base_observations: int
    outside_contig_observations: int
    duplicate_molecule_observations: int
    conflicting_molecule_observations: int


@dataclass(frozen=True, slots=True)
class ComponentPhaseEvidence:
    phase: PhaseBlockEvidence
    placements: tuple[ComponentPlacementEvidence, ...]
    diagnostics: ComponentPhaseDiagnostics


def _link_evidence(left, right, markers):
    shared = markers[left].keys() & markers[right].keys()
    pair_molecules = defaultdict(list)
    for molecule in sorted(shared):
        pair_molecules[(markers[left][molecule], markers[right][molecule])].append(molecule)
    supported = {pair: tuple(molecules) for pair, molecules in pair_molecules.items()
                 if len(molecules) >= 2}
    if not supported:
        return None

    left_alleles = set(markers[left].values())
    right_alleles = set(markers[right].values())
    complete = ({a for a, _ in supported} == left_alleles
                and {b for _, b in supported} == right_alleles)
    ordered = tuple((a, b, supported[a, b]) for a, b in sorted(supported))
    return PhaseLinkEvidence(
        left,
        right,
        tuple((a, b, len(molecules)) for a, b, molecules in ordered),
        ordered,
        complete,
    )


def _one_to_one(link):
    if not link.complete_allele_coverage:
        return False
    left_degree = Counter(a for a, _, _ in link.allele_pair_molecules)
    right_degree = Counter(b for _, b, _ in link.allele_pair_molecules)
    return set(left_degree.values()) == {1} and set(right_degree.values()) == {1}


def _split_consistent_links(links):
    parent = {}

    def root(node):
        parent.setdefault(node, node)
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left, right):
        left_root, right_root = root(left), root(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for link in links:
        for left_allele, right_allele, _ in link.allele_pair_molecules:
            union((link.left_position_0based, left_allele),
                  (link.right_position_0based, right_allele))

    alleles_by_root = defaultdict(lambda: defaultdict(set))
    for position, allele in parent:
        alleles_by_root[root((position, allele))][position].add(allele)
    conflicting_roots = {
        component for component, positions in alleles_by_root.items()
        if any(len(alleles) > 1 for alleles in positions.values())
    }

    consistent, conflicting = [], []
    for link in links:
        roots = {
            root((position, allele))
            for left_allele, right_allele, _ in link.allele_pair_molecules
            for position, allele in (
                (link.left_position_0based, left_allele),
                (link.right_position_0based, right_allele),
            )
        }
        (conflicting if roots & conflicting_roots else consistent).append(link)
    return consistent, conflicting


def _connected_blocks(positions, links):
    parent = {position: position for position in positions}

    def root(position):
        while parent[position] != position:
            parent[position] = parent[parent[position]]
            position = parent[position]
        return position

    for link in links:
        left_root = root(link.left_position_0based)
        right_root = root(link.right_position_0based)
        if left_root != right_root:
            parent[right_root] = left_root
    groups = defaultdict(list)
    for position in positions:
        groups[root(position)].append(position)
    return tuple(sorted((tuple(group) for group in groups.values()), key=lambda group: group[0]))


def phase_blocks(columns):
    """Summarize same-molecule allele co-occurrence at robust variable sites.

    Positions are 0-based. Marker alleles require at least two robust molecules;
    lower-support alleles are retained in ``excluded_alleles`` for auditing but
    cannot veto a link. Every pair of variable sites is examined. ``links`` holds
    only complete one-to-one allele matchings; branching or incomplete supported
    pair evidence is retained in ``ambiguous_links`` and does not join blocks.
    Blocks express connected linkage evidence, not complete haplotypes.
    Input columns are already deduplicated by molecule; conflicting observations
    are represented by ``None`` and excluded by ``variable_neighbours``.
    """
    variable = variable_neighbours(columns)
    markers = {}
    excluded = []
    for position, observations in variable.items():
        counts = Counter(observations.values())
        eligible = {allele for allele, count in counts.items() if count >= 2}
        markers[position] = {molecule: allele for molecule, allele in observations.items()
                             if allele in eligible}
        for allele in sorted(set(counts) - eligible):
            excluded.append((position, allele, tuple(sorted(
                molecule for molecule, observed in observations.items() if observed == allele))))

    positions = tuple(sorted(markers))
    marker_alleles = tuple(
        (position, allele, tuple(sorted(
            molecule for molecule, observed in markers[position].items() if observed == allele)))
        for position in positions
        for allele in sorted(set(markers[position].values()))
    )
    candidate_links, ambiguous, unlinked = [], [], []
    for left, right in combinations(positions, 2):
        evidence = _link_evidence(left, right, markers)
        if evidence is None:
            unlinked.append((left, right))
        elif _one_to_one(evidence):
            candidate_links.append(evidence)
        else:
            ambiguous.append(evidence)
    links, conflicting = _split_consistent_links(candidate_links)
    ambiguous.extend(conflicting)
    ambiguous.sort(key=lambda link: (link.left_position_0based, link.right_position_0based))
    return PhaseBlockEvidence(
        positions,
        _connected_blocks(positions, links),
        tuple(links),
        tuple(ambiguous),
        tuple(unlinked),
        marker_alleles,
        tuple(excluded),
    )


def phase_placements(pool, target_length, raw_placements, profile):
    """Project 0-based raw placements into truth-blind phase evidence.

    Raw intervals and target positions are 0-based and half-open. The adapter
    consumes only immutable reads, molecule IDs, placements and the supplied
    damage profile. Multiple physical placements of one read make its entire
    molecule ambiguous; no reference or simulation truth is accepted.
    """
    if target_length < 1:
        raise ValueError("phase target length must be at least 1")
    raw_placements = tuple(raw_placements)
    mappings = defaultdict(set)
    for placement in raw_placements:
        if not 0 <= placement.read_index < len(pool.records):
            raise ValueError("raw placement read index is outside the pool")
        mappings[placement.read_index].add((placement.offset, placement.reverse))
    ambiguous_molecules = {
        pool.records[read_index].molecule_id
        for read_index, physical in mappings.items() if len(physical) > 1
    }

    placement_evidence = []
    columns = [{} for _ in range(target_length)]
    missing_quality = low_quality = ambiguous_bases = outside = 0
    duplicate = conflicting = usable = 0
    for placement in raw_placements:
        source = pool.records[placement.read_index]
        read = source.raw
        read_length = len(read.sequence)
        if not 0 <= placement.read_start < placement.read_stop <= read_length:
            raise ValueError("raw placement interval is outside the original read")
        ambiguous = source.molecule_id in ambiguous_molecules
        placement_evidence.append(ComponentPlacementEvidence(
            placement.read_index,
            source.molecule_id,
            placement.offset,
            placement.reverse,
            placement.read_start,
            placement.read_stop,
            "ambiguous" if ambiguous else "consistent",
        ))
        if ambiguous:
            continue
        for raw_position in range(placement.read_start, placement.read_stop):
            position = placement.offset + (
                read_length - 1 - raw_position if placement.reverse else raw_position
            )
            if not 0 <= position < target_length:
                outside += 1
                continue
            if read.qualities is None:
                missing_quality += 1
                continue
            quality = read.qualities[raw_position]
            if quality < 3:
                low_quality += 1
                continue
            raw_base = read.sequence[raw_position]
            if raw_base not in BASES:
                ambiguous_bases += 1
                continue
            likelihoods = observation_likelihoods(
                raw_base, quality, raw_position, read_length, placement.reverse, profile,
            )
            allele = 3 - BASES.index(raw_base) if placement.reverse else BASES.index(raw_base)
            observation = (allele, quality, likelihoods)
            column = columns[position]
            if source.molecule_id in column:
                duplicate += 1
                previous = column[source.molecule_id]
                if previous is None:
                    continue
                if previous[0] != allele:
                    column[source.molecule_id] = None
                    conflicting += 1
                    continue
                rank = lambda item: (
                    item[1], -max(p for base, p in enumerate(item[2]) if base != item[0])
                )
                if rank(previous) >= rank(observation):
                    continue
            column[source.molecule_id] = observation
            usable += 1

    observed_molecules = {
        molecule for column in columns for molecule, observation in column.items()
        if observation is not None
    }
    covered_sites = sum(any(value is not None for value in column.values()) for column in columns)
    diagnostics = ComponentPhaseDiagnostics(
        len(raw_placements),
        sum(item.status == "consistent" for item in placement_evidence),
        sum(item.status == "ambiguous" for item in placement_evidence),
        len(ambiguous_molecules),
        usable,
        len(observed_molecules),
        covered_sites,
        target_length - covered_sites,
        missing_quality,
        low_quality,
        ambiguous_bases,
        outside,
        duplicate,
        conflicting,
    )
    return ComponentPhaseEvidence(phase_blocks(columns), tuple(placement_evidence), diagnostics)


def phase_component(pool, record, profile):
    """Project one derived component's raw placements into phase evidence."""
    if record.raw_placements is None:
        raise ValueError("component raw placements are required")
    return phase_placements(
        pool, len(record.current.sequence), record.raw_placements, profile,
    )
