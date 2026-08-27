"""Empirical terminal-damage profiles from matched graph alternatives."""

import json
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.bidirected import BidirectedDeBruijnGraph, _successor
from anvaya.tip_matching import TipBackboneMatch


@dataclass(slots=True, frozen=True)
class DamageProfileBin:
    """Alternative and backbone observations at one terminal distance."""

    distance: int
    five_prime_ct_alternative: int
    five_prime_c_reference: int
    five_prime_ct_fraction: float | None
    three_prime_ga_alternative: int
    three_prime_g_reference: int
    three_prime_ga_fraction: float | None
    damage_alternative: int
    damage_reference: int
    damage_fraction: float | None
    alternative_edge_contributors: int
    reference_edge_contributors: int
    largest_alternative_edge_contribution: int
    largest_reference_edge_contribution: int


@dataclass(slots=True, frozen=True)
class DamageProfile:
    """Candidate-locus terminal substitution profile for one sample."""

    end_window: int
    matched_paths: int
    eligible_loci: int
    bins: tuple[DamageProfileBin, ...]


def _edge_links_at(
    graph: BidirectedDeBruijnGraph,
    start: int,
    edge_ids: tuple[int, ...],
    edge_index: int,
    expected_end: str,
):
    current = start
    for index, edge_id in enumerate(edge_ids):
        if index == edge_index:
            side = expected_end
            if current != graph.edge_sources[edge_id]:
                side = "right" if side == "left" else "left"
            return tuple(
                link
                for link in graph.molecule_end_links(edge_id)
                if link.end == side
            )
        current = _successor(graph, current, edge_id)
    return ()


def infer_damage_profile(
    graph: BidirectedDeBruijnGraph,
    matches: Sequence[TipBackboneMatch],
) -> DamageProfile:
    """Infer per-distance fractions at unique matched substitution loci."""
    if not graph.molecule_links_collected:
        raise ValueError("damage profiling requires molecule links")

    five_alternative = [0] * graph.end_window
    five_reference = [0] * graph.end_window
    three_alternative = [0] * graph.end_window
    three_reference = [0] * graph.end_window
    seen_loci: set[tuple[int, int, str, str]] = set()
    seen_alternatives: set[tuple[int, str, str]] = set()
    seen_references: set[tuple[int, str, str]] = set()
    alternative_contributors = [set() for _ in range(graph.end_window)]
    reference_contributors = [set() for _ in range(graph.end_window)]
    largest_alternative = [0] * graph.end_window
    largest_reference = [0] * graph.end_window

    for match in matches:
        for change in match.substitutions:
            edge_index = change.position - graph.node_length - 1
            if not (
                0 <= edge_index < len(match.tip_edge_ids)
                and edge_index < len(match.backbone_edge_ids)
            ):
                continue
            if change.reference == "C" and change.alternative == "T":
                expected_end = "left"
                alternative_counts = five_alternative
                reference_counts = five_reference
            elif change.reference == "G" and change.alternative == "A":
                expected_end = "right"
                alternative_counts = three_alternative
                reference_counts = three_reference
            else:
                continue
            locus = (
                match.tip_edge_ids[edge_index],
                match.backbone_edge_ids[edge_index],
                change.reference,
                change.alternative,
            )
            if locus in seen_loci:
                continue
            seen_loci.add(locus)
            alternative_key = (
                match.tip_edge_ids[edge_index],
                change.reference,
                change.alternative,
            )
            if alternative_key not in seen_alternatives:
                seen_alternatives.add(alternative_key)
                links = _edge_links_at(
                    graph,
                    match.branch_handle,
                    match.tip_edge_ids,
                    edge_index,
                    expected_end,
                )
                per_distance = [0] * graph.end_window
                for link in links:
                    alternative_counts[link.distance] += 1
                    per_distance[link.distance] += 1
                for distance, count in enumerate(per_distance):
                    if count:
                        alternative_contributors[distance].add(alternative_key)
                        largest_alternative[distance] = max(
                            largest_alternative[distance], count
                        )
            reference_key = (
                match.backbone_edge_ids[edge_index],
                change.reference,
                change.alternative,
            )
            if reference_key not in seen_references:
                seen_references.add(reference_key)
                links = _edge_links_at(
                    graph,
                    match.branch_handle,
                    match.backbone_edge_ids,
                    edge_index,
                    expected_end,
                )
                per_distance = [0] * graph.end_window
                for link in links:
                    reference_counts[link.distance] += 1
                    per_distance[link.distance] += 1
                for distance, count in enumerate(per_distance):
                    if count:
                        reference_contributors[distance].add(reference_key)
                        largest_reference[distance] = max(
                            largest_reference[distance], count
                        )

    bins = []
    for distance in range(graph.end_window):
        five_total = five_alternative[distance] + five_reference[distance]
        three_total = three_alternative[distance] + three_reference[distance]
        combined_alternative = (
            five_alternative[distance] + three_alternative[distance]
        )
        combined_reference = five_reference[distance] + three_reference[distance]
        combined_total = combined_alternative + combined_reference
        bins.append(
            DamageProfileBin(
                distance=distance,
                five_prime_ct_alternative=five_alternative[distance],
                five_prime_c_reference=five_reference[distance],
                five_prime_ct_fraction=(
                    five_alternative[distance] / five_total
                    if five_total
                    else None
                ),
                three_prime_ga_alternative=three_alternative[distance],
                three_prime_g_reference=three_reference[distance],
                three_prime_ga_fraction=(
                    three_alternative[distance] / three_total
                    if three_total
                    else None
                ),
                damage_alternative=combined_alternative,
                damage_reference=combined_reference,
                damage_fraction=(
                    combined_alternative / combined_total
                    if combined_total
                    else None
                ),
                alternative_edge_contributors=len(
                    alternative_contributors[distance]
                ),
                reference_edge_contributors=len(
                    reference_contributors[distance]
                ),
                largest_alternative_edge_contribution=(
                    largest_alternative[distance]
                ),
                largest_reference_edge_contribution=(
                    largest_reference[distance]
                ),
            )
        )
    return DamageProfile(
        end_window=graph.end_window,
        matched_paths=len(matches),
        eligible_loci=len(seen_loci),
        bins=tuple(bins),
    )


def write_damage_profile(profile: DamageProfile, path: str | Path) -> None:
    """Write a candidate-locus damage profile as JSON."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 2,
        "interpretation": "matched_candidate_locus_fraction",
        "end_window": profile.end_window,
        "matched_paths": profile.matched_paths,
        "eligible_loci": profile.eligible_loci,
        "bins": [asdict(value) for value in profile.bins],
    }
    output_path.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
