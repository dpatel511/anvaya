"""Empirical terminal-damage profiles from matched graph alternatives."""

import json
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.bidirected import BidirectedDeBruijnGraph, _successor
from anvaya.damage_likelihood import (
    CandidateDamageFit,
    fit_candidate_damage_model,
)
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
    equal_locus_contributors: int
    equal_locus_fraction: float | None
    coverage_capped_weight: int
    coverage_capped_fraction: float | None


@dataclass(slots=True, frozen=True)
class DamageProfileLocus:
    """Distance-binned observations for one matched graph locus."""

    alternative_edge_id: int
    reference_edge_id: int
    channel: str
    alternative: tuple[int, ...]
    reference: tuple[int, ...]


@dataclass(slots=True, frozen=True)
class DamageProfile:
    """Candidate-locus terminal substitution profile for one sample."""

    end_window: int
    matched_paths: int
    matched_loci: int
    observed_loci: int
    coverage_cap: int
    bins: tuple[DamageProfileBin, ...]
    loci: tuple[DamageProfileLocus, ...]
    candidate_damage_fit: CandidateDamageFit

    @property
    def eligible_loci(self) -> int:
        """Backward-compatible name for loci with terminal observations."""
        return self.observed_loci


def _summarize_loci(
    loci: Sequence[DamageProfileLocus],
    distance: int,
    coverage_cap: int,
) -> tuple[int, float | None, int, float | None]:
    fractions = []
    capped_numerator = 0.0
    capped_weight = 0
    for locus in loci:
        alternative = locus.alternative[distance]
        reference = locus.reference[distance]
        total = alternative + reference
        if not total:
            continue
        fraction = alternative / total
        fractions.append(fraction)
        weight = min(total, coverage_cap)
        capped_numerator += fraction * weight
        capped_weight += weight
    return (
        len(fractions),
        sum(fractions) / len(fractions) if fractions else None,
        capped_weight,
        capped_numerator / capped_weight if capped_weight else None,
    )


def _edge_distances_at(
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
            offset = graph.node_length if expected_end == "left" else 0
            return tuple(
                distance + offset
                for distance in graph.molecule_end_distances(edge_id, side)
                if distance + offset < graph.end_window
            )
        current = _successor(graph, current, edge_id)
    return ()


def infer_damage_profile(
    graph: BidirectedDeBruijnGraph,
    matches: Sequence[TipBackboneMatch],
    coverage_cap: int = 20,
) -> DamageProfile:
    """Infer per-distance fractions at unique matched substitution loci."""
    if not graph.molecule_links_collected:
        raise ValueError("damage profiling requires molecule links")
    if not isinstance(coverage_cap, int) or isinstance(coverage_cap, bool):
        raise TypeError("coverage_cap must be an integer")
    if coverage_cap < 1:
        raise ValueError("coverage_cap must be at least 1")

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
    loci: list[DamageProfileLocus] = []

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
                channel = "five_prime_ct"
                alternative_counts = five_alternative
                reference_counts = five_reference
            elif change.reference == "G" and change.alternative == "A":
                expected_end = "right"
                channel = "three_prime_ga"
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
            alternative_distances = _edge_distances_at(
                graph,
                match.branch_handle,
                match.tip_edge_ids,
                edge_index,
                expected_end,
            )
            reference_distances = _edge_distances_at(
                graph,
                match.branch_handle,
                match.backbone_edge_ids,
                edge_index,
                expected_end,
            )
            locus_alternative = [0] * graph.end_window
            locus_reference = [0] * graph.end_window
            for distance in alternative_distances:
                locus_alternative[distance] += 1
            for distance in reference_distances:
                locus_reference[distance] += 1
            loci.append(
                DamageProfileLocus(
                    alternative_edge_id=match.tip_edge_ids[edge_index],
                    reference_edge_id=match.backbone_edge_ids[edge_index],
                    channel=channel,
                    alternative=tuple(locus_alternative),
                    reference=tuple(locus_reference),
                )
            )
            alternative_key = (
                match.tip_edge_ids[edge_index],
                change.reference,
                change.alternative,
            )
            if alternative_key not in seen_alternatives:
                seen_alternatives.add(alternative_key)
                for distance, count in enumerate(locus_alternative):
                    if count:
                        alternative_counts[distance] += count
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
                for distance, count in enumerate(locus_reference):
                    if count:
                        reference_counts[distance] += count
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
        (
            equal_locus_contributors,
            equal_locus_fraction,
            capped_weight,
            coverage_capped_fraction,
        ) = _summarize_loci(loci, distance, coverage_cap)
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
                equal_locus_contributors=equal_locus_contributors,
                equal_locus_fraction=equal_locus_fraction,
                coverage_capped_weight=capped_weight,
                coverage_capped_fraction=coverage_capped_fraction,
            )
        )
    locus_values = tuple(loci)
    return DamageProfile(
        end_window=graph.end_window,
        matched_paths=len(matches),
        matched_loci=len(seen_loci),
        observed_loci=sum(
            any(locus.alternative) or any(locus.reference)
            for locus in loci
        ),
        coverage_cap=coverage_cap,
        bins=tuple(bins),
        loci=locus_values,
        candidate_damage_fit=fit_candidate_damage_model(
            tuple((locus.alternative, locus.reference) for locus in locus_values),
            graph.end_window,
        ),
    )


def write_damage_profile(profile: DamageProfile, path: str | Path) -> None:
    """Write a candidate-locus damage profile as JSON."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 5,
        "interpretation": "matched_candidate_locus_fraction",
        "end_window": profile.end_window,
        "matched_paths": profile.matched_paths,
        "matched_loci": profile.matched_loci,
        "observed_loci": profile.observed_loci,
        "eligible_loci": profile.eligible_loci,
        "coverage_cap": profile.coverage_cap,
        "bins": [asdict(value) for value in profile.bins],
        "loci": [asdict(value) for value in profile.loci],
        "candidate_damage_fit": asdict(profile.candidate_damage_fit),
    }
    output_path.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
