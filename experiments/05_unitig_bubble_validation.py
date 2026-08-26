"""Validate unitig-bubble scores on labelled simulated alternatives."""

import argparse
import csv
import random
import statistics
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cleaning import remove_weak_tips
from anvaya.sequences import reverse_complement
from anvaya.unitig_bubbles import find_unitig_bubbles, spell_unitig_path
from anvaya.unitig_graph import build_compacted_unitig_graph
from helpers.simulation import (
    SimulatedRead,
    add_sequencing_errors,
    add_terminal_damage,
    introduce_snp,
    simulate_fragments,
)

REFERENCE_LENGTH = 20_000
READ_LENGTH = 60
TOTAL_COVERAGE = 20
SNP_COUNT = 12
K = 21
MIN_COUNT = 2
END_WINDOW = 5
SEEDS = tuple(range(1201, 1221))
ERROR_RATES = (0.001, 0.005, 0.010)
RARE_COVERAGES = (1, 2, 5, 10)
DAMAGE_PROFILES = {
    "low": (0.15, 0.10, 0.06, 0.03, 0.01),
    "standard": (0.45, 0.30, 0.20, 0.12, 0.06),
    "high": (0.70, 0.50, 0.35, 0.20, 0.10),
}
COVERAGE_THRESHOLDS = (0.10, 0.20, 0.30, 0.50)
SIMILARITY_THRESHOLDS = (0.90, 0.95, 0.98)
STRAND_BALANCE_THRESHOLDS = (0.25, 0.50, 1.00)
TERMINAL_ENRICHMENT_THRESHOLDS = (-0.20, -0.15, -0.10, -0.05, 1.00)


@dataclass(frozen=True)
class CandidateRecord:
    condition: str
    scenario: str
    seed: int
    bubble_index: int
    path_index: int
    label: str
    sequence: str
    nucleotide_length: int
    edge_count: int
    observations: int
    minimum_edge_support: int
    mean_edge_support: float
    terminal_fraction: float
    strand_balance: float | None
    relative_coverage: float
    relative_internal_coverage: float | None
    terminal_enrichment: float
    similarity_to_strongest: float


@dataclass(frozen=True)
class RunRecord:
    condition: str
    scenario: str
    seed: int
    reads: int
    introduced_events: int
    tips_removed: int
    nodes: int
    edges: int
    unitigs: int
    bubbles: int
    alternatives: int
    major_paths: int
    rare_paths: int
    error_paths: int
    damage_paths: int
    artifact_paths: int
    shared_paths: int


@dataclass(frozen=True)
class ThresholdRecord:
    condition: str
    scenario: str
    seed: int
    maximum_coverage: float
    minimum_similarity: float
    maximum_strand_balance: float
    maximum_terminal_enrichment: float
    eligible_paths: int
    selected_paths: int
    selected_major: int
    selected_rare: int
    selected_error: int
    selected_damage: int
    selected_artifact: int
    selected_shared: int


def random_reference(seed: int) -> str:
    """Generate one deterministic unambiguous reference."""
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(REFERENCE_LENGTH))


def rare_reference(reference: str, seed: int) -> tuple[str, tuple[int, ...]]:
    """Introduce evenly separated SNPs into one related-strain reference."""
    generator = random.Random(seed)
    margin = 2 * READ_LENGTH
    usable = len(reference) - 2 * margin
    positions = tuple(
        margin + (index + 1) * usable // (SNP_COUNT + 1)
        for index in range(SNP_COUNT)
    )
    mutated = reference
    for position in positions:
        alternatives = [base for base in "ACGT" if base != mutated[position]]
        mutated = introduce_snp(mutated, position, generator.choice(alternatives))
    return mutated, positions


def contains_sequence(reference: str, sequence: str) -> bool:
    """Return whether either orientation of a path occurs in a reference."""
    return sequence in reference or reverse_complement(sequence) in reference


def label_path(
    condition: str,
    sequence: str,
    major_reference: str,
    minor_reference: str | None,
) -> str:
    """Assign reference-derived ground truth to one weaker bubble path."""
    in_major = contains_sequence(major_reference, sequence)
    in_minor = (
        minor_reference is not None
        and contains_sequence(minor_reference, sequence)
    )
    if in_major and in_minor:
        return "shared"
    if in_major:
        return "major"
    if in_minor:
        return "rare"
    if condition == "error":
        return "error"
    if condition == "damage":
        return "damage"
    return "artifact"


def selected_by_threshold(
    candidate: CandidateRecord,
    maximum_coverage: float,
    minimum_similarity: float,
    maximum_strand_balance: float,
    maximum_terminal_enrichment: float,
) -> bool:
    """Apply one prospective rule without changing the graph."""
    return (
        candidate.relative_coverage <= maximum_coverage
        and candidate.similarity_to_strongest >= minimum_similarity
        and candidate.strand_balance is not None
        and candidate.strand_balance <= maximum_strand_balance
        and candidate.terminal_enrichment <= maximum_terminal_enrichment
    )


def evaluate(
    condition: str,
    scenario: str,
    seed: int,
    reads: list[SimulatedRead],
    major_reference: str,
    minor_reference: str | None,
    introduced_events: int,
) -> tuple[RunRecord, list[CandidateRecord], list[ThresholdRecord]]:
    """Build, clean, compact, and score one labelled condition."""
    graph = build_bidirected_dbg(
        [read.sequence for read in reads],
        K,
        MIN_COUNT,
        END_WINDOW,
    )
    cleaning = remove_weak_tips(graph)
    compacted = build_compacted_unitig_graph(graph)
    bubbles = find_unitig_bubbles(compacted)
    candidates: list[CandidateRecord] = []

    for bubble_index, bubble in enumerate(bubbles, start=1):
        for path_index, path in enumerate(bubble.paths):
            if path_index == bubble.strongest_path_index:
                continue
            sequence = spell_unitig_path(compacted, path.handles)
            candidates.append(
                CandidateRecord(
                    condition=condition,
                    scenario=scenario,
                    seed=seed,
                    bubble_index=bubble_index,
                    path_index=path_index,
                    label=label_path(
                        condition,
                        sequence,
                        major_reference,
                        minor_reference,
                    ),
                    sequence=sequence,
                    nucleotide_length=path.nucleotide_length,
                    edge_count=path.edge_count,
                    observations=path.observations,
                    minimum_edge_support=path.minimum_edge_support,
                    mean_edge_support=path.mean_edge_support,
                    terminal_fraction=path.terminal_fraction,
                    strand_balance=path.strand_balance,
                    relative_coverage=path.relative_coverage,
                    relative_internal_coverage=path.relative_internal_coverage,
                    terminal_enrichment=path.terminal_enrichment,
                    similarity_to_strongest=path.similarity_to_strongest,
                )
            )

    labels = Counter(candidate.label for candidate in candidates)
    run = RunRecord(
        condition=condition,
        scenario=scenario,
        seed=seed,
        reads=len(reads),
        introduced_events=introduced_events,
        tips_removed=cleaning.tips_removed,
        nodes=graph.node_count,
        edges=graph.edge_count,
        unitigs=compacted.unitig_count,
        bubbles=len(bubbles),
        alternatives=len(candidates),
        major_paths=labels["major"],
        rare_paths=labels["rare"],
        error_paths=labels["error"],
        damage_paths=labels["damage"],
        artifact_paths=labels["artifact"],
        shared_paths=labels["shared"],
    )

    thresholds: list[ThresholdRecord] = []
    for maximum_coverage in COVERAGE_THRESHOLDS:
        for minimum_similarity in SIMILARITY_THRESHOLDS:
            for maximum_strand_balance in STRAND_BALANCE_THRESHOLDS:
                for maximum_terminal_enrichment in (
                    TERMINAL_ENRICHMENT_THRESHOLDS
                ):
                    selected = [
                        candidate
                        for candidate in candidates
                        if selected_by_threshold(
                            candidate,
                            maximum_coverage,
                            minimum_similarity,
                            maximum_strand_balance,
                            maximum_terminal_enrichment,
                        )
                    ]
                    selected_labels = Counter(
                        candidate.label for candidate in selected
                    )
                    thresholds.append(
                        ThresholdRecord(
                            condition=condition,
                            scenario=scenario,
                            seed=seed,
                            maximum_coverage=maximum_coverage,
                            minimum_similarity=minimum_similarity,
                            maximum_strand_balance=maximum_strand_balance,
                            maximum_terminal_enrichment=(
                                maximum_terminal_enrichment
                            ),
                            eligible_paths=len(candidates),
                            selected_paths=len(selected),
                            selected_major=selected_labels["major"],
                            selected_rare=selected_labels["rare"],
                            selected_error=selected_labels["error"],
                            selected_damage=selected_labels["damage"],
                            selected_artifact=selected_labels["artifact"],
                            selected_shared=selected_labels["shared"],
                        )
                    )

    return run, candidates, thresholds


def _write_dataclasses(path: Path, records: list[object], fields: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)


def write_results(
    output_dir: Path,
    seeds: tuple[int, ...],
    runs: list[RunRecord],
    candidates: list[CandidateRecord],
    thresholds: list[ThresholdRecord],
) -> None:
    """Write parameters, candidate evidence, and threshold summaries."""
    output_dir.mkdir(parents=True, exist_ok=True)
    parameters = {
        "reference_length": REFERENCE_LENGTH,
        "read_length": READ_LENGTH,
        "total_coverage": TOTAL_COVERAGE,
        "rare_coverages": ",".join(map(str, RARE_COVERAGES)),
        "snp_count": SNP_COUNT,
        "k": K,
        "min_count": MIN_COUNT,
        "end_window": END_WINDOW,
        "seeds": ",".join(map(str, seeds)),
        "damage_profiles": ";".join(
            name + ":" + ",".join(map(str, profile))
            for name, profile in DAMAGE_PROFILES.items()
        ),
        "error_rates": ",".join(map(str, ERROR_RATES)),
        "coverage_thresholds": ",".join(map(str, COVERAGE_THRESHOLDS)),
        "similarity_thresholds": ",".join(map(str, SIMILARITY_THRESHOLDS)),
        "strand_balance_thresholds": ",".join(
            map(str, STRAND_BALANCE_THRESHOLDS)
        ),
        "terminal_enrichment_thresholds": ",".join(
            map(str, TERMINAL_ENRICHMENT_THRESHOLDS)
        ),
    }
    with (output_dir / "parameters.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("parameter", "value"))
        writer.writerows(parameters.items())

    _write_dataclasses(
        output_dir / "per_run.tsv",
        list(runs),
        tuple(RunRecord.__dataclass_fields__),
    )
    _write_dataclasses(
        output_dir / "candidates.tsv",
        list(candidates),
        tuple(CandidateRecord.__dataclass_fields__),
    )
    _write_dataclasses(
        output_dir / "thresholds.tsv",
        list(thresholds),
        tuple(ThresholdRecord.__dataclass_fields__),
    )

    with (output_dir / "summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "condition",
                "scenario",
                "maximum_coverage",
                "minimum_similarity",
                "maximum_strand_balance",
                "maximum_terminal_enrichment",
                "metric",
                "mean",
                "minimum",
                "maximum",
                "n",
            )
        )
        scenarios = sorted({(record.condition, record.scenario) for record in thresholds})
        for condition, scenario in scenarios:
            for maximum_coverage in COVERAGE_THRESHOLDS:
                for minimum_similarity in SIMILARITY_THRESHOLDS:
                    for maximum_strand_balance in STRAND_BALANCE_THRESHOLDS:
                        for maximum_terminal_enrichment in (
                            TERMINAL_ENRICHMENT_THRESHOLDS
                        ):
                            selected = [
                                record
                                for record in thresholds
                                if record.condition == condition
                                and record.scenario == scenario
                                and record.maximum_coverage == maximum_coverage
                                and record.minimum_similarity == minimum_similarity
                                and record.maximum_strand_balance
                                == maximum_strand_balance
                                and record.maximum_terminal_enrichment
                                == maximum_terminal_enrichment
                            ]
                            for metric in tuple(
                                ThresholdRecord.__dataclass_fields__
                            )[7:]:
                                values = [
                                    getattr(record, metric) for record in selected
                                ]
                                writer.writerow(
                                    (
                                        condition,
                                        scenario,
                                        maximum_coverage,
                                        minimum_similarity,
                                        maximum_strand_balance,
                                        maximum_terminal_enrichment,
                                        metric,
                                        statistics.mean(values),
                                        min(values),
                                        max(values),
                                        len(values),
                                    )
                                )


def event_count(reads: list[SimulatedRead], attribute: str) -> int:
    return sum(len(getattr(read, attribute)) for read in reads)


def run_matrix(
    seeds: tuple[int, ...],
) -> tuple[list[RunRecord], list[CandidateRecord], list[ThresholdRecord]]:
    """Run every matched simulation condition."""
    runs: list[RunRecord] = []
    candidates: list[CandidateRecord] = []
    thresholds: list[ThresholdRecord] = []

    for seed in seeds:
        reference = random_reference(seed)
        source = simulate_fragments(
            reference,
            READ_LENGTH,
            TOTAL_COVERAGE,
            seed + 10_000,
            reverse_fraction=0.5,
        )
        conditions: list[tuple[str, str, list[SimulatedRead], str | None, int]] = [
            ("clean", "clean", source, None, 0)
        ]
        for index, error_rate in enumerate(ERROR_RATES):
            errors = add_sequencing_errors(
                source,
                error_rate,
                seed + 20_000 + index,
            )
            conditions.append(
                (
                    "error",
                    f"error_{error_rate:.3f}",
                    errors,
                    None,
                    event_count(errors, "error_positions"),
                )
            )
        for index, (name, profile) in enumerate(DAMAGE_PROFILES.items()):
            damage = add_terminal_damage(
                source,
                profile,
                profile,
                seed + 30_000 + index,
            )
            conditions.append(
                (
                    "damage",
                    f"damage_{name}",
                    damage,
                    None,
                    event_count(damage, "damage_positions"),
                )
            )
        minor_reference, snp_positions = rare_reference(reference, seed + 1_000)
        for rare_coverage in RARE_COVERAGES:
            major_coverage = TOTAL_COVERAGE - rare_coverage
            major_reads = simulate_fragments(
                reference,
                READ_LENGTH,
                major_coverage,
                seed + 40_000 + rare_coverage,
                reverse_fraction=0.5,
            )
            rare_reads = simulate_fragments(
                minor_reference,
                READ_LENGTH,
                rare_coverage,
                seed + 50_000 + rare_coverage,
                reverse_fraction=0.5,
            )
            conditions.append(
                (
                    "rare",
                    f"rare_{rare_coverage}x",
                    major_reads + rare_reads,
                    minor_reference,
                    len(snp_positions),
                )
            )

        for condition, scenario, reads, minor, introduced in conditions:
            run, run_candidates, run_thresholds = evaluate(
                condition,
                scenario,
                seed,
                reads,
                reference,
                minor,
                introduced,
            )
            runs.append(run)
            candidates.extend(run_candidates)
            thresholds.extend(run_thresholds)
            print(
                f"{scenario} seed={seed} bubbles={run.bubbles} "
                f"alternatives={run.alternatives} error={run.error_paths} "
                f"damage={run.damage_paths} rare={run.rare_paths}"
            )

    return runs, candidates, thresholds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("experiments/validation/generated/unitig_bubbles"),
    )
    parser.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    args = parser.parse_args()
    seeds = tuple(args.seeds)
    if not seeds:
        parser.error("at least one seed is required")

    runs, candidates, thresholds = run_matrix(seeds)
    write_results(args.output_dir, seeds, runs, candidates, thresholds)
    print(f"runs={len(runs)}")
    print(f"candidates={len(candidates)}")
    print(f"output={args.output_dir}")


if __name__ == "__main__":
    main()
