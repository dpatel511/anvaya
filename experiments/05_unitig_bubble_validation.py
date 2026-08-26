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
MAJOR_COVERAGE = 15
RARE_COVERAGE = 5
SNP_COUNT = 12
K = 21
MIN_COUNT = 2
SEEDS = tuple(range(1201, 1206))
FIVE_PRIME_CT = (0.45, 0.30, 0.20, 0.12, 0.06)
THREE_PRIME_GA = (0.45, 0.30, 0.20, 0.12, 0.06)
ERROR_RATE = 0.005
COVERAGE_THRESHOLDS = (0.10, 0.20, 0.30, 0.50)
SIMILARITY_THRESHOLDS = (0.90, 0.95, 0.98)


@dataclass(frozen=True)
class CandidateRecord:
    condition: str
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
    relative_coverage: float
    similarity_to_strongest: float


@dataclass(frozen=True)
class RunRecord:
    condition: str
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
    seed: int
    maximum_coverage: float
    minimum_similarity: float
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
) -> bool:
    """Apply one prospective rule without changing the graph."""
    return (
        candidate.relative_coverage <= maximum_coverage
        and candidate.similarity_to_strongest >= minimum_similarity
    )


def evaluate(
    condition: str,
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
                    relative_coverage=path.relative_coverage,
                    similarity_to_strongest=path.similarity_to_strongest,
                )
            )

    labels = Counter(candidate.label for candidate in candidates)
    run = RunRecord(
        condition=condition,
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
            selected = [
                candidate
                for candidate in candidates
                if selected_by_threshold(
                    candidate,
                    maximum_coverage,
                    minimum_similarity,
                )
            ]
            selected_labels = Counter(candidate.label for candidate in selected)
            thresholds.append(
                ThresholdRecord(
                    condition=condition,
                    seed=seed,
                    maximum_coverage=maximum_coverage,
                    minimum_similarity=minimum_similarity,
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
        "major_coverage": MAJOR_COVERAGE,
        "rare_coverage": RARE_COVERAGE,
        "snp_count": SNP_COUNT,
        "k": K,
        "min_count": MIN_COUNT,
        "seeds": ",".join(map(str, seeds)),
        "five_prime_ct": ",".join(map(str, FIVE_PRIME_CT)),
        "three_prime_ga": ",".join(map(str, THREE_PRIME_GA)),
        "error_rate": ERROR_RATE,
        "coverage_thresholds": ",".join(map(str, COVERAGE_THRESHOLDS)),
        "similarity_thresholds": ",".join(map(str, SIMILARITY_THRESHOLDS)),
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
                "maximum_coverage",
                "minimum_similarity",
                "metric",
                "mean",
                "minimum",
                "maximum",
                "n",
            )
        )
        for condition in ("clean", "error", "damage", "rare"):
            for maximum_coverage in COVERAGE_THRESHOLDS:
                for minimum_similarity in SIMILARITY_THRESHOLDS:
                    selected = [
                        record
                        for record in thresholds
                        if record.condition == condition
                        and record.maximum_coverage == maximum_coverage
                        and record.minimum_similarity == minimum_similarity
                    ]
                    for metric in tuple(ThresholdRecord.__dataclass_fields__)[4:]:
                        values = [getattr(record, metric) for record in selected]
                        writer.writerow(
                            (
                                condition,
                                maximum_coverage,
                                minimum_similarity,
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
        minor_reference, snp_positions = rare_reference(reference, seed + 1_000)
        source = simulate_fragments(
            reference, READ_LENGTH, TOTAL_COVERAGE, seed + 10_000
        )
        errors = add_sequencing_errors(source, ERROR_RATE, seed + 20_000)
        damage = add_terminal_damage(
            source,
            FIVE_PRIME_CT,
            THREE_PRIME_GA,
            seed + 30_000,
        )
        major_reads = simulate_fragments(
            reference, READ_LENGTH, MAJOR_COVERAGE, seed + 40_000
        )
        rare_reads = simulate_fragments(
            minor_reference, READ_LENGTH, RARE_COVERAGE, seed + 50_000
        )
        conditions = (
            ("clean", source, None, 0),
            ("error", errors, None, event_count(errors, "error_positions")),
            ("damage", damage, None, event_count(damage, "damage_positions")),
            ("rare", major_reads + rare_reads, minor_reference, len(snp_positions)),
        )

        for condition, reads, minor, introduced in conditions:
            run, run_candidates, run_thresholds = evaluate(
                condition,
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
                f"{condition} seed={seed} bubbles={run.bubbles} "
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
