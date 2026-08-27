"""Validate incomplete-branch detection against controlled graph truth."""

import argparse
import csv
import random
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.bidirected import (
    _positioned_encoded_kmers,
    build_bidirected_dbg,
    spell_edge_path,
)
from anvaya.bubbles import find_simple_bubbles
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.incomplete_branches import (
    find_incomplete_branch_candidates,
    match_incomplete_branch_to_backbone,
)
from anvaya.sequences import reverse_complement
from anvaya.tip_classification import classify_tip_match
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
RARE_COVERAGE = 2
SNP_COUNT = 12
K = 21
MIN_COUNT = 2
END_WINDOW = 5
SEEDS = tuple(range(1501, 1521))
DAMAGE_PROFILE = (0.45, 0.30, 0.20, 0.12, 0.06)
ERROR_RATE = 0.010
MAX_BRANCH_EDGES = 64
MAXIMUM_RELATIVE_SUPPORT = 0.50


@dataclass(frozen=True)
class CandidateRecord:
    condition: str
    seed: int
    candidate_index: int
    truth_label: str
    edge_count: int
    relative_coverage: float
    sequence_identity: float
    ry_identity: float
    classification: str
    damage_score: float
    error_score: float
    variation_score: float
    substitutions: str


@dataclass(frozen=True)
class RunRecord:
    condition: str
    seed: int
    retained_truth_edges: int
    existing_truth_edges: int
    incomplete_truth_edges: int
    combined_truth_edges: int
    incomplete_branches: int
    truth_branches: int
    damage_like_branches: int
    truth_damage_like_branches: int
    false_damage_like_branches: int
    ambiguous_branches: int


RUN_FIELDS = tuple(RunRecord.__dataclass_fields__)
CANDIDATE_FIELDS = tuple(CandidateRecord.__dataclass_fields__)


def random_reference(seed: int) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(REFERENCE_LENGTH))


def rare_reference(reference: str, seed: int) -> str:
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
    return mutated


def canonical_code(sequence: str) -> int:
    return next(_positioned_encoded_kmers(sequence, len(sequence)))[1]


def mutated_kmer_codes(
    source_reads: list[SimulatedRead],
    changed_reads: list[SimulatedRead],
    position_attribute: str,
) -> set[int]:
    codes: set[int] = set()
    for source, changed in zip(source_reads, changed_reads, strict=True):
        positions = getattr(changed, position_attribute)
        for start, changed_code, _ in _positioned_encoded_kmers(
            changed.sequence,
            K,
        ):
            if not any(start <= position < start + K for position in positions):
                continue
            source_code = canonical_code(source.sequence[start : start + K])
            if changed_code != source_code:
                codes.add(changed_code)
    return codes


def graph_edge_codes(graph) -> dict[int, int]:
    return {
        canonical_code(
            spell_edge_path(graph, graph.edge_sources[edge_id], (edge_id,))
        ): edge_id
        for edge_id in range(graph.edge_count)
    }


def contains_sequence(reference: str, sequence: str) -> bool:
    return sequence in reference or reverse_complement(sequence) in reference


def strain_label(sequence: str, major: str, rare: str) -> str:
    in_major = contains_sequence(major, sequence)
    in_rare = contains_sequence(rare, sequence)
    if in_major and in_rare:
        return "shared"
    if in_major:
        return "major"
    if in_rare:
        return "rare"
    return "artifact"


def evaluate(
    condition: str,
    seed: int,
    reads: list[SimulatedRead],
    truth_codes: set[int] | None = None,
    major_reference: str | None = None,
    rare_reference_sequence: str | None = None,
) -> tuple[RunRecord, list[CandidateRecord]]:
    graph = build_bidirected_dbg(
        [read.sequence for read in reads],
        K,
        MIN_COUNT,
        END_WINDOW,
    )
    tips = find_weak_tip_candidates(graph)
    bubbles = find_simple_bubbles(graph)
    branches = find_incomplete_branch_candidates(
        graph,
        max_edges=MAX_BRANCH_EDGES,
        maximum_relative_support=MAXIMUM_RELATIVE_SUPPORT,
    )
    edge_codes = graph_edge_codes(graph)
    retained_truth = {
        edge_codes[code]
        for code in truth_codes or ()
        if code in edge_codes
    }
    existing_edges = {edge for tip in tips for edge in tip.edge_ids}
    existing_edges.update(
        edge
        for bubble in bubbles
        for path in bubble.paths
        for edge in path.edge_ids
    )
    branch_edges = {edge for branch in branches for edge in branch.edge_ids}

    candidates: list[CandidateRecord] = []
    for index, branch in enumerate(branches, start=1):
        match = match_incomplete_branch_to_backbone(graph, branch)
        if match is None:
            continue
        decision = classify_tip_match(graph, match)
        if major_reference is not None and rare_reference_sequence is not None:
            label = strain_label(
                match.tip_sequence,
                major_reference,
                rare_reference_sequence,
            )
        elif set(branch.edge_ids) & retained_truth:
            label = condition
        elif condition == "clean":
            label = "clean"
        else:
            label = "artifact"
        candidates.append(
            CandidateRecord(
                condition=condition,
                seed=seed,
                candidate_index=index,
                truth_label=label,
                edge_count=len(branch.edge_ids),
                relative_coverage=match.relative_coverage,
                sequence_identity=match.sequence_identity,
                ry_identity=match.ry_identity,
                classification=decision.label,
                damage_score=decision.damage_score,
                error_score=decision.error_score,
                variation_score=decision.variation_score,
                substitutions=";".join(
                    f"{change.position}:{change.reference}>{change.alternative}"
                    for change in match.substitutions
                ),
            )
        )

    truth_labels = {"damage", "error", "rare"}
    truth_candidates = [
        candidate
        for candidate in candidates
        if candidate.truth_label in truth_labels
    ]
    return RunRecord(
        condition=condition,
        seed=seed,
        retained_truth_edges=len(retained_truth),
        existing_truth_edges=len(existing_edges & retained_truth),
        incomplete_truth_edges=len(branch_edges & retained_truth),
        combined_truth_edges=len((existing_edges | branch_edges) & retained_truth),
        incomplete_branches=len(candidates),
        truth_branches=len(truth_candidates),
        damage_like_branches=sum(
            candidate.classification == "damage-like"
            for candidate in candidates
        ),
        truth_damage_like_branches=sum(
            candidate.classification == "damage-like"
            and candidate.truth_label == "damage"
            for candidate in candidates
        ),
        false_damage_like_branches=sum(
            candidate.classification == "damage-like"
            and candidate.truth_label != "damage"
            for candidate in candidates
        ),
        ambiguous_branches=sum(
            candidate.classification == "ambiguous"
            for candidate in candidates
        ),
    ), candidates


def run_matrix(
    seeds: tuple[int, ...],
) -> tuple[list[RunRecord], list[CandidateRecord]]:
    runs: list[RunRecord] = []
    candidates: list[CandidateRecord] = []
    for seed in seeds:
        reference = random_reference(seed)
        source = simulate_fragments(
            reference,
            READ_LENGTH,
            TOTAL_COVERAGE,
            seed + 10_000,
            reverse_fraction=0.5,
        )
        damage = add_terminal_damage(
            source,
            DAMAGE_PROFILE,
            DAMAGE_PROFILE,
            seed + 20_000,
        )
        errors = add_sequencing_errors(
            source,
            ERROR_RATE,
            seed + 30_000,
        )
        related = rare_reference(reference, seed + 50_000)
        major_reads = simulate_fragments(
            reference,
            READ_LENGTH,
            TOTAL_COVERAGE - RARE_COVERAGE,
            seed + 60_000,
            reverse_fraction=0.5,
        )
        rare_reads = simulate_fragments(
            related,
            READ_LENGTH,
            RARE_COVERAGE,
            seed + 70_000,
            reverse_fraction=0.5,
        )
        conditions = (
            ("clean", source, None, None, None),
            (
                "damage",
                damage,
                mutated_kmer_codes(source, damage, "damage_positions"),
                None,
                None,
            ),
            (
                "error",
                errors,
                mutated_kmer_codes(source, errors, "error_positions"),
                None,
                None,
            ),
            ("rare", major_reads + rare_reads, None, reference, related),
        )
        for condition, reads, truth, major, rare in conditions:
            run, found = evaluate(
                condition,
                seed,
                reads,
                truth_codes=truth,
                major_reference=major,
                rare_reference_sequence=rare,
            )
            runs.append(run)
            candidates.extend(found)
    return runs, candidates


def _write_dataclasses(
    path: Path,
    records: list[object],
    fields: tuple[str, ...],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)


def write_results(
    output_dir: Path,
    seeds: tuple[int, ...],
    runs: list[RunRecord],
    candidates: list[CandidateRecord],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_dataclasses(output_dir / "per_run.tsv", runs, RUN_FIELDS)
    _write_dataclasses(
        output_dir / "candidates.tsv",
        candidates,
        CANDIDATE_FIELDS,
    )
    with (output_dir / "summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("condition", "metric", "mean", "stdev", "minimum", "maximum", "n")
        )
        for condition in ("clean", "damage", "error", "rare"):
            selected = [run for run in runs if run.condition == condition]
            for metric in RUN_FIELDS[2:]:
                values = [getattr(run, metric) for run in selected]
                writer.writerow(
                    (
                        condition,
                        metric,
                        statistics.mean(values),
                        statistics.stdev(values) if len(values) > 1 else 0.0,
                        min(values),
                        max(values),
                        len(values),
                    )
                )
    with (output_dir / "parameters.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("parameter", "value"))
        writer.writerows(
            (
                ("reference_length", REFERENCE_LENGTH),
                ("read_length", READ_LENGTH),
                ("total_coverage", TOTAL_COVERAGE),
                ("rare_coverage", RARE_COVERAGE),
                ("k", K),
                ("min_count", MIN_COUNT),
                ("end_window", END_WINDOW),
                ("damage_profile", ",".join(map(str, DAMAGE_PROFILE))),
                ("error_rate", ERROR_RATE),
                ("max_branch_edges", MAX_BRANCH_EDGES),
                ("maximum_relative_support", MAXIMUM_RELATIVE_SUPPORT),
                ("seeds", ",".join(map(str, seeds))),
            )
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "experiments/validation/generated/incomplete_branch_validation"
        ),
    )
    parser.add_argument(
        "--seed-count",
        type=int,
        default=len(SEEDS),
        choices=range(1, len(SEEDS) + 1),
    )
    args = parser.parse_args()
    seeds = SEEDS[: args.seed_count]
    runs, candidates = run_matrix(seeds)
    write_results(args.output_dir, seeds, runs, candidates)
    print(f"runs={len(runs)}")
    print(f"candidates={len(candidates)}")
    print(f"output={args.output_dir}")


if __name__ == "__main__":
    main()
