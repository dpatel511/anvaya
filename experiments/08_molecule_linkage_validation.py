"""Ablate molecule linkage across damage, error, and strain matrices."""

import argparse
import csv
import importlib.util
import math
import statistics
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.tip_classification import classify_tip_match
from anvaya.tip_matching import match_tip_to_backbone


def _load_tip_experiment():
    path = Path(__file__).with_name("06_tip_matching_validation.py")
    spec = importlib.util.spec_from_file_location("tip_matching_validation", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


previous = _load_tip_experiment()


@dataclass(frozen=True)
class CandidateRecord:
    condition: str
    scenario: str
    seed: int
    candidate_type: str
    candidate_index: int
    truth_label: str
    substitution_count: int
    joint_molecule_observations: int
    joint_molecule_fraction: float
    aggregate_label: str
    linked_label: str
    aggregate_damage_score: float
    linked_damage_score: float


@dataclass(frozen=True)
class RunRecord:
    condition: str
    scenario: str
    seed: int
    candidates: int
    multi_substitution_candidates: int
    label_changes: int
    aggregate_build_seconds: float
    linked_build_seconds: float
    molecule_link_observations: int
    molecule_link_bytes: int


def clustered_rare_reference(reference: str) -> str:
    """Introduce one nearby C>T and G>A pair as genuine linked variation."""
    margin = 2 * previous.READ_LENGTH
    minimum_gap = previous.READ_LENGTH - 2 * previous.END_WINDOW
    maximum_gap = previous.READ_LENGTH - 1
    for first in range(margin, len(reference) - margin - maximum_gap):
        if reference[first] != "C":
            continue
        for second in range(first + minimum_gap, first + maximum_gap + 1):
            if reference[second] == "G":
                changed = previous.introduce_snp(reference, first, "T")
                return previous.introduce_snp(changed, second, "A")
    raise RuntimeError("reference lacks a suitable clustered C/G pair")


def add_paired_terminal_errors(reads, probability: float, seed: int):
    """Add correlated ordinary C>T/G>A errors at opposite read ends."""
    generator = previous.random.Random(seed)
    changed = []
    for read in reads:
        bases = list(read.sequence)
        positions = set(read.error_positions)
        if generator.random() < probability:
            left = next(
                (
                    index
                    for index in range(previous.END_WINDOW)
                    if bases[index] == "C"
                ),
                None,
            )
            right = next(
                (
                    len(bases) - 1 - distance
                    for distance in range(previous.END_WINDOW)
                    if bases[len(bases) - 1 - distance] == "G"
                ),
                None,
            )
            if left is not None and right is not None:
                bases[left] = "T"
                bases[right] = "A"
                positions.update((left, right))
        changed.append(
            previous.replace(
                read,
                sequence="".join(bases),
                error_positions=tuple(sorted(positions)),
            )
        )
    return changed


def _condition_inputs(seed: int):
    reference = previous.random_reference(seed)
    source = previous.simulate_fragments(
        reference,
        previous.READ_LENGTH,
        previous.TOTAL_COVERAGE,
        seed + 10_000,
        reverse_fraction=0.5,
    )
    conditions = [("clean", "clean", source, None, None, None)]
    zero_profile = (0.0,) * previous.END_WINDOW
    for profile_name, profile in previous.DAMAGE_PROFILES.items():
        for end_name, five_prime, three_prime in (
            ("5p", profile, zero_profile),
            ("3p", zero_profile, profile),
            ("both", profile, profile),
        ):
            damaged = previous.add_terminal_damage(
                source, five_prime, three_prime, seed + 20_000 + len(conditions)
            )
            conditions.append(
                (
                    "damage",
                    f"damage_{profile_name}_{end_name}",
                    damaged,
                    previous.mutated_kmer_codes(
                        source, damaged, "damage_positions"
                    ),
                    None,
                    None,
                )
            )
    for error_rate in previous.ERROR_RATES:
        errors = previous.add_sequencing_errors(
            source,
            error_rate,
            seed + 30_000 + int(error_rate * 1000),
        )
        conditions.append(
            (
                "error",
                f"error_{error_rate:.3f}",
                errors,
                previous.mutated_kmer_codes(source, errors, "error_positions"),
                None,
                None,
            )
        )
    terminal_errors = previous.add_terminal_errors(
        source, previous.TERMINAL_ERROR_RATE, seed + 40_000
    )
    conditions.append(
        (
            "error",
            "error_terminal",
            terminal_errors,
            previous.mutated_kmer_codes(
                source, terminal_errors, "error_positions"
            ),
            None,
            None,
        )
    )
    paired_terminal_errors = add_paired_terminal_errors(
        source, previous.TERMINAL_ERROR_RATE, seed + 45_000
    )
    conditions.append(
        (
            "error",
            "error_terminal_paired",
            paired_terminal_errors,
            previous.mutated_kmer_codes(
                source, paired_terminal_errors, "error_positions"
            ),
            None,
            None,
        )
    )
    related, _ = previous.rare_reference(reference, seed + 50_000)
    clustered = clustered_rare_reference(reference)
    for rare_coverage in previous.RARE_COVERAGES:
        major_reads = previous.simulate_fragments(
            reference,
            previous.READ_LENGTH,
            previous.TOTAL_COVERAGE - rare_coverage,
            seed + 60_000 + rare_coverage,
            reverse_fraction=0.5,
        )
        for scenario, rare_reference, offset in (
            (f"rare_{rare_coverage}x", related, 0),
            (f"rare_clustered_{rare_coverage}x", clustered, 100),
        ):
            rare_reads = previous.simulate_fragments(
                rare_reference,
                previous.READ_LENGTH,
                rare_coverage,
                seed + 70_000 + rare_coverage + offset,
                reverse_fraction=0.5,
            )
            conditions.append(
                (
                    "rare",
                    scenario,
                    major_reads + rare_reads,
                    None,
                    reference,
                    rare_reference,
                )
            )
    return conditions


def _matched_candidates(graph):
    for index, tip in enumerate(find_weak_tip_candidates(graph), start=1):
        match = match_tip_to_backbone(graph, tip)
        if match is not None:
            yield "tip", index, tip.edge_ids, match


def evaluate(condition, scenario, seed, reads, truth_codes, major, rare):
    sequences = [read.sequence for read in reads]
    started = time.perf_counter()
    aggregate_graph = build_bidirected_dbg(
        sequences, previous.K, previous.MIN_COUNT, previous.END_WINDOW
    )
    aggregate_seconds = time.perf_counter() - started
    started = time.perf_counter()
    graph = build_bidirected_dbg(
        sequences,
        previous.K,
        previous.MIN_COUNT,
        previous.END_WINDOW,
        track_molecule_links=True,
    )
    linked_seconds = time.perf_counter() - started
    if (
        aggregate_graph.node_count != graph.node_count
        or aggregate_graph.edge_count != graph.edge_count
        or aggregate_graph.observations != graph.observations
    ):
        raise AssertionError("molecule tracking changed graph construction")

    edge_codes = previous.graph_edge_codes(graph)
    retained_truth = {
        edge_codes[code] for code in truth_codes or () if code in edge_codes
    }
    candidates = []
    for candidate_type, index, edge_ids, match in _matched_candidates(graph):
        aggregate = classify_tip_match(
            graph, match, include_molecule_linkage=False
        )
        linked = classify_tip_match(graph, match)
        if major is not None and rare is not None:
            truth_label = previous.strain_label(match.tip_sequence, major, rare)
        elif set(edge_ids) & retained_truth:
            truth_label = condition
        elif condition == "clean":
            truth_label = "clean"
        else:
            truth_label = "artifact"
        candidates.append(
            CandidateRecord(
                condition=condition,
                scenario=scenario,
                seed=seed,
                candidate_type=candidate_type,
                candidate_index=index,
                truth_label=truth_label,
                substitution_count=len(match.substitutions),
                joint_molecule_observations=(
                    linked.evidence.joint_molecule_observations or 0
                ),
                joint_molecule_fraction=(
                    linked.evidence.joint_molecule_fraction or 0.0
                ),
                aggregate_label=aggregate.label,
                linked_label=linked.label,
                aggregate_damage_score=aggregate.damage_score,
                linked_damage_score=linked.damage_score,
            )
        )
    return RunRecord(
        condition=condition,
        scenario=scenario,
        seed=seed,
        candidates=len(candidates),
        multi_substitution_candidates=sum(
            candidate.substitution_count > 1 for candidate in candidates
        ),
        label_changes=sum(
            candidate.aggregate_label != candidate.linked_label
            for candidate in candidates
        ),
        aggregate_build_seconds=aggregate_seconds,
        linked_build_seconds=linked_seconds,
        molecule_link_observations=len(graph.molecule_link_tokens),
        molecule_link_bytes=(
            len(graph.molecule_link_offsets) * graph.molecule_link_offsets.itemsize
            + len(graph.molecule_link_tokens) * graph.molecule_link_tokens.itemsize
        ),
    ), candidates


def exact_mcnemar_pvalue(aggregate_correct: int, linked_correct: int) -> float:
    """Two-sided exact McNemar p-value from discordant pair counts."""
    discordant = aggregate_correct + linked_correct
    if discordant == 0:
        return 1.0
    smaller = min(aggregate_correct, linked_correct)
    tail = sum(math.comb(discordant, index) for index in range(smaller + 1))
    return min(1.0, 2.0 * tail / (2**discordant))


def summarize(candidates):
    rows = []
    for truth_label in sorted({candidate.truth_label for candidate in candidates}):
        selected = [c for c in candidates if c.truth_label == truth_label]
        target = truth_label == "damage"
        aggregate_correct_only = sum(
            (c.aggregate_label == "damage-like") == target
            and (c.linked_label == "damage-like") != target
            for c in selected
        )
        linked_correct_only = sum(
            (c.aggregate_label == "damage-like") != target
            and (c.linked_label == "damage-like") == target
            for c in selected
        )
        rows.append(
            (
                truth_label,
                len(selected),
                sum(c.substitution_count > 1 for c in selected),
                sum(c.aggregate_label == "damage-like" for c in selected),
                sum(c.linked_label == "damage-like" for c in selected),
                sum(c.aggregate_label != c.linked_label for c in selected),
                statistics.mean(
                    c.linked_damage_score - c.aggregate_damage_score
                    for c in selected
                ),
                aggregate_correct_only,
                linked_correct_only,
                exact_mcnemar_pvalue(
                    aggregate_correct_only, linked_correct_only
                ),
            )
        )
    return rows


def _write_records(path, records):
    fields = tuple(records[0].__dataclass_fields__) if records else ()
    with path.open("w", encoding="utf-8", newline="") as handle:
        if not fields:
            return
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "experiments/validation/generated/molecule_linkage_validation"
        ),
    )
    parser.add_argument(
        "--seed-count",
        type=int,
        default=len(previous.SEEDS),
        choices=range(1, len(previous.SEEDS) + 1),
    )
    arguments = parser.parse_args()
    seeds = previous.SEEDS[: arguments.seed_count]
    runs = []
    candidates = []
    for seed in seeds:
        for values in _condition_inputs(seed):
            run, found = evaluate(*values[:2], seed, *values[2:])
            runs.append(run)
            candidates.extend(found)

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    _write_records(arguments.output_dir / "per_run.tsv", runs)
    _write_records(arguments.output_dir / "candidates.tsv", candidates)
    with (arguments.output_dir / "summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "truth_label",
                "candidates",
                "multi_substitution_candidates",
                "aggregate_damage_like",
                "linked_damage_like",
                "label_changes",
                "mean_damage_score_change",
                "aggregate_correct_only",
                "linked_correct_only",
                "mcnemar_exact_pvalue",
            )
        )
        writer.writerows(summarize(candidates))
    print(f"runs={len(runs)}")
    print(f"candidates={len(candidates)}")
    print(f"label_changes={sum(run.label_changes for run in runs)}")
    print(f"output={arguments.output_dir}")


if __name__ == "__main__":
    main()
