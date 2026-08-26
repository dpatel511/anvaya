"""Validate terminal-event evidence against controlled simulated ground truth."""

import argparse
import csv
import random
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.bidirected import (
    _edge_support_total,
    _positioned_encoded_kmers,
    build_bidirected_dbg,
    extract_bidirected_unitigs,
    spell_edge_path,
    summarize_bidirected_graph,
)
from anvaya.bubbles import find_simple_bubbles
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.tip_matching import match_tip_to_backbone
from helpers.simulation import (
    SimulatedRead,
    add_sequencing_errors,
    add_terminal_damage,
    simulate_fragments,
)

REFERENCE_LENGTH = 20_000
READ_LENGTH = 60
COVERAGE = 20
K = 21
MIN_COUNT = 2
END_WINDOW = 5
SEEDS = tuple(range(901, 911))
FIVE_PRIME_CT = (0.45, 0.30, 0.20, 0.12, 0.06)
THREE_PRIME_GA = (0.45, 0.30, 0.20, 0.12, 0.06)
ERROR_RATE = 0.005


@dataclass(frozen=True)
class RunResult:
    condition: str
    seed: int
    introduced_events: int
    retained_truth_edges: int
    candidate_edges: int
    truth_candidate_edges: int
    nontruth_candidate_edges: int
    candidate_recall: float | None
    candidate_precision: float | None
    terminal_only_candidate_edges: int
    tips: int
    truth_tips: int
    matched_tips: int
    matched_truth_tips: int
    ry_exact_matches: int
    damage_compatible_matches: int
    damage_compatible_truth_matches: int
    bubbles: int
    branching_nodes: int
    unitigs: int
    longest_unitig: int
    terminal_fraction: float
    candidate_terminal_fraction: float | None


FIELDS = tuple(RunResult.__dataclass_fields__)


def random_reference(seed: int) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(REFERENCE_LENGTH))


def canonical_code(sequence: str) -> int:
    return next(_positioned_encoded_kmers(sequence, len(sequence)))[1]


def graph_edge_codes(graph) -> dict[int, int]:
    return {
        canonical_code(
            spell_edge_path(graph, graph.edge_sources[edge_id], (edge_id,))
        ): edge_id
        for edge_id in range(graph.edge_count)
    }


def mutated_kmer_codes(
    source_reads: list[SimulatedRead],
    changed_reads: list[SimulatedRead],
    position_attribute: str,
) -> set[int]:
    codes: set[int] = set()
    for source, changed in zip(source_reads, changed_reads, strict=True):
        positions = getattr(changed, position_attribute)
        for start, changed_code, _ in _positioned_encoded_kmers(changed.sequence, K):
            if not any(start <= position < start + K for position in positions):
                continue
            source_code = canonical_code(source.sequence[start : start + K])
            if changed_code != source_code:
                codes.add(changed_code)
    return codes


def candidate_edges(tips, bubbles) -> set[int]:
    edges = {edge for tip in tips for edge in tip.edge_ids}
    edges.update(
        edge
        for bubble in bubbles
        for path in bubble.paths
        for edge in path.edge_ids
    )
    return edges


def terminal_fraction(graph, edge_ids: set[int] | None = None) -> float | None:
    selected = range(graph.edge_count) if edge_ids is None else edge_ids
    terminal = 0
    observations = 0
    for edge_id in selected:
        support = graph.end_support(edge_id)
        terminal += sum(support.left) + sum(support.right) + support.ambiguous_terminal
        observations += _edge_support_total(graph, edge_id)
    return terminal / observations if observations else None


def terminal_only_edges(graph, edge_ids: set[int]) -> set[int]:
    selected: set[int] = set()
    for edge_id in edge_ids:
        support = graph.end_support(edge_id)
        terminal = sum(support.left) + sum(support.right) + support.ambiguous_terminal
        if terminal and support.internal == 0:
            selected.add(edge_id)
    return selected


def evaluate(condition, seed, reads, truth_codes, position_attribute) -> RunResult:
    graph = build_bidirected_dbg(
        [read.sequence for read in reads], K, MIN_COUNT, END_WINDOW
    )
    tips = find_weak_tip_candidates(graph)
    bubbles = find_simple_bubbles(graph)
    candidates = candidate_edges(tips, bubbles)
    codes = graph_edge_codes(graph)
    retained_truth = {codes[code] for code in truth_codes if code in codes}
    true_positives = candidates & retained_truth
    false_positives = candidates - retained_truth
    terminal_candidates = terminal_only_edges(graph, candidates)
    truth_tips = 0
    matched_tips = 0
    matched_truth_tips = 0
    ry_exact_matches = 0
    damage_compatible_matches = 0
    damage_compatible_truth_matches = 0
    for tip in tips:
        is_truth = bool(set(tip.edge_ids) & retained_truth)
        truth_tips += is_truth
        match = match_tip_to_backbone(graph, tip)
        if match is None:
            continue
        matched_tips += 1
        matched_truth_tips += is_truth
        ry_exact_matches += match.ry_identity == 1.0
        damage_compatible_matches += match.damage_compatible
        damage_compatible_truth_matches += match.damage_compatible and is_truth
    unitig_sequences = extract_bidirected_unitigs(graph)
    summary = summarize_bidirected_graph(graph, unitig_sequences)
    introduced = (
        sum(len(getattr(read, position_attribute)) for read in reads)
        if position_attribute
        else 0
    )
    return RunResult(
        condition, seed, introduced, len(retained_truth), len(candidates),
        len(true_positives), len(false_positives),
        len(true_positives) / len(retained_truth) if retained_truth else None,
        len(true_positives) / len(candidates) if candidates else None,
        len(terminal_candidates),
        len(tips), truth_tips, matched_tips, matched_truth_tips,
        ry_exact_matches, damage_compatible_matches,
        damage_compatible_truth_matches, len(bubbles),
        summary.branching_nodes, summary.unitigs,
        max(map(len, unitig_sequences)), terminal_fraction(graph) or 0.0,
        terminal_fraction(graph, candidates),
    )


def write_results(output_dir: Path, results: list[RunResult]) -> None:
    parameters = {
        "reference_length": REFERENCE_LENGTH,
        "read_length": READ_LENGTH,
        "coverage": COVERAGE,
        "k": K,
        "min_count": MIN_COUNT,
        "end_window": END_WINDOW,
        "seeds": ",".join(map(str, SEEDS)),
        "five_prime_ct": ",".join(map(str, FIVE_PRIME_CT)),
        "three_prime_ga": ",".join(map(str, THREE_PRIME_GA)),
        "error_rate": ERROR_RATE,
    }
    with (output_dir / "parameters.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("parameter", "value"))
        writer.writerows(parameters.items())

    with (output_dir / "per_run.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(result) for result in results)

    with (output_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("condition", "metric", "mean", "stdev", "minimum", "maximum", "n"))
        for condition in ("clean", "damage", "error"):
            selected = [result for result in results if result.condition == condition]
            for metric in FIELDS[2:]:
                values = [getattr(result, metric) for result in selected]
                numeric = [value for value in values if value is not None]
                if numeric:
                    writer.writerow((
                        condition, metric, statistics.mean(numeric),
                        statistics.stdev(numeric) if len(numeric) > 1 else 0.0,
                        min(numeric), max(numeric), len(numeric),
                    ))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("experiments/validation/generated/terminal_evidence"),
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    results: list[RunResult] = []
    for seed in SEEDS:
        source = simulate_fragments(
            random_reference(seed), READ_LENGTH, COVERAGE, seed + 10_000
        )
        damaged = add_terminal_damage(
            source, FIVE_PRIME_CT, THREE_PRIME_GA, seed + 20_000
        )
        errors = add_sequencing_errors(source, ERROR_RATE, seed + 30_000)
        conditions = (
            ("clean", source, set(), None),
            ("damage", damaged, mutated_kmer_codes(source, damaged, "damage_positions"), "damage_positions"),
            ("error", errors, mutated_kmer_codes(source, errors, "error_positions"), "error_positions"),
        )
        for condition, reads, truth, attribute in conditions:
            result = evaluate(condition, seed, reads, truth, attribute)
            results.append(result)
            print(
                f"{condition} seed={seed} tips={result.tips} bubbles={result.bubbles} "
                f"candidates={result.candidate_edges} truth={result.retained_truth_edges} "
                f"matched={result.truth_candidate_edges} "
                f"tip_matches={result.matched_tips} "
                f"damage_compatible={result.damage_compatible_matches} "
                f"terminal_only={result.terminal_only_candidate_edges}"
            )

    write_results(args.output_dir, results)
    print(f"runs={len(results)}")
    print(f"output={args.output_dir}")


if __name__ == "__main__":
    main()
