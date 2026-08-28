"""Stress-test calibrated confidence on truth-labelled graph events."""

import argparse
import csv
import importlib.util
import math
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg, spell_edge_path
from anvaya.bubbles import find_simple_bubbles
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.damage_profile import infer_damage_profile
from anvaya.event_calibration import (
    CalibrationExample,
    calibrate_score_batch,
    fit_conformal_calibration,
    scores_from_likelihood,
)
from anvaya.event_likelihood import (
    fit_cross_fitted_damage_models,
    score_matched_event,
)
from anvaya.incomplete_branches import (
    find_incomplete_branch_candidates,
    match_incomplete_branch_to_backbone,
)
from anvaya.tip_matching import match_competing_paths, match_tip_to_backbone


def _load_experiment(filename: str, name: str):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


previous = _load_experiment("06_tip_matching_validation.py", "tip_validation")
linkage = _load_experiment("08_molecule_linkage_validation.py", "link_validation")


@dataclass(slots=True, frozen=True)
class Condition:
    truth: str
    scenario: str
    reads: list
    quality: int
    truth_codes: set[int] | None = None
    major_reference: str | None = None
    rare_reference: str | None = None


@dataclass(slots=True, frozen=True)
class GraphEventRecord:
    phase: str
    seed: int
    scenario: str
    truth: str
    event_type: str
    event_index: int
    likelihood_status: str
    observations: int
    alternative_observations: int
    reference_observations: int
    ranked_explanation: str
    decision: str
    decision_reasons: str


@dataclass(slots=True, frozen=True)
class GraphRunRecord:
    phase: str
    seed: int
    scenario: str
    truth: str
    tips: int
    incomplete_branches: int
    bubbles: int
    retained_truth_edges: int
    truth_events: int
    matched_truth_events: int
    scored_truth_events: int


@dataclass(slots=True, frozen=True)
class GraphFunnelRecord:
    phase: str
    truth: str
    retained_truth_edges: int
    detected_truth_events: int
    matched_truth_events: int
    scored_truth_events: int


def _quality_for_error_rate(rate: float) -> int:
    return max(2, round(-10.0 * math.log10(rate)))


def _damage_compatible_rare_reference(reference: str, snps: int = 12) -> str:
    """Create separated genuine C>T/G>A variants that challenge damage models."""
    margin = 2 * previous.READ_LENGTH
    candidates = [
        index
        for index in range(margin, len(reference) - margin)
        if reference[index] in "CG"
    ]
    if len(candidates) < snps:
        raise RuntimeError("reference lacks damage-compatible rare-strain loci")
    selected = [
        candidates[(index + 1) * len(candidates) // (snps + 1)]
        for index in range(snps)
    ]
    changed = reference
    for position in selected:
        alternative = "T" if changed[position] == "C" else "A"
        changed = previous.introduce_snp(changed, position, alternative)
    return changed


def _sequence_kmer_codes(sequence: str) -> set[int]:
    return {
        previous.canonical_code(sequence[start : start + previous.K])
        for start in range(len(sequence) - previous.K + 1)
    }


def _conditions(seed: int, reference_length: int, total_coverage: int):
    reference = previous.random_reference(seed)[:reference_length]
    source = previous.simulate_fragments(
        reference,
        previous.READ_LENGTH,
        total_coverage,
        seed + 10_000,
        reverse_fraction=0.5,
    )
    conditions = []
    for offset, (name, profile) in enumerate(previous.DAMAGE_PROFILES.items()):
        reads = previous.add_terminal_damage(
            source,
            profile,
            profile,
            seed + 20_000 + offset,
        )
        conditions.append(
            Condition(
                "damage",
                f"damage_{name}",
                reads,
                35,
                previous.mutated_kmer_codes(source, reads, "damage_positions"),
            )
        )

    for offset, rate in enumerate((0.005, 0.01)):
        reads = previous.add_sequencing_errors(source, rate, seed + 30_000 + offset)
        conditions.append(
            Condition(
                "sequencing_error",
                f"independent_error_{rate:.3f}",
                reads,
                _quality_for_error_rate(rate),
                previous.mutated_kmer_codes(source, reads, "error_positions"),
            )
        )
    systematic = linkage.add_paired_terminal_errors(
        source, previous.TERMINAL_ERROR_RATE, seed + 40_000
    )
    conditions.append(
        Condition(
            "sequencing_error",
            "systematic_terminal_error",
            systematic,
            35,
            previous.mutated_kmer_codes(source, systematic, "error_positions"),
        )
    )

    related = _damage_compatible_rare_reference(reference)
    rare_truth_codes = _sequence_kmer_codes(related) - _sequence_kmer_codes(reference)
    for rare_coverage in (1, 2, 5):
        major_reads = previous.simulate_fragments(
            reference,
            previous.READ_LENGTH,
            total_coverage,
            seed + 60_000 + rare_coverage,
            reverse_fraction=0.5,
        )
        rare_reads = previous.simulate_fragments(
            related,
            previous.READ_LENGTH,
            rare_coverage,
            seed + 70_000 + rare_coverage,
            reverse_fraction=0.5,
        )
        conditions.append(
            Condition(
                "variation",
                f"rare_{rare_coverage}x",
                major_reads + rare_reads,
                35,
                rare_truth_codes,
                major_reference=reference,
                rare_reference=related,
            )
        )
    return conditions


def _build_graph(reads, quality: int):
    return build_bidirected_dbg(
        [read.sequence for read in reads],
        previous.K,
        previous.MIN_COUNT,
        previous.END_WINDOW,
        track_molecule_links=True,
        read_qualities=[(quality,) * len(read.sequence) for read in reads],
    )


def _models_for_seed(conditions: list[Condition]):
    profile_condition = next(
        condition for condition in conditions if condition.scenario == "damage_high"
    )
    graph = _build_graph(profile_condition.reads, profile_condition.quality)
    matches = [
        match_tip_to_backbone(graph, tip)
        for tip in find_weak_tip_candidates(graph)
    ]
    profile = infer_damage_profile(
        graph, [match for match in matches if match is not None]
    )
    return fit_cross_fitted_damage_models(profile)


def _truth_label(condition: Condition, edge_ids, sequence, retained_truth):
    if condition.truth == "variation":
        assert condition.major_reference is not None
        assert condition.rare_reference is not None
        return (
            "variation"
            if previous.strain_label(
                sequence,
                condition.major_reference,
                condition.rare_reference,
            )
            == "rare"
            else None
        )
    return condition.truth if set(edge_ids) & retained_truth else None


def _evaluate_condition(phase, seed, condition, models):
    graph = _build_graph(condition.reads, condition.quality)
    edge_codes = previous.graph_edge_codes(graph)
    retained_truth = {
        edge_codes[code]
        for code in condition.truth_codes or ()
        if code in edge_codes
    }
    tips = find_weak_tip_candidates(graph)
    branches = find_incomplete_branch_candidates(graph)
    bubbles = find_simple_bubbles(graph)
    events = []
    detected_truth_events = 0
    matched_truth_events = 0
    for event_type, candidates in (("tip", tips), ("incomplete_branch", branches)):
        for index, candidate in enumerate(candidates, start=1):
            sequence = spell_edge_path(graph, candidate.start, candidate.edge_ids)
            truth = _truth_label(
                condition,
                candidate.edge_ids,
                sequence,
                retained_truth,
            )
            if truth is None:
                continue
            detected_truth_events += 1
            match = (
                match_tip_to_backbone(graph, candidate)
                if event_type == "tip"
                else match_incomplete_branch_to_backbone(graph, candidate)
            )
            if match is None:
                continue
            matched_truth_events += 1
            likelihood = score_matched_event(graph, match, models)
            events.append((event_type, index, truth, likelihood))
    bubble_path_index = 0
    for bubble in bubbles:
        reference_index = max(
            range(len(bubble.paths)),
            key=lambda index: (
                bubble.paths[index].minimum_edge_support,
                bubble.paths[index].observations,
                -index,
            ),
        )
        reference_path = bubble.paths[reference_index]
        for path_index, path in enumerate(bubble.paths):
            bubble_path_index += 1
            if path_index == reference_index:
                continue
            sequence = spell_edge_path(graph, bubble.start, path.edge_ids)
            truth = _truth_label(
                condition,
                path.edge_ids,
                sequence,
                retained_truth,
            )
            if truth is None:
                continue
            detected_truth_events += 1
            match = match_competing_paths(
                graph,
                bubble.start,
                path.edge_ids,
                reference_path.edge_ids,
            )
            if match is None:
                continue
            matched_truth_events += 1
            likelihood = score_matched_event(graph, match, models)
            events.append(("bubble_path", bubble_path_index, truth, likelihood))
    run = GraphRunRecord(
        phase,
        seed,
        condition.scenario,
        condition.truth,
        len(tips),
        len(branches),
        len(bubbles),
        len(retained_truth),
        detected_truth_events,
        matched_truth_events,
        sum(likelihood.status == "scored" for *_, likelihood in events),
    )
    return run, events


def _write(path: Path, records, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration-seeds", type=int, default=5)
    parser.add_argument("--validation-seeds", type=int, default=2)
    parser.add_argument("--minimum-per-class", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=0.1)
    parser.add_argument("--reference-length", type=int, default=20_000)
    parser.add_argument("--coverage", type=int, default=20)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "experiments/validation/generated/graph_event_calibration_validation"
        ),
    )
    arguments = parser.parse_args()
    if arguments.calibration_seeds < 1 or arguments.validation_seeds < 1:
        parser.error("seed counts must be positive")

    runs = []
    calibration_examples = []
    validation_events = []
    phases = (
        ("calibration", range(2100, 2100 + arguments.calibration_seeds)),
        ("validation", range(3100, 3100 + arguments.validation_seeds)),
    )
    for phase, seeds in phases:
        for seed in seeds:
            conditions = _conditions(seed, arguments.reference_length, arguments.coverage)
            models = _models_for_seed(conditions)
            for condition in conditions:
                run, events = _evaluate_condition(
                    phase, seed, condition, models
                )
                runs.append(run)
                for event_type, index, truth, likelihood in events:
                    scores = scores_from_likelihood(likelihood)
                    if scores is None:
                        continue
                    if phase == "calibration":
                        calibration_examples.append(
                            CalibrationExample(
                                truth,
                                scores,
                                f"{seed}:{condition.scenario}",
                            )
                        )
                    else:
                        validation_events.append(
                            (seed, condition.scenario, event_type, index, truth, likelihood, scores)
                        )

    counts = Counter(example.truth for example in calibration_examples)
    missing = [
        label
        for label in ("damage", "sequencing_error", "variation")
        if counts[label] < arguments.minimum_per_class
    ]
    if missing:
        raise RuntimeError(
            "insufficient scored graph events for "
            + ", ".join(missing)
            + f"; counts={dict(counts)}"
        )
    model = fit_conformal_calibration(
        calibration_examples,
        minimum_per_class=arguments.minimum_per_class,
    )
    confidence = calibrate_score_batch(
        [event[-1] for event in validation_events],
        model,
        alpha=arguments.alpha,
    )
    records = []
    for event, decision in zip(validation_events, confidence, strict=True):
        seed, scenario, event_type, index, truth, likelihood, _ = event
        records.append(
            GraphEventRecord(
                "validation",
                seed,
                scenario,
                truth,
                event_type,
                index,
                likelihood.status,
                likelihood.observations,
                likelihood.alternative_observations,
                likelihood.reference_observations,
                likelihood.best_explanation or "",
                decision.decision,
                ";".join(decision.reasons),
            )
        )

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    _write(
        arguments.output_dir / "per_run.tsv",
        runs,
        tuple(GraphRunRecord.__dataclass_fields__),
    )
    _write(
        arguments.output_dir / "validation_events.tsv",
        records,
        tuple(GraphEventRecord.__dataclass_fields__),
    )
    funnel_records = []
    for phase in ("calibration", "validation"):
        for truth in ("damage", "sequencing_error", "variation"):
            selected = [
                run
                for run in runs
                if run.phase == phase and run.truth == truth
            ]
            funnel_records.append(
                GraphFunnelRecord(
                    phase,
                    truth,
                    sum(run.retained_truth_edges for run in selected),
                    sum(run.truth_events for run in selected),
                    sum(run.matched_truth_events for run in selected),
                    sum(run.scored_truth_events for run in selected),
                )
            )
    _write(
        arguments.output_dir / "evidence_funnel.tsv",
        funnel_records,
        tuple(GraphFunnelRecord.__dataclass_fields__),
    )
    with (arguments.output_dir / "summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("truth", "decision", "events"))
        summary = Counter((record.truth, record.decision) for record in records)
        for key, count in sorted(summary.items()):
            writer.writerow((*key, count))
    print(f"calibration_counts={dict(counts)}")
    print(f"validation_events={len(records)}")
    for record in funnel_records:
        print(
            "evidence_funnel="
            f"{record.phase}:{record.truth}:"
            f"{record.retained_truth_edges}:"
            f"{record.detected_truth_events}:"
            f"{record.matched_truth_events}:"
            f"{record.scored_truth_events}"
        )
    print(f"output={arguments.output_dir}")


if __name__ == "__main__":
    main()
