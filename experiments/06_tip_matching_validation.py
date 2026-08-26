"""Validate weak-tip matching across damage, error, and strain scenarios."""

import argparse
import csv
import random
import statistics
from collections import Counter
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from anvaya.bidirected import (
    _positioned_encoded_kmers,
    build_bidirected_dbg,
    spell_edge_path,
)
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.classification import format_substitutions
from anvaya.sequences import reverse_complement
from anvaya.tip_matching import match_tip_to_backbone
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
SEEDS = tuple(range(1501, 1521))
ERROR_RATES = (0.001, 0.005, 0.010)
RARE_COVERAGES = (1, 2, 5, 10)
DAMAGE_PROFILES = {
    "low": (0.15, 0.10, 0.06, 0.03, 0.01),
    "standard": (0.45, 0.30, 0.20, 0.12, 0.06),
    "high": (0.70, 0.50, 0.35, 0.20, 0.10),
}
TERMINAL_ERROR_RATE = 0.05


@dataclass(frozen=True)
class CandidateRecord:
    condition: str
    scenario: str
    seed: int
    tip_index: int
    label: str
    matched: bool
    ry_exact: bool
    damage_compatible: bool
    relative_coverage: float | None
    sequence_identity: float | None
    ry_identity: float | None
    substitutions: str


@dataclass(frozen=True)
class RunRecord:
    condition: str
    scenario: str
    seed: int
    reads: int
    introduced_events: int
    tips: int
    truth_tips: int
    matched_tips: int
    matched_truth_tips: int
    ry_exact_matches: int
    damage_compatible_matches: int
    damage_compatible_truth_matches: int
    rare_tips: int
    rare_matched_tips: int
    rare_damage_compatible_tips: int


RUN_FIELDS = tuple(RunRecord.__dataclass_fields__)
CANDIDATE_FIELDS = tuple(CandidateRecord.__dataclass_fields__)


def random_reference(seed: int) -> str:
    """Generate one deterministic unambiguous reference."""
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(REFERENCE_LENGTH))


def rare_reference(reference: str, seed: int) -> tuple[str, tuple[int, ...]]:
    """Introduce reproducible, evenly separated SNPs into a related strain."""
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


def add_terminal_errors(
    reads: list[SimulatedRead],
    error_rate: float,
    seed: int,
) -> list[SimulatedRead]:
    """Apply ordinary substitutions only within the terminal evidence window."""
    if not 0.0 <= error_rate <= 1.0:
        raise ValueError("error_rate must be between 0 and 1")
    generator = random.Random(seed)
    changed: list[SimulatedRead] = []
    for read in reads:
        bases = list(read.sequence)
        positions = set(read.error_positions)
        terminal_positions = set(range(min(END_WINDOW, len(bases))))
        terminal_positions.update(
            range(max(0, len(bases) - END_WINDOW), len(bases))
        )
        for position in sorted(terminal_positions):
            if generator.random() >= error_rate:
                continue
            base = bases[position]
            bases[position] = generator.choice(
                [candidate for candidate in "ACGT" if candidate != base]
            )
            positions.add(position)
        changed.append(
            replace(
                read,
                sequence="".join(bases),
                error_positions=tuple(sorted(positions)),
            )
        )
    return changed


def canonical_code(sequence: str) -> int:
    return next(_positioned_encoded_kmers(sequence, len(sequence)))[1]


def mutated_kmer_codes(
    source_reads: list[SimulatedRead],
    changed_reads: list[SimulatedRead],
    position_attribute: str,
) -> set[int]:
    """Return changed k-mers that overlap introduced events."""
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


def graph_edge_codes(graph) -> dict[int, int]:
    return {
        canonical_code(
            spell_edge_path(graph, graph.edge_sources[edge_id], (edge_id,))
        ): edge_id
        for edge_id in range(graph.edge_count)
    }


def contains_sequence(reference: str, sequence: str) -> bool:
    return sequence in reference or reverse_complement(sequence) in reference


def strain_label(
    sequence: str,
    major_reference: str,
    rare_reference_sequence: str,
) -> str:
    """Label a tip by exact membership in either related-strain reference."""
    in_major = contains_sequence(major_reference, sequence)
    in_rare = contains_sequence(rare_reference_sequence, sequence)
    if in_major and in_rare:
        return "shared"
    if in_major:
        return "major"
    if in_rare:
        return "rare"
    return "artifact"


def evaluate(
    condition: str,
    scenario: str,
    seed: int,
    reads: list[SimulatedRead],
    introduced_events: int,
    truth_codes: set[int] | None = None,
    major_reference: str | None = None,
    rare_reference_sequence: str | None = None,
) -> tuple[RunRecord, list[CandidateRecord]]:
    """Evaluate all weak tips without changing the graph."""
    graph = build_bidirected_dbg(
        [read.sequence for read in reads],
        K,
        MIN_COUNT,
        END_WINDOW,
    )
    tips = find_weak_tip_candidates(graph)
    retained_truth: set[int] = set()
    if truth_codes is not None:
        edge_codes = graph_edge_codes(graph)
        retained_truth = {
            edge_codes[code] for code in truth_codes if code in edge_codes
        }

    candidates: list[CandidateRecord] = []
    for tip_index, tip in enumerate(tips, start=1):
        match = match_tip_to_backbone(graph, tip)
        sequence = (
            match.tip_sequence
            if match is not None
            else spell_edge_path(graph, tip.start, tip.edge_ids)
        )
        if major_reference is not None and rare_reference_sequence is not None:
            label = strain_label(sequence, major_reference, rare_reference_sequence)
        elif retained_truth and set(tip.edge_ids) & retained_truth:
            label = condition
        elif condition == "clean":
            label = "clean"
        else:
            label = "artifact"
        candidates.append(
            CandidateRecord(
                condition=condition,
                scenario=scenario,
                seed=seed,
                tip_index=tip_index,
                label=label,
                matched=match is not None,
                ry_exact=match is not None and match.ry_identity == 1.0,
                damage_compatible=(
                    match is not None and match.damage_compatible
                ),
                relative_coverage=(
                    None if match is None else match.relative_coverage
                ),
                sequence_identity=(
                    None if match is None else match.sequence_identity
                ),
                ry_identity=None if match is None else match.ry_identity,
                substitutions=(
                    "" if match is None else format_substitutions(match.substitutions)
                ),
            )
        )

    counts = Counter(candidate.label for candidate in candidates)
    truth_labels = {"damage", "error"}
    truth_candidates = [
        candidate for candidate in candidates if candidate.label in truth_labels
    ]
    return RunRecord(
        condition=condition,
        scenario=scenario,
        seed=seed,
        reads=len(reads),
        introduced_events=introduced_events,
        tips=len(candidates),
        truth_tips=len(truth_candidates),
        matched_tips=sum(candidate.matched for candidate in candidates),
        matched_truth_tips=sum(candidate.matched for candidate in truth_candidates),
        ry_exact_matches=sum(candidate.ry_exact for candidate in candidates),
        damage_compatible_matches=sum(
            candidate.damage_compatible for candidate in candidates
        ),
        damage_compatible_truth_matches=sum(
            candidate.damage_compatible for candidate in truth_candidates
        ),
        rare_tips=counts["rare"],
        rare_matched_tips=sum(
            candidate.matched and candidate.label == "rare"
            for candidate in candidates
        ),
        rare_damage_compatible_tips=sum(
            candidate.damage_compatible and candidate.label == "rare"
            for candidate in candidates
        ),
    ), candidates


def event_count(reads: list[SimulatedRead], attribute: str) -> int:
    return sum(len(getattr(read, attribute)) for read in reads)


def run_matrix(
    seeds: tuple[int, ...],
) -> tuple[list[RunRecord], list[CandidateRecord]]:
    """Run matched bidirectional damage, error, and strain conditions."""
    runs: list[RunRecord] = []
    candidates: list[CandidateRecord] = []
    zero_profile = (0.0,) * END_WINDOW

    for seed in seeds:
        reference = random_reference(seed)
        source = simulate_fragments(
            reference,
            READ_LENGTH,
            TOTAL_COVERAGE,
            seed + 10_000,
            reverse_fraction=0.5,
        )
        conditions: list[
            tuple[str, str, list[SimulatedRead], str | None, set[int] | None]
        ] = [("clean", "clean", source, None, None)]

        for profile_name, profile in DAMAGE_PROFILES.items():
            for end_name, five_prime, three_prime in (
                ("5p", profile, zero_profile),
                ("3p", zero_profile, profile),
                ("both", profile, profile),
            ):
                damaged = add_terminal_damage(
                    source,
                    five_prime,
                    three_prime,
                    seed + 20_000 + len(conditions),
                )
                conditions.append(
                    (
                        "damage",
                        f"damage_{profile_name}_{end_name}",
                        damaged,
                        "damage_positions",
                        mutated_kmer_codes(
                            source, damaged, "damage_positions"
                        ),
                    )
                )

        for error_rate in ERROR_RATES:
            errors = add_sequencing_errors(
                source,
                error_rate,
                seed + 30_000 + int(error_rate * 1000),
            )
            conditions.append(
                (
                    "error",
                    f"error_{error_rate:.3f}",
                    errors,
                    "error_positions",
                    mutated_kmer_codes(source, errors, "error_positions"),
                )
            )

        terminal_errors = add_terminal_errors(
            source,
            TERMINAL_ERROR_RATE,
            seed + 40_000,
        )
        conditions.append(
            (
                "error",
                "error_terminal",
                terminal_errors,
                "error_positions",
                mutated_kmer_codes(
                    source, terminal_errors, "error_positions"
                ),
            )
        )

        for condition, scenario, reads, attribute, truth in conditions:
            run, found = evaluate(
                condition,
                scenario,
                seed,
                reads,
                event_count(reads, attribute) if attribute else 0,
                truth_codes=truth,
            )
            runs.append(run)
            candidates.extend(found)

        related, _ = rare_reference(reference, seed + 50_000)
        for rare_coverage in RARE_COVERAGES:
            major_reads = simulate_fragments(
                reference,
                READ_LENGTH,
                TOTAL_COVERAGE - rare_coverage,
                seed + 60_000 + rare_coverage,
                reverse_fraction=0.5,
            )
            rare_reads = simulate_fragments(
                related,
                READ_LENGTH,
                rare_coverage,
                seed + 70_000 + rare_coverage,
                reverse_fraction=0.5,
            )
            run, found = evaluate(
                "rare",
                f"rare_{rare_coverage}x",
                seed,
                major_reads + rare_reads,
                SNP_COUNT,
                major_reference=reference,
                rare_reference_sequence=related,
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
    parameters = {
        "reference_length": REFERENCE_LENGTH,
        "read_length": READ_LENGTH,
        "total_coverage": TOTAL_COVERAGE,
        "k": K,
        "min_count": MIN_COUNT,
        "end_window": END_WINDOW,
        "seeds": ",".join(map(str, seeds)),
        "error_rates": ",".join(map(str, ERROR_RATES)),
        "terminal_error_rate": TERMINAL_ERROR_RATE,
        "rare_coverages": ",".join(map(str, RARE_COVERAGES)),
    }
    parameters.update(
        {
            f"damage_{name}": ",".join(map(str, profile))
            for name, profile in DAMAGE_PROFILES.items()
        }
    )
    with (output_dir / "parameters.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("parameter", "value"))
        writer.writerows(parameters.items())

    _write_dataclasses(output_dir / "per_run.tsv", runs, RUN_FIELDS)
    _write_dataclasses(output_dir / "candidates.tsv", candidates, CANDIDATE_FIELDS)

    with (output_dir / "summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("condition", "scenario", "metric", "mean", "stdev", "minimum", "maximum", "n")
        )
        scenarios = sorted({(run.condition, run.scenario) for run in runs})
        for condition, scenario in scenarios:
            selected = [
                run
                for run in runs
                if run.condition == condition and run.scenario == scenario
            ]
            for metric in RUN_FIELDS[3:]:
                values = [getattr(run, metric) for run in selected]
                writer.writerow(
                    (
                        condition,
                        scenario,
                        metric,
                        statistics.mean(values),
                        statistics.stdev(values) if len(values) > 1 else 0.0,
                        min(values),
                        max(values),
                        len(values),
                    )
                )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("experiments/validation/generated/tip_matching_matrix"),
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
