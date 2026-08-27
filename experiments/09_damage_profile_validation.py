"""Validate inferred candidate-locus profiles against simulated damage shape."""

import argparse
import csv
import importlib.util
import math
import statistics
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.damage_profile import infer_damage_profile
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
class ProfileRecord:
    condition: str
    seed: int
    matched_paths: int
    eligible_loci: int
    damage_observations: int
    damage_correlation: float | None
    damage_decreasing: bool | None


def pearson(left, right):
    """Return Pearson correlation, or None when either vector is constant."""
    if len(left) != len(right) or len(left) < 2:
        return None
    left_mean = statistics.mean(left)
    right_mean = statistics.mean(right)
    numerator = sum(
        (a - left_mean) * (b - right_mean)
        for a, b in zip(left, right, strict=True)
    )
    left_scale = math.sqrt(sum((value - left_mean) ** 2 for value in left))
    right_scale = math.sqrt(sum((value - right_mean) ** 2 for value in right))
    if left_scale == 0 or right_scale == 0:
        return None
    return numerator / (left_scale * right_scale)


def _profile_values(profile, alternative_field, reference_field, fraction_field):
    values = []
    observations = 0
    for value in profile.bins:
        alternative = getattr(value, alternative_field)
        reference = getattr(value, reference_field)
        fraction = getattr(value, fraction_field)
        observations += alternative + reference
        values.append(fraction)
    return values, observations


def evaluate(condition, seed, reads, expected):
    graph = build_bidirected_dbg(
        [read.sequence for read in reads],
        previous.K,
        previous.MIN_COUNT,
        previous.END_WINDOW,
        track_molecule_links=True,
    )
    matches = []
    for tip in find_weak_tip_candidates(graph):
        match = match_tip_to_backbone(graph, tip)
        if match is not None:
            matches.append(match)
    profile = infer_damage_profile(graph, matches)
    combined, observations = _profile_values(
        profile,
        "damage_alternative",
        "damage_reference",
        "damage_fraction",
    )

    def compare(values):
        observed = [value for value in values if value is not None]
        expected_observed = [
            expected[index]
            for index, value in enumerate(values)
            if value is not None
        ]
        return (
            pearson(expected_observed, observed),
            observed[0] > observed[-1] if len(observed) >= 2 else None,
        )

    correlation, decreasing = compare(combined)
    return ProfileRecord(
        condition=condition,
        seed=seed,
        matched_paths=profile.matched_paths,
        eligible_loci=profile.eligible_loci,
        damage_observations=observations,
        damage_correlation=correlation,
        damage_decreasing=decreasing,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "experiments/validation/generated/damage_profile_validation"
        ),
    )
    parser.add_argument(
        "--seed-count",
        type=int,
        default=len(previous.SEEDS),
        choices=range(1, len(previous.SEEDS) + 1),
    )
    arguments = parser.parse_args()
    records = []
    for seed in previous.SEEDS[: arguments.seed_count]:
        reference = previous.random_reference(seed)
        source = previous.simulate_fragments(
            reference,
            previous.READ_LENGTH,
            previous.TOTAL_COVERAGE,
            seed + 10_000,
            reverse_fraction=0.5,
        )
        for index, (name, expected) in enumerate(
            previous.DAMAGE_PROFILES.items()
        ):
            damaged = previous.add_terminal_damage(
                source,
                expected,
                expected,
                seed + 80_000 + index,
            )
            records.append(evaluate(name, seed, damaged, expected))

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    fields = tuple(ProfileRecord.__dataclass_fields__)
    with (arguments.output_dir / "per_run.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)
    with (arguments.output_dir / "summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "condition",
                "mean_eligible_loci",
                "mean_damage_correlation",
                "damage_decreasing_runs",
                "runs",
            )
        )
        for condition in previous.DAMAGE_PROFILES:
            selected = [r for r in records if r.condition == condition]
            correlations = [
                r.damage_correlation
                for r in selected
                if r.damage_correlation is not None
            ]
            writer.writerow(
                (
                    condition,
                    statistics.mean(r.eligible_loci for r in selected),
                    statistics.mean(correlations) if correlations else "",
                    sum(r.damage_decreasing is True for r in selected),
                    len(selected),
                )
            )
    print(f"runs={len(records)}")
    print(f"output={arguments.output_dir}")


if __name__ == "__main__":
    main()
