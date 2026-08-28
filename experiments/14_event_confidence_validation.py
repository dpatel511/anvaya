"""Validate conformal event decisions on independently seeded samples."""

import argparse
import csv
import random
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from anvaya.event_calibration import (
    CalibrationExample,
    calibrate_score_batch,
    fit_conformal_calibration,
    scores_from_likelihood,
    write_calibration_model,
)
from anvaya.event_likelihood import NucleotideObservation, compare_event_likelihoods


DAMAGE_CURVE = (0.25, 0.16, 0.10, 0.06, 0.03)
CALIBRATION_SEEDS = tuple(range(1000))
TEST_SEEDS = tuple(range(10000, 10100))


@dataclass(slots=True, frozen=True)
class Condition:
    name: str
    truth: str
    probabilities: tuple[float, ...]
    quality: int


def _simulate(condition: Condition, seed: int) -> list[NucleotideObservation]:
    generator = random.Random(seed)
    observations = []
    for distance, latent in enumerate(condition.probabilities):
        error = 10.0 ** (-condition.quality / 10.0)
        observed = latent * (1.0 - error) + (1.0 - latent) * error / 3.0
        for _ in range(100):
            observations.append(
                NucleotideObservation(
                    generator.random() < observed,
                    distance,
                    condition.quality,
                )
            )
    return observations


def _scores(condition: Condition, seed: int):
    likelihood = compare_event_likelihoods(
        _simulate(condition, seed),
        DAMAGE_CURVE,
    )
    result = scores_from_likelihood(likelihood)
    if result is None:
        raise RuntimeError(f"simulation was not scoreable: {condition.name}/{seed}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("experiments/validation/generated/event_confidence_validation"),
    )
    arguments = parser.parse_args()
    calibration_conditions = (
        Condition("damage", "damage", DAMAGE_CURVE, 35),
        Condition(
            "profile_misspecification",
            "damage",
            (0.20, 0.18, 0.14, 0.09, 0.05),
            35,
        ),
        Condition("independent_error", "sequencing_error", (0.0,) * 5, 12),
        Condition("systematic_error", "sequencing_error", DAMAGE_CURVE, 35),
        Condition("variation", "variation", (0.12,) * 5, 35),
    )
    examples = [
        CalibrationExample(
            condition.truth,
            _scores(condition, seed),
            f"{condition.name}-{seed}",
        )
        for condition in calibration_conditions
        for seed in CALIBRATION_SEEDS
    ]
    model = fit_conformal_calibration(examples, minimum_per_class=1000)

    test_conditions = (
        Condition("damage", "damage", DAMAGE_CURVE, 35),
        Condition("independent_error", "sequencing_error", (0.0,) * 5, 12),
        Condition("systematic_error", "sequencing_error", DAMAGE_CURVE, 35),
        Condition("variation", "variation", (0.12,) * 5, 35),
        Condition("profile_misspecification", "damage", (0.20, 0.18, 0.14, 0.09, 0.05), 35),
    )
    output_rows = []
    for condition in test_conditions:
        scores = [_scores(condition, seed) for seed in TEST_SEEDS]
        confidence = calibrate_score_batch(scores, model, alpha=0.01)
        for seed, result in zip(TEST_SEEDS, confidence, strict=True):
            output_rows.append(
                {
                    "condition": condition.name,
                    "truth": condition.truth,
                    "seed": seed,
                    "decision": result.decision,
                    "prediction_set": ";".join(result.prediction_set),
                    "damage_pvalue": result.damage_pvalue,
                    "error_pvalue": result.error_pvalue,
                    "variation_pvalue": result.variation_pvalue,
                    "damage_qvalue": result.damage_qvalue,
                    "variation_qvalue": result.variation_qvalue,
                }
            )

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    write_calibration_model(model, arguments.output_dir / "calibration-model.json")
    with (arguments.output_dir / "per_run.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            tuple(output_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)
    with (arguments.output_dir / "summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("condition", "truth", "decision", "runs"))
        for condition in test_conditions:
            counts = Counter(
                row["decision"]
                for row in output_rows
                if row["condition"] == condition.name
            )
            for decision, count in sorted(counts.items()):
                writer.writerow((condition.name, condition.truth, decision, count))
    print(f"calibration_examples={len(examples)}")
    print(f"test_runs={len(output_rows)}")
    print(f"output={arguments.output_dir}")


if __name__ == "__main__":
    main()
