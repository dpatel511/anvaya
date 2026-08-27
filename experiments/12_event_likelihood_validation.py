"""Validate damage, sequencing-error, and variation likelihood rankings."""

import argparse
import csv
import random
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.event_likelihood import (
    NucleotideObservation,
    compare_event_likelihoods,
)


SEEDS = tuple(range(20))
DAMAGE_CURVE = (0.25, 0.16, 0.10, 0.06, 0.03)


@dataclass(frozen=True)
class Record:
    condition: str
    seed: int
    expected: str
    predicted: str | None
    observations: int
    alternative_observations: int
    damage_log_likelihood: float | None
    error_log_likelihood: float | None
    variation_penalized_log_likelihood: float | None
    variation_frequency: float | None
    log_likelihood_margin: float | None


def _observed_alternative_probability(latent: float, quality: int) -> float:
    error = 10.0 ** (-quality / 10.0)
    return latent * (1.0 - error) + (1.0 - latent) * error / 3.0


def simulate(
    seed: int,
    latent_probabilities: tuple[float, ...],
    quality: int,
) -> list[NucleotideObservation]:
    generator = random.Random(seed)
    observations = []
    for distance, latent in enumerate(latent_probabilities):
        observed_probability = _observed_alternative_probability(
            latent,
            quality,
        )
        for _ in range(100):
            observations.append(
                NucleotideObservation(
                    alternative=generator.random() < observed_probability,
                    distance=distance,
                    quality=quality,
                )
            )
    return observations


def evaluate(
    condition: str,
    expected: str,
    seed: int,
    generating_probabilities: tuple[float, ...],
    quality: int,
) -> Record:
    observations = simulate(seed, generating_probabilities, quality)
    result = compare_event_likelihoods(observations, DAMAGE_CURVE)
    return Record(
        condition=condition,
        seed=seed,
        expected=expected,
        predicted=result.best_explanation,
        observations=result.observations,
        alternative_observations=result.alternative_observations,
        damage_log_likelihood=result.damage_log_likelihood,
        error_log_likelihood=result.error_log_likelihood,
        variation_penalized_log_likelihood=(
            result.variation_penalized_log_likelihood
        ),
        variation_frequency=result.variation_frequency,
        log_likelihood_margin=result.log_likelihood_margin,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "experiments/validation/generated/event_likelihood_validation"
        ),
    )
    arguments = parser.parse_args()
    conditions = (
        ("damage", "damage", DAMAGE_CURVE, 35),
        ("sequencing_error", "sequencing_error", (0.0,) * 5, 12),
        ("variation", "variation", (0.12,) * 5, 35),
        (
            "correlated_damage_shaped_error",
            "sequencing_error",
            DAMAGE_CURVE,
            35,
        ),
    )
    records = [
        evaluate(condition, expected, seed, probabilities, quality)
        for condition, expected, probabilities, quality in conditions
        for seed in SEEDS
    ]

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    with (arguments.output_dir / "per_run.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            tuple(Record.__dataclass_fields__),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)

    with (arguments.output_dir / "summary.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("condition", "expected", "prediction", "runs"))
        for condition, expected, _, _ in conditions:
            counts = Counter(
                record.predicted
                for record in records
                if record.condition == condition
            )
            for prediction, count in sorted(counts.items()):
                writer.writerow((condition, expected, prediction, count))

    print(f"runs={len(records)}")
    print(f"output={arguments.output_dir}")


if __name__ == "__main__":
    main()
