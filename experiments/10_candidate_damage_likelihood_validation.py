"""Validate candidate damage likelihoods against non-damage controls."""

import argparse
import csv
import importlib.util
import statistics
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


def _load_experiment(filename: str, module_name: str):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


profiles = _load_experiment(
    "09_damage_profile_validation.py",
    "damage_profile_validation",
)
previous = profiles.previous
linkage = _load_experiment(
    "08_molecule_linkage_validation.py",
    "molecule_linkage_validation",
)


@dataclass(frozen=True)
class ControlRecord:
    condition: str
    seed: int
    observed_loci: int
    fit_status: str
    likelihood_ratio_statistic: float | None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "experiments/validation/generated/"
            "candidate_damage_likelihood_validation"
        ),
    )
    parser.add_argument(
        "--seed-count",
        type=int,
        default=len(previous.SEEDS),
        choices=range(1, len(previous.SEEDS) + 1),
    )
    arguments = parser.parse_args()
    zero = (0.0,) * previous.END_WINDOW
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
        conditions = (
            ("clean", source),
            (
                "terminal_error",
                previous.add_terminal_errors(
                    source,
                    previous.TERMINAL_ERROR_RATE,
                    seed + 40_000,
                ),
            ),
            (
                "correlated_terminal_error",
                linkage.add_paired_terminal_errors(
                    source,
                    0.15,
                    seed + 50_000,
                ),
            ),
        )
        for condition, reads in conditions:
            result = profiles.evaluate(condition, seed, reads, zero)
            records.append(
                ControlRecord(
                    condition=condition,
                    seed=seed,
                    observed_loci=result.eligible_loci,
                    fit_status=result.fit_status,
                    likelihood_ratio_statistic=(
                        result.likelihood_ratio_statistic
                    ),
                )
            )

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    with (arguments.output_dir / "per_run.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            tuple(ControlRecord.__dataclass_fields__),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)
    with (arguments.output_dir / "summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "condition",
                "mean_observed_loci",
                "fitted_runs",
                "mean_likelihood_ratio_statistic",
                "runs",
            )
        )
        for condition in (
            "clean",
            "terminal_error",
            "correlated_terminal_error",
        ):
            selected = [r for r in records if r.condition == condition]
            statistics_values = [
                r.likelihood_ratio_statistic
                for r in selected
                if r.likelihood_ratio_statistic is not None
            ]
            writer.writerow(
                (
                    condition,
                    statistics.mean(r.observed_loci for r in selected),
                    sum(r.fit_status == "fitted" for r in selected),
                    (
                        statistics.mean(statistics_values)
                        if statistics_values
                        else ""
                    ),
                    len(selected),
                )
            )
    print(f"runs={len(records)}")
    print(f"output={arguments.output_dir}")


if __name__ == "__main__":
    main()
