"""Fit deterministic held-out damage folds from an existing profile JSON."""

import argparse
import csv
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.damage_likelihood import fit_candidate_damage_model


@dataclass(frozen=True)
class FoldRecord:
    fold: int
    held_out_loci: int
    training_loci: int
    status: str
    reasons: str
    fitted_probabilities: str
    likelihood_ratio_statistic: float | None
    elapsed_seconds: float


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/validation/generated/"
            "cross_fit_profile_validation/folds.tsv"
        ),
    )
    arguments = parser.parse_args()
    if arguments.folds < 2:
        raise ValueError("folds must be at least 2")

    payload = json.loads(arguments.profile.read_text(encoding="utf-8"))
    end_window = int(payload["end_window"])
    loci = sorted(
        payload["loci"],
        key=lambda locus: (
            locus["alternative_edge_id"],
            locus["reference_edge_id"],
            locus["channel"],
        ),
    )
    records = []
    for fold in range(arguments.folds):
        training = [
            (tuple(locus["alternative"]), tuple(locus["reference"]))
            for index, locus in enumerate(loci)
            if index % arguments.folds != fold
        ]
        started = time.perf_counter()
        fit = fit_candidate_damage_model(training, end_window)
        records.append(
            FoldRecord(
                fold=fold + 1,
                held_out_loci=sum(
                    index % arguments.folds == fold
                    for index in range(len(loci))
                ),
                training_loci=len(training),
                status=fit.status,
                reasons=";".join(fit.reasons),
                fitted_probabilities=";".join(
                    f"{value:.6f}" for value in fit.fitted_probabilities
                ),
                likelihood_ratio_statistic=fit.likelihood_ratio_statistic,
                elapsed_seconds=time.perf_counter() - started,
            )
        )

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            tuple(FoldRecord.__dataclass_fields__),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)
    print(f"folds={len(records)}")
    print(f"fitted={sum(record.status == 'fitted' for record in records)}")
    print(f"output={arguments.output}")


if __name__ == "__main__":
    main()
