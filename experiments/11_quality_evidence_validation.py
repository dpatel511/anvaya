"""Validate quality propagation and sequencing-error likelihood reporting."""

import argparse
import csv
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.tip_classification import classify_tip_match
from anvaya.tip_matching import match_tip_to_backbone


@dataclass(frozen=True)
class Record:
    condition: str
    alternative_support: int
    quality: int | None
    quality_observations: int
    missing_quality_observations: int
    mean_base_quality: float | None
    sequencing_error_log_likelihood: float | None
    per_observation_log_likelihood: float | None


def evaluate(
    condition: str,
    alternative_support: int,
    quality: int | None,
) -> Record:
    backbone = "AAGCCCAAA"
    alternative = "AAGCCTAAA"
    backbone_support = 30
    reads = [backbone] * backbone_support + [alternative] * alternative_support
    read_qualities = (
        None
        if quality is None
        else [(35,) * len(backbone)] * backbone_support
        + [(quality,) * len(alternative)] * alternative_support
    )
    graph = build_bidirected_dbg(
        reads,
        5,
        end_window=11,
        track_molecule_links=True,
        read_qualities=read_qualities,
    )
    tip = find_weak_tip_candidates(graph)[0]
    match = match_tip_to_backbone(graph, tip)
    if match is None:
        raise RuntimeError("controlled weak tip did not match its backbone")
    evidence = classify_tip_match(graph, match).evidence
    likelihood = evidence.sequencing_error_log_likelihood
    return Record(
        condition=condition,
        alternative_support=alternative_support,
        quality=quality,
        quality_observations=evidence.quality_observations,
        missing_quality_observations=evidence.missing_quality_observations,
        mean_base_quality=evidence.mean_base_quality,
        sequencing_error_log_likelihood=likelihood,
        per_observation_log_likelihood=(
            likelihood / evidence.quality_observations
            if likelihood is not None and evidence.quality_observations
            else None
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "experiments/validation/generated/quality_evidence_validation"
        ),
    )
    arguments = parser.parse_args()
    conditions = (
        ("low_quality_error", 8),
        ("moderate_quality", 20),
        ("high_quality_damage_or_variant", 35),
        ("missing_quality", None),
    )
    records = [
        evaluate(condition, support, quality)
        for condition, quality in conditions
        for support in range(1, 6)
    ]

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    output = arguments.output_dir / "matrix.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            tuple(Record.__dataclass_fields__),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)

    for condition, quality in conditions:
        selected = [record for record in records if record.condition == condition]
        if quality is None:
            if any(record.sequencing_error_log_likelihood is not None for record in selected):
                raise AssertionError("missing qualities produced an error likelihood")
            continue
        expected = math.log((10.0 ** (-quality / 10.0)) / 3.0)
        if any(
            not math.isclose(record.per_observation_log_likelihood or 0.0, expected)
            for record in selected
        ):
            raise AssertionError("reported likelihood does not match Phred model")

    print(f"runs={len(records)}")
    print(f"output={output}")


if __name__ == "__main__":
    main()
