"""Validate read threading with short, variable, damaged ancient-like reads."""

import argparse
import csv
import importlib.util
import math
import random
import sys
from dataclasses import asdict
from pathlib import Path

from anvaya.sequences import reverse_complement
from helpers.simulation import SimulatedRead, add_sequencing_errors, add_terminal_damage


def _load_read_threading_experiment():
    path = Path(__file__).with_name("16_read_threading_validation.py")
    spec = importlib.util.spec_from_file_location("read_threading_validation", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


threading = _load_read_threading_experiment()
previous = threading.previous

MINIMUM_LENGTH = 35
MAXIMUM_LENGTH = 120
MEDIAN_LENGTH = 55
LOG_LENGTH_SIGMA = 0.35


def simulate_ancient_fragments(
    reference: str,
    coverage: float,
    seed: int,
) -> list[SimulatedRead]:
    """Sample a truncated log-normal fragment distribution at target coverage."""
    generator = random.Random(seed)
    target_bases = math.ceil(coverage * len(reference))
    sampled_bases = 0
    reads: list[SimulatedRead] = []
    while sampled_bases < target_bases:
        length = round(generator.lognormvariate(math.log(MEDIAN_LENGTH), LOG_LENGTH_SIGMA))
        length = min(MAXIMUM_LENGTH, max(MINIMUM_LENGTH, length))
        start = generator.randrange(len(reference) - length + 1)
        end = start + length
        sequence = reference[start:end]
        is_reverse = generator.random() < 0.5
        if is_reverse:
            sequence = reverse_complement(sequence)
        reads.append(
            SimulatedRead(
                name=f"ancient-{len(reads)}",
                sequence=sequence,
                source_sequence=sequence,
                reference_start=start,
                reference_end=end,
                is_reverse=is_reverse,
            )
        )
        sampled_bases += length
    return reads


def run_matrix(seeds: tuple[int, ...]):
    """Run short-fragment damage, error, strain, and contamination conditions."""
    records = []
    profile = previous.DAMAGE_PROFILES["standard"]
    for seed in seeds:
        reference = previous.random_reference(seed)
        source = simulate_ancient_fragments(reference, 20, seed + 10_000)
        damaged = add_terminal_damage(source, profile, profile, seed + 20_000)
        damaged_error = add_sequencing_errors(damaged, 0.005, seed + 30_000)
        rare_reference, _ = previous.rare_reference(reference, seed + 40_000)
        contaminant = previous.random_reference(seed + 1_000_000)
        mixed = (
            simulate_ancient_fragments(reference, 15, seed + 50_000)
            + simulate_ancient_fragments(rare_reference, 5, seed + 60_000)
            + simulate_ancient_fragments(contaminant, 5, seed + 70_000)
        )
        mixed = add_terminal_damage(mixed, profile, profile, seed + 80_000)
        for condition, reads, truth in (
            ("ancient_short", source, (reference,)),
            ("ancient_damage", damaged, (reference,)),
            ("ancient_damage_error", damaged_error, (reference,)),
            (
                "ancient_strain_contamination",
                mixed,
                (reference, rare_reference, contaminant),
            ),
        ):
            records.append(threading._evaluate(condition, seed, reads, truth))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seeds", nargs="+", type=int, default=range(1701, 1706))
    arguments = parser.parse_args()
    records = run_matrix(tuple(arguments.seeds))
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0].__dataclass_fields__)
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)
    links = sum(record.reciprocal_links for record in records)
    true_links = sum(record.true_links for record in records)
    ry_true_links = sum(record.ry_true_links for record in records)
    guarded_links = sum(record.repeat_guarded_links for record in records)
    print(f"runs={len(records)}")
    print(f"reciprocal_links={links}")
    print(f"true_links={true_links}")
    print(f"false_links={links - true_links}")
    print(f"precision={true_links / links if links else 1.0:.6f}")
    print(f"ry_true_links={ry_true_links}")
    print(f"ry_false_links={links - ry_true_links}")
    print(f"ry_precision={ry_true_links / links if links else 1.0:.6f}")
    print(f"repeat_guarded_links={guarded_links}")
    print(
        "false_contig_delta="
        f"{sum(record.false_extended_contigs - record.false_baseline_contigs for record in records)}"
    )
    print(f"output={arguments.output}")


if __name__ == "__main__":
    main()
