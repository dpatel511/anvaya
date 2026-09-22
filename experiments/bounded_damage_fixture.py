"""Generate deterministic bounded FASTQ inputs for operational profiling."""
import argparse
import json
import random
from pathlib import Path

from anvaya.sequences import reverse_complement


def write_profile(prefix, damage):
    header = [f"{source}>{target}" for source in "ACGT" for target in "ACGT"
              if source != target]
    for end, transition in (("5p", "C>T"), ("3p", "G>A")):
        with Path(f"{prefix}{end}.prof").open("w", encoding="utf-8") as handle:
            handle.write("\t".join(header) + "\n")
            for position in range(5):
                rate = damage * .65 ** position
                handle.write("\t".join(
                    str(rate if label == transition else 0.0) for label in header
                ) + "\n")


def write_fastq(path, count, seed, damage=.25, error_rate=.002):
    reference_rng = random.Random(seed)
    genome_length = max(2400, round(count * 115 / 10))
    genome = "".join(reference_rng.choice("ACGT") for _ in range(genome_length))
    placement_rng = random.Random(seed + 1)
    damage_rng = random.Random(seed + 2)
    error_rng = random.Random(seed + 3)
    bases = 0
    with path.open("w", encoding="utf-8") as handle:
        for index in range(count):
            length = placement_rng.randint(31, 200)
            start = placement_rng.randrange(genome_length - length + 1)
            sequence = genome[start:start + length]
            if placement_rng.randrange(2):
                sequence = reverse_complement(sequence)
            observed = list(sequence)
            qualities = [35] * length
            for position, base in enumerate(sequence):
                distance = position if base == "C" else length - 1 - position
                if (
                    base in "CG" and distance < 5
                    and damage_rng.random() < damage * .65 ** distance
                ):
                    observed[position] = "T" if base == "C" else "A"
            for position, base in enumerate(observed):
                if error_rng.random() < error_rate:
                    observed[position] = error_rng.choice(
                        [candidate for candidate in "ACGT" if candidate != base]
                    )
                    qualities[position] = 10
            sequence = "".join(observed)
            bases += len(sequence)
            quality = "".join(chr(value + 33) for value in qualities)
            handle.write(f"@profile_{count}_{index}\n{sequence}\n+\n{quality}\n")
    return {"reads": count, "bases": bases, "genome_length": genome_length}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reads", nargs="+", type=int, default=(250, 500, 1000))
    parser.add_argument("--seed", type=int, default=198001)
    arguments = parser.parse_args()
    if arguments.output.exists():
        parser.error("output already exists")
    if any(count < 1 for count in arguments.reads):
        parser.error("read counts must be positive")
    arguments.output.mkdir(parents=True)
    prefix = arguments.output / "profile"
    write_profile(prefix, .25)
    fixtures = []
    for offset, count in enumerate(arguments.reads):
        path = arguments.output / f"reads-{count}.fastq"
        fixtures.append({
            **write_fastq(path, count, arguments.seed + offset * 10),
            "path": str(path),
        })
    (arguments.output / "manifest.json").write_text(json.dumps({
        "purpose": "operational profiling only; not biological validation",
        "seed": arguments.seed,
        "damage": .25,
        "sequencing_error_rate": .002,
        "read_length_range": [31, 200],
        "fixtures": fixtures,
        "resource_stop_budget": {
            "wall_seconds_at_1000_reads": 120,
            "peak_rss_kib_at_1000_reads": 2097152
        }
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
