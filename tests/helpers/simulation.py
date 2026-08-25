"""Controlled read simulation with explicit ground truth."""

import math
import random
from dataclasses import dataclass, replace
from numbers import Real
from collections.abc import Sequence

from anvaya.sequences import normalize_dna


@dataclass(frozen=True)
class SimulatedRead:
    """A simulated read with its source coordinates and introduced events."""

    name: str
    sequence: str
    source_sequence: str
    reference_start: int
    reference_end: int
    damage_positions: tuple[int, ...] = ()
    error_positions: tuple[int, ...] = ()


def _validate_probability(value: Real, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a number")
    probability = float(value)
    if not 0.0 <= probability <= 1.0:
        raise ValueError(f"{name} must be between 0 and 1")
    return probability


def simulate_fragments(
    reference: str,
    read_length: int,
    coverage: Real,
    seed: int,
) -> list[SimulatedRead]:
    """Sample fixed-length forward reads from a reference with replacement."""
    normalized = normalize_dna(reference)
    if "N" in normalized:
        raise ValueError("simulation reference must not contain N")
    if not isinstance(read_length, int) or isinstance(read_length, bool):
        raise TypeError("read_length must be an integer")
    if read_length < 1 or read_length > len(normalized):
        raise ValueError("read_length must be between 1 and the reference length")
    if isinstance(coverage, bool) or not isinstance(coverage, Real):
        raise TypeError("coverage must be a number")
    if coverage <= 0:
        raise ValueError("coverage must be greater than 0")
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise TypeError("seed must be an integer")

    read_count = math.ceil(float(coverage) * len(normalized) / read_length)
    random_generator = random.Random(seed)
    reads: list[SimulatedRead] = []

    for index in range(read_count):
        start = random_generator.randrange(len(normalized) - read_length + 1)
        end = start + read_length
        source = normalized[start:end]
        reads.append(
            SimulatedRead(
                name=f"read-{index}",
                sequence=source,
                source_sequence=source,
                reference_start=start,
                reference_end=end,
            )
        )

    return reads


def add_sequencing_errors(
    reads: Sequence[SimulatedRead],
    error_rate: Real,
    seed: int,
) -> list[SimulatedRead]:
    """Apply independent random base substitutions to simulated reads."""
    probability = _validate_probability(error_rate, "error_rate")
    random_generator = random.Random(seed)
    mutated_reads: list[SimulatedRead] = []

    for read in reads:
        bases = list(read.sequence)
        positions: list[int] = []
        for position, base in enumerate(bases):
            if random_generator.random() < probability:
                bases[position] = random_generator.choice(
                    [candidate for candidate in "ACGT" if candidate != base]
                )
                positions.append(position)
        mutated_reads.append(
            replace(
                read,
                sequence="".join(bases),
                error_positions=tuple(sorted(set(read.error_positions) | set(positions))),
            )
        )

    return mutated_reads


def add_terminal_damage(
    reads: Sequence[SimulatedRead],
    five_prime_ct: Sequence[Real],
    three_prime_ga: Sequence[Real],
    seed: int,
) -> list[SimulatedRead]:
    """Apply position-specific 5' C-to-T and 3' G-to-A damage."""
    five_prime_rates = [
        _validate_probability(rate, "five_prime_ct rate") for rate in five_prime_ct
    ]
    three_prime_rates = [
        _validate_probability(rate, "three_prime_ga rate") for rate in three_prime_ga
    ]
    random_generator = random.Random(seed)
    damaged_reads: list[SimulatedRead] = []

    for read in reads:
        bases = list(read.sequence)
        positions: set[int] = set(read.damage_positions)

        for distance, probability in enumerate(five_prime_rates):
            if distance >= len(bases):
                break
            if bases[distance] == "C" and random_generator.random() < probability:
                bases[distance] = "T"
                positions.add(distance)

        for distance, probability in enumerate(three_prime_rates):
            if distance >= len(bases):
                break
            position = len(bases) - 1 - distance
            if bases[position] == "G" and random_generator.random() < probability:
                bases[position] = "A"
                positions.add(position)

        damaged_reads.append(
            replace(
                read,
                sequence="".join(bases),
                damage_positions=tuple(sorted(positions)),
            )
        )

    return damaged_reads


def introduce_snp(reference: str, position: int, alternate: str) -> str:
    """Return a reference carrying one specified nucleotide substitution."""
    normalized = normalize_dna(reference)
    if not isinstance(position, int) or isinstance(position, bool):
        raise TypeError("position must be an integer")
    if not 0 <= position < len(normalized):
        raise ValueError("position is outside the reference")

    alternate_base = normalize_dna(alternate)
    if len(alternate_base) != 1 or alternate_base == "N":
        raise ValueError("alternate must be one unambiguous DNA base")
    if alternate_base == normalized[position]:
        raise ValueError("alternate must differ from the reference base")

    return normalized[:position] + alternate_base + normalized[position + 1 :]
