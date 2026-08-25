"""Development smoke test for damage-induced graph fragmentation."""

import random

from helpers.simulation import SimulatedRead, add_terminal_damage
from anvaya.assembly import assemble
from anvaya.graph import build_dbg
from anvaya.metrics import summarize_graph


SEED = 42
REFERENCE_LENGTH = 500
READ_LENGTH = 60
STEP = 10
K = 21
FIVE_PRIME_CT = [0.50, 0.30, 0.15, 0.08, 0.04]
THREE_PRIME_GA = [0.50, 0.30, 0.15, 0.08, 0.04]


def main() -> None:
    random_generator = random.Random(SEED)
    reference = "".join(
        random_generator.choice("ACGT") for _ in range(REFERENCE_LENGTH)
    )

    starts = list(range(0, len(reference) - READ_LENGTH + 1, STEP))
    if starts[-1] != len(reference) - READ_LENGTH:
        starts.append(len(reference) - READ_LENGTH)

    clean = [
        SimulatedRead(
            name=f"read-{index}",
            sequence=reference[start : start + READ_LENGTH],
            source_sequence=reference[start : start + READ_LENGTH],
            reference_start=start,
            reference_end=start + READ_LENGTH,
        )
        for index, start in enumerate(starts)
    ]
    damaged = add_terminal_damage(
        clean,
        five_prime_ct=FIVE_PRIME_CT,
        three_prime_ga=THREE_PRIME_GA,
        seed=SEED + 1,
    )

    clean_sequences = [read.sequence for read in clean]
    damaged_sequences = [read.sequence for read in damaged]
    clean_unitigs = assemble(clean_sequences, K)
    damaged_unitigs = assemble(damaged_sequences, K)

    print("Clean:", summarize_graph(build_dbg(clean_sequences, K)))
    print("Damaged:", summarize_graph(build_dbg(damaged_sequences, K)))
    print("Clean exact reconstruction:", clean_unitigs == [reference])
    print("Longest clean unitig:", max(map(len, clean_unitigs)))
    print("Longest damaged unitig:", max(map(len, damaged_unitigs)))
    print("Damage events:", sum(len(read.damage_positions) for read in damaged))


if __name__ == "__main__":
    main()
