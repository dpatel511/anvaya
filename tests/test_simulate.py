import unittest

from helpers.simulation import (
    SimulatedRead,
    add_sequencing_errors,
    add_terminal_damage,
    introduce_snp,
    simulate_fragments,
)


class SimulateFragmentsTests(unittest.TestCase):
    def test_generates_requested_coverage_and_coordinates(self) -> None:
        reads = simulate_fragments("ACGT" * 25, read_length=20, coverage=2, seed=7)

        self.assertEqual(len(reads), 10)
        for read in reads:
            self.assertEqual(len(read.sequence), 20)
            self.assertEqual(read.reference_end - read.reference_start, 20)
            self.assertEqual(read.sequence, read.source_sequence)

    def test_is_reproducible_for_fixed_seed(self) -> None:
        first = simulate_fragments("ACGT" * 25, 20, 2, seed=11)
        second = simulate_fragments("ACGT" * 25, 20, 2, seed=11)

        self.assertEqual(first, second)


class SequencingErrorTests(unittest.TestCase):
    def test_error_rate_one_changes_every_base(self) -> None:
        read = SimulatedRead("read-0", "AAAA", "AAAA", 0, 4)

        mutated = add_sequencing_errors([read], error_rate=1, seed=3)[0]

        self.assertNotIn("A", mutated.sequence)
        self.assertEqual(mutated.error_positions, (0, 1, 2, 3))
        self.assertEqual(mutated.source_sequence, "AAAA")

    def test_error_rate_zero_preserves_read(self) -> None:
        read = SimulatedRead("read-0", "ACGT", "ACGT", 0, 4)

        self.assertEqual(add_sequencing_errors([read], 0, seed=3), [read])


class TerminalDamageTests(unittest.TestCase):
    def test_applies_position_specific_terminal_damage(self) -> None:
        read = SimulatedRead("read-0", "CCAAAAGG", "CCAAAAGG", 0, 8)

        damaged = add_terminal_damage(
            [read],
            five_prime_ct=[1, 0],
            three_prime_ga=[1, 0],
            seed=5,
        )[0]

        self.assertEqual(damaged.sequence, "TCAAAAGA")
        self.assertEqual(damaged.damage_positions, (0, 7))
        self.assertEqual(damaged.source_sequence, "CCAAAAGG")

    def test_non_target_bases_are_unchanged(self) -> None:
        read = SimulatedRead("read-0", "AAAA", "AAAA", 0, 4)

        self.assertEqual(
            add_terminal_damage([read], [1, 1], [1, 1], seed=5),
            [read],
        )


class IntroduceSnpTests(unittest.TestCase):
    def test_introduces_requested_snp(self) -> None:
        self.assertEqual(introduce_snp("ACGT", 1, "T"), "ATGT")

    def test_rejects_unchanged_allele(self) -> None:
        with self.assertRaises(ValueError):
            introduce_snp("ACGT", 1, "C")


if __name__ == "__main__":
    unittest.main()
