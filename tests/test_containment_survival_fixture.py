"""Regression for containment changing damage-correction path survival."""
import unittest

from anvaya.damage_string_graph import assemble
from anvaya.raw_consensus import DamageProfile
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


class ContainmentSurvivalFixtureTests(unittest.TestCase):
    def test_contained_read_remains_available_as_damage_correction_partner(self):
        # Reference coordinates are 0-based and intervals are half-open.
        outer = (
            "CTACTATCATTCTCAGCCAGCGACTATCATACTGTCAAACCAGCGACGTGGGGAGGAGATCGTC"
            "GCGCCGTGCATACGGAATGCTTACCCCGATTGTTCGTCATAAGGGGCACGAACGATGCGTGGT"
            "CTAAACCCCCGGTACGAACTCTTGACGAATCTGAAGTCATACGTTCGCTGTCAGGAT"
        )
        contained = (
            "ATCATACTGTCAAACCAGCGACGTGGGGAGGAGATCGTCGCGCCGTGCATACGGAATGCTTACC"
            "CCGATTGTTCGTCATAAGGGGCACGAACGATGCGTGGTCTAAACCCCCGGTACGAACTCTTGAC"
            "GAATCTGAAG"
        )
        outer_partner = (
            "ATTCCTTTGTCTACTATCATTCTCAGCCAGCGACTATCATACTGTCAAACCAGCGACGTGGGGAG"
            "GAGATCGTCGCGCCGTGCATACGGAATGCTTACCCCGATTGTTCGTCATAAGGGGCACGAACGA"
            "TGCGTGGTCTAAACCCCCGGTACGAACTCTTGACGAATCTGAAGTCATAC"
        )
        damaged_partner = (
            "TTGTCAAGAGTTCGTACCGGGGGTTTAGACCACGCATCGTTCGTGCCCCTTATGACGAACAATC"
            "GGGGTAAGCATTCCGTATGCACGGCGCGACGATCTCCTCCCCACGTCGCTGGTTTGACAGTATGA"
            "TAGTCGC"
        )
        reads = [
            Read("outer", outer, (35,) * len(outer)),
            Read("contained", contained, (35,) * len(contained)),
            Read("outer_partner", outer_partner, (35,) * len(outer_partner)),
            Read("damaged_partner", damaged_partner, (35,) * len(damaged_partner)),
        ]
        reference = reverse_complement(outer)
        partner_reference = reverse_complement(outer_partner)
        self.assertEqual(reference[15:], partner_reference[:-10])
        reference += partner_reference[-10:]
        self.assertEqual(reverse_complement(contained), reference[21:159])
        expected_partner = reference[29:165]
        self.assertEqual(damaged_partner[0], expected_partner[0])
        self.assertEqual(damaged_partner[1], "T")
        self.assertEqual(expected_partner[1], "C")
        self.assertEqual(damaged_partner[2:], expected_partner[2:])

        rates = tuple(0.25 * 0.65 ** position for position in range(5))
        pool, diagnostics = assemble(
            reads, damage_profile=DamageProfile(rates, rates),
        )

        by_placements = {
            frozenset(placement.read_index for placement in record.raw_placements):
            record.current.sequence
            for record in pool.active_derived
        }
        self.assertEqual(set(by_placements), {frozenset({0, 2}), frozenset({1, 3})})
        self.assertEqual(by_placements[frozenset({0, 2})], reference)
        self.assertEqual(by_placements[frozenset({1, 3})], reference[21:165])
        self.assertEqual(diagnostics["contained_nodes"], 0)
        self.assertEqual(diagnostics["ambiguous_containments"], 1)
        self.assertEqual(diagnostics["directed_edges"], 4)
        self.assertEqual(diagnostics["linear_paths"], 2)


if __name__ == "__main__":
    unittest.main()
