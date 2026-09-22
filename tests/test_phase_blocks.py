"""Synthetic controls for truth-blind, damage-aware phase-link evidence."""

import unittest

from anvaya.phase_blocks import phase_blocks
from anvaya.raw_consensus import DamageProfile, observation_likelihoods


BASES = "ACGT"


def observation(base, position=10, length=30, profile=DamageProfile((), ()), reverse=False):
    quality = 38
    raw_base = BASES[3 - BASES.index(base)] if reverse else base
    return (
        BASES.index(base), quality,
        observation_likelihoods(raw_base, quality, position, length, reverse, profile),
    )


class PhaseBlockTests(unittest.TestCase):
    def test_linked_allele_pairs_form_one_block_and_retain_molecules(self):
        columns = [{} for _ in range(11)]
        columns[2] = {m: observation("A" if m < 2 else "C") for m in range(4)}
        columns[10] = {m: observation("G" if m < 2 else "T") for m in range(4)}
        evidence = phase_blocks(columns)
        self.assertEqual(evidence.blocks, ((2, 10),))
        self.assertEqual(evidence.links[0].allele_pair_support, ((0, 2, 2), (1, 3, 2)))
        self.assertEqual(
            evidence.links[0].allele_pair_molecules,
            ((0, 2, (0, 1)), (1, 3, (2, 3))),
        )

    def test_unlinked_haplotypes_form_separate_blocks(self):
        columns = [{} for _ in range(11)]
        columns[2] = {m: observation("A" if m < 2 else "C") for m in range(4)}
        columns[10] = {m: observation("G" if m < 6 else "T") for m in range(4, 8)}
        evidence = phase_blocks(columns)
        self.assertEqual(evidence.blocks, ((2,), (10,)))
        self.assertEqual(evidence.unlinked_pairs_0based, ((2, 10),))
        self.assertEqual(
            evidence.marker_alleles,
            ((2, 0, (0, 1)), (2, 1, (2, 3)),
             (10, 2, (4, 5)), (10, 3, (6, 7))),
        )

    def test_singleton_third_allele_does_not_veto_supported_link(self):
        columns = [{} for _ in range(11)]
        columns[2] = {m: observation("A" if m < 2 else "C") for m in range(4)}
        columns[2][8] = observation("G")
        columns[10] = {m: observation("G" if m < 2 else "T") for m in range(4)}
        columns[10][8] = observation("G")
        evidence = phase_blocks(columns)
        self.assertEqual(evidence.blocks, ((2, 10),))
        self.assertEqual(evidence.excluded_alleles, ((2, 2, (8,)),))

    def test_all_four_pairs_are_ambiguous_not_a_resolved_link(self):
        columns = [{} for _ in range(2)]
        pairs = (("A", "G"), ("A", "T"), ("C", "G"), ("C", "T"))
        for pair_index, (left, right) in enumerate(pairs):
            for molecule in range(2 * pair_index, 2 * pair_index + 2):
                columns[0][molecule] = observation(left)
                columns[1][molecule] = observation(right)
        evidence = phase_blocks(columns)
        self.assertEqual(evidence.blocks, ((0,), (1,)))
        self.assertEqual(evidence.links, ())
        self.assertEqual(len(evidence.ambiguous_links), 1)
        self.assertEqual(
            evidence.ambiguous_links[0].allele_pair_support,
            ((0, 2, 2), (0, 3, 2), (1, 2, 2), (1, 3, 2)),
        )

    def test_nonadjacent_link_is_retained(self):
        columns = [{} for _ in range(3)]
        columns[0] = {m: observation("A" if m < 2 else "C") for m in range(4)}
        columns[1] = {m: observation("G" if m < 6 else "T") for m in range(4, 8)}
        columns[2] = {m: observation("G" if m < 2 else "T") for m in range(4)}
        evidence = phase_blocks(columns)
        self.assertEqual(evidence.blocks, ((0, 2), (1,)))
        self.assertEqual(
            [(link.left_position_0based, link.right_position_0based) for link in evidence.links],
            [(0, 2)],
        )

    def test_three_site_switching_retains_pair_specific_molecules(self):
        columns = [{} for _ in range(3)]
        columns[0] = {m: observation("A" if m < 2 else "C") for m in range(4)}
        columns[1] = {m: observation("G" if m < 2 else "T") for m in range(4)}
        columns[1].update({m: observation("G" if m < 6 else "T") for m in range(4, 8)})
        columns[2] = {m: observation("C" if m < 6 else "A") for m in range(4, 8)}
        evidence = phase_blocks(columns)
        self.assertEqual(evidence.blocks, ((0, 1, 2),))
        self.assertEqual(evidence.links[0].allele_pair_molecules[0][2], (0, 1))
        self.assertEqual(evidence.links[1].allele_pair_molecules[0][2], (4, 5))
        self.assertIn((0, 2), evidence.unlinked_pairs_0based)

    def test_inconsistent_three_site_switching_is_ambiguous(self):
        columns = [{} for _ in range(3)]
        pairs = (
            (0, 1, 0, (("A", "G"), ("C", "T"))),
            (1, 2, 4, (("G", "C"), ("T", "A"))),
            (0, 2, 8, (("A", "A"), ("C", "C"))),
        )
        for left, right, first_molecule, allele_pairs in pairs:
            for pair_index, (left_allele, right_allele) in enumerate(allele_pairs):
                for molecule in range(first_molecule + 2 * pair_index,
                                      first_molecule + 2 * pair_index + 2):
                    columns[left][molecule] = observation(left_allele)
                    columns[right][molecule] = observation(right_allele)
        evidence = phase_blocks(columns)
        self.assertEqual(evidence.blocks, ((0,), (1,), (2,)))
        self.assertEqual(evidence.links, ())
        self.assertEqual(
            [(link.left_position_0based, link.right_position_0based)
             for link in evidence.ambiguous_links],
            [(0, 1), (0, 2), (1, 2)],
        )

    def test_conflicting_molecule_and_singleton_error_are_excluded(self):
        columns = [{} for _ in range(2)]
        columns[0] = {m: observation("A" if m < 2 else "C") for m in range(4)}
        columns[1] = {m: observation("G" if m < 2 else "T") for m in range(4)}
        columns[0][8] = None
        columns[1][8] = observation("A")
        evidence = phase_blocks(columns)
        self.assertEqual(evidence.blocks, ((0, 1),))
        self.assertEqual(evidence.excluded_alleles, ((1, 0, (8,)),))
        self.assertNotIn(8, evidence.links[0].allele_pair_molecules[0][2])

    def test_reverse_observations_use_contig_oriented_alleles(self):
        columns = [{} for _ in range(2)]
        columns[0] = {m: observation("A" if m < 2 else "C", reverse=bool(m % 2))
                      for m in range(4)}
        columns[1] = {m: observation("G" if m < 2 else "T", reverse=bool(m % 2))
                      for m in range(4)}
        evidence = phase_blocks(columns)
        self.assertEqual(evidence.blocks, ((0, 1),))
        self.assertEqual(evidence.links[0].allele_pair_support, ((0, 2, 2), (1, 3, 2)))

    def test_terminal_damage_does_not_create_variable_site(self):
        profile = DamageProfile((0.4,), ())
        column = {
            0: observation("C", 0, profile=profile),
            1: observation("C", 0, profile=profile),
            2: observation("T", 0, profile=profile),
            3: observation("T", 0, profile=profile),
        }
        self.assertEqual(phase_blocks([column]).variable_positions_0based, ())

    def test_singleton_error_does_not_create_variable_site(self):
        column = {m: observation("A") for m in range(4)}
        column[4] = observation("C")
        self.assertEqual(phase_blocks([column]).variable_positions_0based, ())


if __name__ == "__main__":
    unittest.main()
