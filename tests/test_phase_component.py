"""Bounded component-adapter fixtures with truth kept outside the adapter."""

import unittest
from dataclasses import replace

from anvaya.damage_string_graph import assemble
from anvaya.overlap_progressive import ProgressiveSequencePool
from anvaya.phase_blocks import phase_component
from anvaya.raw_consensus import DamageProfile, RawPlacement
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


def component_pool(reads, molecule_ids, contig, placements):
    pool = ProgressiveSequencePool.from_reads(reads, molecule_ids)
    record = replace(
        pool.records[0].corrected(Read("component", contig)),
        raw_placements=tuple(placements),
    )
    return pool.replace_record(record), record


class PhaseComponentAdapterTests(unittest.TestCase):
    def test_graph_output_preserves_private_chain_coordinates(self):
        reference = (
            "AACATAGACTTTCCGGTCATTCGCTTCCACATGTAACAACCTATGTAACGGGGCTCGCGC"
            "TCGGCCCATTCGATTGAATTCGACCGACGCTCTCCTTGGTACGTTGCCAAATCGTCGTGA"
        )
        starts = (0, 30, 60)
        reads = [Read("forward-left", reference[:60], (38,) * 60),
                 Read("reverse-middle", reverse_complement(reference[30:90]), (38,) * 60),
                 Read("forward-right", reference[60:], (38,) * 60)]

        pool, diagnostics = assemble(reads, [20, 21, 22], maximum_reads=10)

        self.assertEqual(len(pool.active_derived), 1)
        record = pool.active_derived[0]
        self.assertEqual(record.current.sequence, reference)
        self.assertEqual(pool.raw_evidence, tuple(reads))
        self.assertEqual(
            tuple((placement.read_index, placement.offset, placement.reverse)
                  for placement in record.raw_placements),
            ((0, starts[0], False), (1, starts[1], True), (2, starts[2], False)),
        )
        self.assertEqual(record.contributing_molecules, frozenset({20, 21, 22}))
        self.assertEqual(diagnostics["ambiguous_containments"], 0)
        adapted = phase_component(pool, record, DamageProfile((), ()))
        self.assertEqual(adapted.diagnostics.consistent_placements, 3)
        self.assertEqual(adapted.diagnostics.ambiguous_placements, 0)
        self.assertEqual(adapted.diagnostics.covered_sites, len(reference))

    def test_reverse_containment_has_consistent_prelayout_coordinates(self):
        parent = (
            "AACATAGACTTTCCGGTCATTCGCTTCCACATGTAACAACCTATGTAACGGGGCTCGCGC"
        )
        child = reverse_complement(parent[11:34])
        reads = [
            Read("parent", parent, (38,) * len(parent)),
            Read("reverse-child", child, (38,) * len(child)),
        ]

        _, diagnostics = assemble(reads, [1, 2], maximum_reads=10)

        self.assertEqual(diagnostics["contained_nodes"], 1)
        self.assertEqual(diagnostics["prelayout_phase_components"], 1)
        self.assertEqual(diagnostics["prelayout_phase_coordinate_consistent"], 1)
        self.assertEqual(diagnostics["prelayout_phase_coordinate_unresolved"], 0)

    def test_linked_branch_recovers_private_molecule_truth_on_both_strands(self):
        first = "AAAAAAAAAGAA"
        second = "AACAAAAAATAA"
        oriented = [first, first, second, second]
        reads = [
            Read(str(index), reverse_complement(sequence) if index % 2 else sequence, (38,) * 12)
            for index, sequence in enumerate(oriented)
        ]
        placements = [RawPlacement(index, 0, bool(index % 2), 0, 12) for index in range(4)]
        pool, record = component_pool(reads, list(range(4)), first, placements)

        evidence = phase_component(pool, record, DamageProfile((), ()))

        self.assertEqual(record.current.sequence, first)
        self.assertEqual(pool.raw_evidence, tuple(reads))
        self.assertEqual(evidence.phase.variable_positions_0based, (2, 9))
        self.assertEqual(evidence.phase.blocks, ((2, 9),))
        self.assertEqual(
            evidence.phase.links[0].allele_pair_molecules,
            ((0, 2, (0, 1)), (1, 3, (2, 3))),
        )
        self.assertEqual(evidence.diagnostics.distinct_observed_molecules, 4)
        self.assertEqual(evidence.diagnostics.ambiguous_placements, 0)

    def test_unlinked_branch_remains_separate_and_reports_missing_sites(self):
        reads = [Read(str(index), base, (38,)) for index, base in enumerate("AACC GGTT".replace(" ", ""))]
        placements = [RawPlacement(index, 2 if index < 4 else 9, False, 0, 1)
                      for index in range(8)]
        pool, record = component_pool(reads, list(range(8)), "A" * 12, placements)

        evidence = phase_component(pool, record, DamageProfile((), ()))

        self.assertEqual(evidence.phase.blocks, ((2,), (9,)))
        self.assertEqual(evidence.phase.unlinked_pairs_0based, ((2, 9),))
        self.assertEqual(evidence.diagnostics.covered_sites, 2)
        self.assertEqual(evidence.diagnostics.uncovered_sites, 10)

    def test_repeat_placement_is_ambiguous_and_contributes_no_observation(self):
        reads = [Read("repeat", "ACGT", (38,) * 4)]
        placements = (
            RawPlacement(0, 0, False, 0, 4),
            RawPlacement(0, 4, False, 0, 4),
        )
        pool, record = component_pool(reads, [7], "ACGTACGT", placements)

        evidence = phase_component(pool, record, DamageProfile((), ()))

        self.assertEqual(record.current.sequence, "ACGTACGT")
        self.assertEqual(pool.raw_evidence, tuple(reads))
        self.assertEqual([item.status for item in evidence.placements], ["ambiguous", "ambiguous"])
        self.assertEqual(evidence.diagnostics.ambiguous_molecules, 1)
        self.assertEqual(evidence.diagnostics.ambiguous_placements, 2)
        self.assertEqual(evidence.diagnostics.distinct_observed_molecules, 0)
        self.assertEqual(evidence.diagnostics.uncovered_sites, 8)
        self.assertEqual(evidence.phase.variable_positions_0based, ())

    def test_missing_low_quality_ambiguous_and_outside_bases_are_reported(self):
        reads = [
            Read("missing", "AAAA", None),
            Read("low", "CCCC", (0,) * 4),
            Read("ambiguous", "NNNN", (38,) * 4),
            Read("outside", "GGGG", (38,) * 4),
        ]
        placements = (
            RawPlacement(0, 0, False, 0, 4),
            RawPlacement(1, 0, False, 0, 4),
            RawPlacement(2, 0, False, 0, 4),
            RawPlacement(3, 3, False, 0, 4),
        )
        pool, record = component_pool(reads, list(range(4)), "A" * 4, placements)

        evidence = phase_component(pool, record, DamageProfile((), ()))

        self.assertEqual(evidence.diagnostics.missing_quality_observations, 4)
        self.assertEqual(evidence.diagnostics.low_quality_observations, 4)
        self.assertEqual(evidence.diagnostics.ambiguous_base_observations, 4)
        self.assertEqual(evidence.diagnostics.outside_contig_observations, 3)
        self.assertEqual(evidence.diagnostics.usable_observations, 1)
        self.assertEqual(evidence.diagnostics.covered_sites, 1)

    def test_duplicate_molecule_conflict_is_excluded(self):
        reads = [Read("a", "A", (38,)), Read("c", "C", (38,))]
        placements = (RawPlacement(0, 0, False, 0, 1), RawPlacement(1, 0, False, 0, 1))
        pool, record = component_pool(reads, [3, 3], "A", placements)

        evidence = phase_component(pool, record, DamageProfile((), ()))

        self.assertEqual(evidence.diagnostics.duplicate_molecule_observations, 1)
        self.assertEqual(evidence.diagnostics.conflicting_molecule_observations, 1)
        self.assertEqual(evidence.diagnostics.distinct_observed_molecules, 0)
        self.assertEqual(evidence.diagnostics.covered_sites, 0)


if __name__ == "__main__":
    unittest.main()
