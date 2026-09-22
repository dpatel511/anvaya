import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))

import branch_transition_opportunity as audit
from anvaya.damage_string_graph import assemble
from anvaya.overlap_assembly import _MasterOverlapEdge
from anvaya.raw_consensus import DamageProfile
from anvaya.reads import Read


class BranchTransitionOpportunityTests(unittest.TestCase):
    def test_anchor_placements_match_exhaustive_oracle(self):
        rng = random.Random(73)
        sequence = "".join(rng.choice("ACGT") for _ in range(120))
        reads = [
            Read("first", sequence[10:90], (35,) * 80),
            Read("second", sequence[18:105], (35,) * 87),
            Read("outside", sequence[:35], (35,) * 35),
        ]
        profile = DamageProfile((), ())
        indexed = audit._placements(
            sequence, 35, 70, reads, profile, audit._raw_anchor_index(reads),
        )
        exhaustive = set()
        for read_index, read in enumerate(reads):
            for reverse in (False, True):
                for offset in range(len(sequence) - len(read.sequence) + 1):
                    if (
                        offset < 35
                        and offset + len(read.sequence) > 70
                        and audit._compatible(read, reverse, sequence, offset, profile)
                    ):
                        exhaustive.add((read_index, reverse, offset))
        self.assertEqual(indexed, exhaustive)

    def test_two_independent_molecules_resolve_one_transition(self):
        rng = random.Random(91)
        reference = "".join(rng.choice("ACGT") for _ in range(120))
        contigs = [
            Read("source", reference[0:55]),
            Read("centre", reference[35:80]),
            Read("target", reference[60:115]),
            Read("alternative", reference[60:80] + "A" * 35),
        ]
        incoming = _MasterOverlapEdge((0, False), (1, False), 35, 20, ())
        supported = _MasterOverlapEdge((1, False), (2, False), 25, 20, ())
        alternative = _MasterOverlapEdge((1, False), (3, False), 25, 20, ())
        reduced = {
            (incoming.source, incoming.target): incoming,
            (supported.source, supported.target): supported,
            (alternative.source, alternative.target): alternative,
        }
        bridge = reference[25:95]
        reads = [
            Read("bridge_1", bridge, (35,) * len(bridge)),
            Read("bridge_2", bridge, (35,) * len(bridge)),
        ]
        result = audit.audit_reduced_graph(
            contigs, reduced, set(), reads, [10, 11], DamageProfile((), ()),
        )
        self.assertEqual(result["possible_transitions"], 2)
        self.assertEqual(result["transitions_with_2_molecules"], 1)
        self.assertEqual(result["uniquely_resolved_transitions"], 1)
        self.assertGreater(result["projection_only_n50"], 0)

    def test_reduced_graph_callback_does_not_change_assembly(self):
        reads = [
            Read("a", "A" * 20 + "CGTACGTACGTACGTACGTACGTACGTACGT"),
            Read("b", "CGTACGTACGTACGTACGTACGTACGTACGT" + "C" * 20),
        ]
        baseline, baseline_diagnostics = assemble(reads)
        calls = []
        audited, audited_diagnostics = assemble(
            reads,
            reduced_graph_audit=lambda contigs, edges, excluded: calls.append(
                (len(contigs), len(edges), len(excluded))
            ),
        )
        self.assertEqual(calls and len(calls), 1)
        self.assertEqual(baseline, audited)
        self.assertEqual(baseline_diagnostics, audited_diagnostics)


if __name__ == "__main__":
    unittest.main()
