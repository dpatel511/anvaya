import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))

import paired_branch_opportunity as audit
from anvaya.overlap_assembly import _MasterOverlapEdge
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


class PairedBranchOpportunityTests(unittest.TestCase):
    def test_two_independent_pairs_resolve_one_transition(self):
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
        left = [Read(f"pair_{i}/1", reference[10:40]) for i in range(2)]
        right = [
            Read(f"pair_{i}/2", reverse_complement(reference[75:105]))
            for i in range(2)
        ]
        result = audit.audit_paired_transitions(
            contigs, reduced, set(), left, right,
        )
        self.assertEqual(result["possible_transitions"], 2)
        self.assertEqual(result["transitions_with_2_pairs"], 1)
        self.assertEqual(result["uniquely_resolved_transitions"], 1)

    def test_pair_name_normalizes_slash_suffix_only(self):
        self.assertEqual(audit._pair_name("molecule/1 extra"), "molecule")
        self.assertEqual(audit._pair_name("molecule/2"), "molecule")
        self.assertEqual(audit._pair_name("molecule.1"), "molecule.1")


if __name__ == "__main__":
    unittest.main()
