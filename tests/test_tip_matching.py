import unittest

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.tip_matching import match_tip_to_backbone


class TipBackboneMatchingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.graph = build_bidirected_dbg(
            ["AAGCCCAAA"] * 10 + ["AAGCCTAAA"],
            5,
            end_window=3,
        )
        self.tip = find_weak_tip_candidates(self.graph)[0]

    def test_matches_tip_to_strong_linear_backbone(self) -> None:
        match = match_tip_to_backbone(self.graph, self.tip)

        self.assertIsNotNone(match)
        assert match is not None
        self.assertEqual(match.tip_sequence, "AGCCTAAA")
        self.assertEqual(match.backbone_sequence, "AGCCCAAA")
        self.assertEqual(match.backbone_minimum_edge_support, 10)
        self.assertAlmostEqual(match.relative_coverage, 0.1)
        self.assertAlmostEqual(match.sequence_identity, 0.875)
        self.assertEqual(match.ry_identity, 1.0)

    def test_reports_terminal_damage_compatible_substitution(self) -> None:
        match = match_tip_to_backbone(self.graph, self.tip)

        self.assertIsNotNone(match)
        assert match is not None
        self.assertEqual(len(match.substitutions), 1)
        change = match.substitutions[0]
        self.assertEqual(
            (change.position, change.reference, change.alternative),
            (5, "C", "T"),
        )
        self.assertTrue(match.damage_compatible)

    def test_does_not_modify_graph(self) -> None:
        degrees = bytes(self.graph.out_degrees)
        active_edges = self.graph.active_edge_count
        observations = self.graph.observations

        match_tip_to_backbone(self.graph, self.tip)

        self.assertEqual(bytes(self.graph.out_degrees), degrees)
        self.assertEqual(self.graph.active_edge_count, active_edges)
        self.assertEqual(self.graph.observations, observations)

    def test_rejects_non_damage_substitution(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAT"] * 10 + ["AAGCCTAAA"],
            5,
            end_window=3,
        )
        tip = find_weak_tip_candidates(graph)[0]

        match = match_tip_to_backbone(graph, tip)

        self.assertIsNotNone(match)
        assert match is not None
        self.assertFalse(match.damage_compatible)

    def test_reports_three_prime_g_to_a_substitution(self) -> None:
        graph = build_bidirected_dbg(
            ["ACATACACG"] * 10 + ["ACATACACA"],
            5,
            end_window=3,
        )
        tip = find_weak_tip_candidates(graph)[0]

        match = match_tip_to_backbone(graph, tip)

        self.assertIsNotNone(match)
        assert match is not None
        self.assertEqual(
            [
                (change.position, change.reference, change.alternative)
                for change in match.substitutions
            ],
            [(5, "G", "A")],
        )
        self.assertTrue(match.damage_compatible)

    def test_accepts_multiple_oriented_damage_substitutions(self) -> None:
        graph = build_bidirected_dbg(
            ["TCTCGTGAAGCC"] * 10 + ["TCTTATGAAGCC"],
            5,
            end_window=4,
        )
        tip = find_weak_tip_candidates(graph)[0]

        match = match_tip_to_backbone(graph, tip)

        self.assertIsNotNone(match)
        assert match is not None
        self.assertEqual(
            [
                (change.position, change.reference, change.alternative)
                for change in match.substitutions
            ],
            [(5, "C", "T"), (6, "G", "A")],
        )
        self.assertTrue(match.damage_compatible)


if __name__ == "__main__":
    unittest.main()
