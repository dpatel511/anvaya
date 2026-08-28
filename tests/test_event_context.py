import unittest

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.event_context import extract_event_graph_context
from anvaya.tip_matching import match_tip_to_backbone


class EventGraphContextTests(unittest.TestCase):
    def test_extracts_local_path_and_evidence_features(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAA"] * 10 + ["AAGCCTAAA"] * 2,
            5,
            end_window=11,
            track_molecule_links=True,
            read_qualities=[(35,) * 9] * 12,
        )
        candidate = find_weak_tip_candidates(graph)[0]
        match = match_tip_to_backbone(graph, candidate)
        self.assertIsNotNone(match)
        assert match is not None

        context = extract_event_graph_context(graph, match, "tip")

        self.assertEqual(context.event_type, "tip")
        self.assertEqual(context.path_edge_count, len(match.tip_edge_ids))
        self.assertEqual(context.substitution_count, len(match.substitutions))
        self.assertEqual(
            context.alternative_path_observations,
            match.tip_observations,
        )
        self.assertEqual(
            context.reference_path_observations,
            match.backbone_observations,
        )
        self.assertGreater(context.alternative_minimum_edge_support, 0)
        self.assertEqual(
            context.reference_minimum_edge_support,
            match.backbone_minimum_edge_support,
        )
        self.assertAlmostEqual(context.relative_coverage, match.relative_coverage)


if __name__ == "__main__":
    unittest.main()
