import unittest

from anvaya.bidirected import build_bidirected_dbg
from anvaya.bubbles import find_simple_bubbles
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.incomplete_branches import (
    find_incomplete_branch_candidates,
    match_incomplete_branch_to_backbone,
)
from anvaya.sequences import reverse_complement
from anvaya.tip_classification import classify_tip_match


class IncompleteBranchTests(unittest.TestCase):
    def setUp(self) -> None:
        self.graph = build_bidirected_dbg(
            ["AAGCCCAAA"] * 10 + ["AAGCCTAAA", "AAGCCTAAC"],
            5,
            end_window=3,
        )

    def test_finds_weak_branch_ending_at_another_branch(self) -> None:
        candidates = find_incomplete_branch_candidates(self.graph)

        self.assertEqual(len(candidates), 1)
        candidate = candidates[0]
        self.assertEqual(candidate.minimum_edge_support, 2)
        self.assertEqual(candidate.competing_minimum_edge_support, 10)
        self.assertGreater(self.graph.out_degrees[candidate.end], 1)

    def test_matches_branch_to_equal_length_backbone(self) -> None:
        candidate = find_incomplete_branch_candidates(self.graph)[0]

        match = match_incomplete_branch_to_backbone(self.graph, candidate)

        self.assertIsNotNone(match)
        assert match is not None
        self.assertEqual(match.tip_sequence, "AGCCTAA")
        self.assertEqual(match.backbone_sequence, "AGCCCAA")
        self.assertEqual(match.relative_coverage, 0.2)
        self.assertTrue(match.damage_compatible)
        self.assertEqual(classify_tip_match(self.graph, match).label, "damage-like")

    def test_is_distinct_from_tip_and_bubble_detection(self) -> None:
        candidate_edges = set(
            find_incomplete_branch_candidates(self.graph)[0].edge_ids
        )
        tip_edges = {
            edge_id
            for tip in find_weak_tip_candidates(self.graph)
            for edge_id in tip.edge_ids
        }

        self.assertTrue(candidate_edges.isdisjoint(tip_edges))
        self.assertEqual(find_simple_bubbles(self.graph), [])

    def test_does_not_report_a_simple_dead_end_tip(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAA"] * 10 + ["AAGCCTAAA"],
            5,
            end_window=3,
        )

        self.assertEqual(find_incomplete_branch_candidates(graph), [])

    def test_does_not_report_a_rejoining_bubble(self) -> None:
        graph = build_bidirected_dbg(
            ["GCTTGTTCCGGA"] * 10 + ["GCTTATTCCGGA"] * 2,
            4,
            end_window=3,
        )

        self.assertEqual(find_incomplete_branch_candidates(graph), [])

    def test_detection_threshold_is_broader_than_tip_cleaning(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAA"] * 10
            + ["AAGCCTAAA", "AAGCCTAAC"] * 2,
            5,
            end_window=3,
        )

        self.assertEqual(len(find_incomplete_branch_candidates(graph)), 1)
        self.assertEqual(
            find_incomplete_branch_candidates(
                graph,
                maximum_relative_support=0.3,
            ),
            [],
        )

    def test_reverse_complements_do_not_duplicate_candidate(self) -> None:
        major = "AAGCCCAAA"
        weak_left = "AAGCCTAAA"
        weak_right = "AAGCCTAAC"
        graph = build_bidirected_dbg(
            [major, reverse_complement(major)] * 10
            + [
                weak_left,
                reverse_complement(weak_left),
                weak_right,
                reverse_complement(weak_right),
            ],
            5,
            end_window=3,
        )

        self.assertEqual(len(find_incomplete_branch_candidates(graph)), 1)

    def test_detection_does_not_modify_graph(self) -> None:
        degrees = bytes(self.graph.out_degrees)
        observations = self.graph.observations
        active_edges = self.graph.active_edge_count

        find_incomplete_branch_candidates(self.graph)

        self.assertEqual(bytes(self.graph.out_degrees), degrees)
        self.assertEqual(self.graph.observations, observations)
        self.assertEqual(self.graph.active_edge_count, active_edges)

    def test_rejects_invalid_limit(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least 1"):
            find_incomplete_branch_candidates(self.graph, max_edges=0)
        with self.assertRaisesRegex(TypeError, "integer"):
            find_incomplete_branch_candidates(
                self.graph,
                max_edges=1.5,  # type: ignore[arg-type]
            )

    def test_rejects_invalid_relative_support(self) -> None:
        with self.assertRaisesRegex(ValueError, "between 0 and 1"):
            find_incomplete_branch_candidates(
                self.graph,
                maximum_relative_support=1.1,
            )
        with self.assertRaisesRegex(TypeError, "numeric"):
            find_incomplete_branch_candidates(
                self.graph,
                maximum_relative_support="0.5",  # type: ignore[arg-type]
            )


if __name__ == "__main__":
    unittest.main()
