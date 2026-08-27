import math
import unittest

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.tip_classification import (
    TipClassificationThresholds,
    classify_tip_match,
    collect_tip_evidence,
)
from anvaya.tip_matching import match_tip_to_backbone


class TipClassificationTests(unittest.TestCase):
    def _classify(
        self,
        backbone: str,
        tip_read: str,
        end_window: int = 3,
    ):
        graph = build_bidirected_dbg(
            [backbone] * 10 + [tip_read],
            5,
            end_window=end_window,
        )
        tip = find_weak_tip_candidates(graph)[0]
        match = match_tip_to_backbone(graph, tip)
        self.assertIsNotNone(match)
        assert match is not None
        return graph, match, classify_tip_match(graph, match)

    def test_reports_exact_five_prime_damage_cycle(self) -> None:
        _, _, decision = self._classify(
            "AAGCCCAAA", "AAGCCTAAA", end_window=11
        )

        self.assertEqual(decision.label, "ambiguous")
        self.assertGreater(decision.damage_score, decision.error_score)
        self.assertEqual(decision.evidence.compatible_substitutions, 1)
        self.assertEqual(
            decision.evidence.substitution_terminal_observations,
            1,
        )
        self.assertEqual(decision.evidence.mean_damage_distance, 5.0)

    def test_classifies_three_prime_damage_supported_tip(self) -> None:
        _, _, decision = self._classify("ACATACACG", "ACATACACA")

        self.assertEqual(decision.label, "damage-like")
        self.assertEqual(decision.substitutions[0].expected_end, "right")
        self.assertEqual(decision.evidence.mean_damage_distance, 0.0)

    def test_reports_molecule_linked_substitution_support(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAA"] * 10 + ["AAGCCTAAA"],
            5,
            end_window=11,
            track_molecule_links=True,
        )
        tip = find_weak_tip_candidates(graph)[0]
        match = match_tip_to_backbone(graph, tip)
        self.assertIsNotNone(match)
        assert match is not None

        decision = classify_tip_match(graph, match)
        self.assertTrue(decision.evidence.molecule_links_collected)
        self.assertEqual(decision.evidence.joint_molecule_observations, 1)
        self.assertEqual(decision.evidence.joint_molecule_fraction, 1.0)
        self.assertNotIn("molecule_linked_substitutions", decision.reasons)

    def test_reports_quality_based_sequencing_error_likelihood(self) -> None:
        reads = ["AAGCCCAAA"] * 10 + ["AAGCCTAAA"]
        graph = build_bidirected_dbg(
            reads,
            5,
            end_window=11,
            track_molecule_links=True,
            read_qualities=[(40,) * 9] * 10 + [(10,) * 9],
        )
        tip = find_weak_tip_candidates(graph)[0]
        match = match_tip_to_backbone(graph, tip)
        self.assertIsNotNone(match)
        assert match is not None

        decision = classify_tip_match(graph, match)

        expected = math.log((10.0 ** -1.0) / 3.0)
        self.assertEqual(decision.evidence.quality_observations, 1)
        self.assertEqual(decision.evidence.missing_quality_observations, 0)
        self.assertEqual(decision.evidence.mean_base_quality, 10.0)
        self.assertAlmostEqual(
            decision.evidence.sequencing_error_log_likelihood,
            expected,
        )

    def test_missing_quality_does_not_invent_error_likelihood(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAA"] * 10 + ["AAGCCTAAA"],
            5,
            end_window=11,
            track_molecule_links=True,
        )
        tip = find_weak_tip_candidates(graph)[0]
        match = match_tip_to_backbone(graph, tip)
        self.assertIsNotNone(match)
        assert match is not None

        decision = classify_tip_match(graph, match)

        self.assertEqual(decision.evidence.quality_observations, 0)
        self.assertEqual(decision.evidence.missing_quality_observations, 1)
        self.assertIsNone(decision.evidence.mean_base_quality)
        self.assertIsNone(
            decision.evidence.sequencing_error_log_likelihood
        )

    def test_intersects_molecules_across_multiple_substitutions(self) -> None:
        graph = build_bidirected_dbg(
            ["TCTCGTGAAGCC"] * 10 + ["TCTTATGAAGCC"],
            5,
            end_window=7,
            track_molecule_links=True,
        )
        tip = find_weak_tip_candidates(graph)[0]
        match = match_tip_to_backbone(graph, tip)
        self.assertIsNotNone(match)
        assert match is not None

        decision = classify_tip_match(graph, match)
        aggregate_decision = classify_tip_match(
            graph,
            match,
            include_molecule_linkage=False,
        )

        self.assertEqual(len(decision.substitutions), 2)
        self.assertEqual(
            [change.expected_molecule_count for change in decision.substitutions],
            [0, 1],
        )
        self.assertEqual(decision.evidence.joint_molecule_observations, 0)
        self.assertEqual(decision.evidence.joint_molecule_fraction, 0.0)
        self.assertNotIn("molecule_linked_substitutions", decision.reasons)
        self.assertLess(
            decision.damage_score,
            aggregate_decision.damage_score,
        )
        self.assertNotIn(
            "molecule_linked_substitutions",
            aggregate_decision.reasons,
        )

    def test_rejects_invalid_molecule_linkage_flag(self) -> None:
        graph, match, _ = self._classify("AAGCCCAAA", "AAGCCTAAA")

        with self.assertRaisesRegex(TypeError, "must be a boolean"):
            classify_tip_match(
                graph,
                match,
                include_molecule_linkage=1,  # type: ignore[arg-type]
            )

    def test_leaves_mixed_substitution_tip_ambiguous(self) -> None:
        _, _, decision = self._classify("AAGCCCAAT", "AAGCCTAAA")

        self.assertEqual(decision.label, "ambiguous")
        self.assertEqual(
            decision.evidence.substitution_other_terminal_observations,
            1,
        )

    def test_does_not_call_terminally_supported_ordinary_change_error_like(
        self,
    ) -> None:
        _, _, decision = self._classify("AAGCCCAAA", "AAGCCGAAA")

        self.assertGreater(decision.error_score, 0.8)
        self.assertEqual(decision.label, "ambiguous")
        self.assertEqual(
            decision.reasons,
            ("terminal_depletion_not_observed",),
        )

    def test_thresholds_can_keep_a_supported_tip_ambiguous(self) -> None:
        graph, match, _ = self._classify("AAGCCCAAA", "AAGCCTAAA")

        decision = classify_tip_match(
            graph,
            match,
            TipClassificationThresholds(minimum_damage_score=0.95),
        )

        self.assertEqual(decision.label, "ambiguous")

    def test_evidence_collection_does_not_modify_graph(self) -> None:
        graph, match, _ = self._classify("AAGCCCAAA", "AAGCCTAAA")
        degrees = bytes(graph.out_degrees)
        observations = graph.observations

        collect_tip_evidence(graph, match)

        self.assertEqual(bytes(graph.out_degrees), degrees)
        self.assertEqual(graph.observations, observations)


if __name__ == "__main__":
    unittest.main()
