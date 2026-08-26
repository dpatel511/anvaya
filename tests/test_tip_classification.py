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
    ):
        graph = build_bidirected_dbg(
            [backbone] * 10 + [tip_read],
            5,
            end_window=3,
        )
        tip = find_weak_tip_candidates(graph)[0]
        match = match_tip_to_backbone(graph, tip)
        self.assertIsNotNone(match)
        assert match is not None
        return graph, match, classify_tip_match(graph, match)

    def test_classifies_five_prime_damage_supported_tip(self) -> None:
        _, _, decision = self._classify("AAGCCCAAA", "AAGCCTAAA")

        self.assertEqual(decision.label, "damage-like")
        self.assertGreater(decision.damage_score, decision.error_score)
        self.assertEqual(decision.evidence.compatible_substitutions, 1)
        self.assertEqual(
            decision.evidence.substitution_terminal_observations,
            1,
        )
        self.assertEqual(decision.evidence.mean_damage_distance, 1.0)

    def test_classifies_three_prime_damage_supported_tip(self) -> None:
        _, _, decision = self._classify("ACATACACG", "ACATACACA")

        self.assertEqual(decision.label, "damage-like")
        self.assertEqual(decision.substitutions[0].expected_end, "right")
        self.assertEqual(decision.evidence.mean_damage_distance, 0.0)

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
