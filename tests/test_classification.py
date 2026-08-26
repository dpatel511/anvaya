import unittest
from dataclasses import dataclass, replace

from anvaya.classification import classify_unitig_path, substitutions


@dataclass(frozen=True)
class Evidence:
    left_terminal_observations: int = 0
    right_terminal_observations: int = 0
    relative_coverage: float = 0.2
    strand_balance: float | None = 0.5
    terminal_enrichment: float = -0.2
    similarity_to_strongest: float = 0.98


class PathClassificationTests(unittest.TestCase):
    def test_labels_dominant_path(self) -> None:
        decision = classify_unitig_path(
            "ACGT",
            "ACGT",
            Evidence(),
            dominant=True,
        )

        self.assertEqual(decision.label, "dominant")

    def test_labels_terminally_supported_damage_substitution(self) -> None:
        evidence = Evidence(
            left_terminal_observations=3,
            terminal_enrichment=0.2,
        )

        decision = classify_unitig_path("CCAA", "TCAA", evidence)

        self.assertEqual(decision.label, "damage-like")
        self.assertTrue(decision.damage_compatible)

    def test_requires_correct_end_for_damage_substitution(self) -> None:
        evidence = Evidence(
            right_terminal_observations=3,
            terminal_enrichment=0.2,
        )

        decision = classify_unitig_path("CCAA", "TCAA", evidence)

        self.assertFalse(decision.damage_compatible)
        self.assertEqual(decision.label, "ambiguous")

    def test_does_not_call_supported_variation_damage_like(self) -> None:
        evidence = Evidence(
            left_terminal_observations=3,
            relative_coverage=0.6,
            strand_balance=0.8,
            terminal_enrichment=0.2,
        )

        decision = classify_unitig_path("CCAA", "TCAA", evidence)

        self.assertTrue(decision.damage_compatible)
        self.assertEqual(decision.label, "variation-like")

    def test_labels_low_coverage_terminally_depleted_path_as_error_like(self) -> None:
        decision = classify_unitig_path("ACGT", "ATGT", Evidence())

        self.assertEqual(decision.label, "error-like")
        self.assertIn("terminal_depleted", decision.reasons)

    def test_labels_balanced_supported_path_as_variation_like(self) -> None:
        evidence = replace(
            Evidence(),
            relative_coverage=0.5,
            strand_balance=0.8,
            terminal_enrichment=0.0,
        )

        decision = classify_unitig_path("ACGT", "ATGT", evidence)

        self.assertEqual(decision.label, "variation-like")

    def test_keeps_conflicting_evidence_ambiguous(self) -> None:
        evidence = replace(Evidence(), terminal_enrichment=0.0)

        decision = classify_unitig_path("ACGT", "ATGT", evidence)

        self.assertEqual(decision.label, "ambiguous")

    def test_does_not_call_substitutions_for_different_lengths(self) -> None:
        self.assertEqual(substitutions("ACGT", "ACGTT"), ())


if __name__ == "__main__":
    unittest.main()
