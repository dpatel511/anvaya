import importlib.util
import sys
import unittest
from pathlib import Path


def _load_experiment():
    path = (
        Path(__file__).resolve().parents[1]
        / "experiments"
        / "05_unitig_bubble_validation.py"
    )
    specification = importlib.util.spec_from_file_location(
        "unitig_bubble_validation",
        path,
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("could not load validation experiment")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


experiment = _load_experiment()


class UnitigBubbleValidationTests(unittest.TestCase):
    def test_labels_major_rare_shared_and_artifact_paths(self) -> None:
        major = "ACGTACTGGAACCTTA"
        rare = "ACGTGCTGGAACCTTA"

        self.assertEqual(
            experiment.label_path("rare", "ACGTACTG", major, rare),
            "major",
        )
        self.assertEqual(
            experiment.label_path("rare", "ACGTGCTG", major, rare),
            "rare",
        )
        self.assertEqual(
            experiment.label_path("rare", "GGAACC", major, rare),
            "shared",
        )
        self.assertEqual(
            experiment.label_path("clean", "ACACAC", major, None),
            "artifact",
        )

    def test_labels_nonreference_paths_by_introduced_process(self) -> None:
        reference = "AAAACCCCGGGGTTTT"

        self.assertEqual(
            experiment.label_path("error", "ACACAC", reference, None),
            "error",
        )
        self.assertEqual(
            experiment.label_path("damage", "ACACAC", reference, None),
            "damage",
        )

    def test_threshold_requires_both_coverage_and_similarity(self) -> None:
        candidate = experiment.CandidateRecord(
            condition="error",
            scenario="error_0.005",
            seed=1,
            bubble_index=1,
            path_index=0,
            label="error",
            classification="error-like",
            classification_reasons="low_local_coverage",
            damage_compatible=False,
            sequence="AACT",
            nucleotide_length=4,
            edge_count=1,
            observations=2,
            minimum_edge_support=2,
            mean_edge_support=2.0,
            terminal_fraction=0.2,
            strand_balance=0.1,
            relative_coverage=0.2,
            relative_internal_coverage=0.2,
            terminal_enrichment=-0.1,
            similarity_to_strongest=0.95,
        )

        self.assertTrue(
            experiment.selected_by_threshold(candidate, 0.2, 0.95, 0.25, -0.05)
        )
        self.assertFalse(
            experiment.selected_by_threshold(candidate, 0.1, 0.95, 0.25, -0.05)
        )
        self.assertFalse(
            experiment.selected_by_threshold(candidate, 0.2, 0.98, 0.25, -0.05)
        )
        self.assertFalse(
            experiment.selected_by_threshold(candidate, 0.2, 0.95, 0.05, -0.05)
        )
        self.assertFalse(
            experiment.selected_by_threshold(candidate, 0.2, 0.95, 0.25, -0.15)
        )

    def test_related_reference_has_spaced_reproducible_snps(self) -> None:
        reference = "A" * experiment.REFERENCE_LENGTH

        first, first_positions = experiment.rare_reference(reference, 7)
        second, second_positions = experiment.rare_reference(reference, 7)

        self.assertEqual(first, second)
        self.assertEqual(first_positions, second_positions)
        self.assertEqual(len(first_positions), experiment.SNP_COUNT)
        self.assertEqual(
            sum(left != right for left, right in zip(reference, first, strict=True)),
            experiment.SNP_COUNT,
        )


if __name__ == "__main__":
    unittest.main()
