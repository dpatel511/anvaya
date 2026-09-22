import importlib.util
import sys
import unittest
from pathlib import Path


root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / "experiments"))
path = root / "experiments/variable_length_validation.py"
spec = importlib.util.spec_from_file_location("variable_length_validation_test", path)
validation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validation)


class VariableLengthValidationTests(unittest.TestCase):
    def test_default_plan_is_development(self):
        self.assertEqual(
            validation.resolve_seed_plan(),
            (validation.DEVELOPMENT_SEEDS, "development"),
        )

    def test_legacy_held_out_seeds_are_consumed_diagnosis(self):
        self.assertEqual(
            validation.resolve_seed_plan(held_out=True),
            (validation.HELD_OUT_SEEDS, "consumed_validation_diagnosis"),
        )

    def test_consumed_seed_cannot_be_claimed_as_fresh_validation(self):
        for seed in (96012, 197101, 197102, 197103, 197104, 20261401, 20261402):
            with self.subTest(seed=seed), self.assertRaisesRegex(
                ValueError, "consumed validation seeds"
            ):
                validation.resolve_seed_plan([seed], seed_role="validation")

    def test_fresh_validation_requires_explicit_seeds(self):
        with self.assertRaisesRegex(ValueError, "requires explicit"):
            validation.resolve_seed_plan(seed_role="validation")

    def test_fresh_validation_requires_preregistered_seeds(self):
        with self.assertRaisesRegex(ValueError, "not preregistered"):
            validation.resolve_seed_plan([97001], seed_role="validation")

    def test_error_only_role_is_explicit(self):
        self.assertEqual(
            validation.resolve_seed_plan(
                [97001], error_only=True, seed_role="development"
            ),
            ((97001,), "development_error_only"),
        )

    def test_fixture_can_return_private_damage_and_error_truth(self):
        references, reads, origins, events, errors = validation.fixture(
            96012, "strains", .25, .002, include_events=True,
        )
        self.assertEqual(len(reads), len(origins))
        self.assertGreater(len(events), 0)
        self.assertGreater(len(errors), 0)
        for read_index, position, latent, observed in events:
            self.assertEqual(reads[read_index].sequence[position], observed)
            self.assertEqual(reads[read_index].qualities[position], 35)
            self.assertIn((latent, observed), (("C", "T"), ("G", "A")))
        self.assertEqual(len(references), 2)

    def test_fixture_coverage_controls_read_count(self):
        low = validation.fixture(199001, "unique", .25, .002, coverage=5)
        high = validation.fixture(199001, "unique", .25, .002, coverage=20)
        self.assertEqual(len(high[1]), 4 * len(low[1]))
        with self.assertRaisesRegex(ValueError, "positive"):
            validation.fixture(199001, "unique", .25, .002, coverage=0)

    def test_case_identities_are_deterministic_and_content_sensitive(self):
        first = validation.evaluate(199001, *validation.SCENARIOS[0], coverage=2)
        second = validation.evaluate(199001, *validation.SCENARIOS[0], coverage=2)
        for field in (
            "input_sha256", "reference_sha256", "truth_sha256",
            "damage_profile_sha256",
        ):
            self.assertEqual(first[field], second[field])
            self.assertEqual(len(first[field]), 64)
        changed = validation.evaluate(199002, *validation.SCENARIOS[0], coverage=2)
        self.assertNotEqual(first["input_sha256"], changed["input_sha256"])
        self.assertNotEqual(first["reference_sha256"], changed["reference_sha256"])

    def test_manifest_declares_complete_matrix_and_source_identity(self):
        manifest = validation.experiment_manifest(
            (1, 2), "development", validation.SCENARIOS[:2], (2.0, 20.0),
        )
        self.assertEqual(manifest["schema_version"], 2)
        self.assertEqual(len(manifest["expected_cases"]), 8)
        self.assertEqual(manifest["modes"], ["exact", "damage_aware"])
        self.assertEqual(manifest["stages"], ["layout", "consensus"])
        self.assertEqual(len(manifest["generator_sha256"]), 64)
        self.assertEqual(len(manifest["evaluator_sha256"]), 64)
        self.assertEqual(manifest["generator_path"], "experiments/variable_length_validation.py")
        self.assertEqual(manifest["evaluator_path"], "experiments/controlled_assembly.py")

if __name__ == "__main__":
    unittest.main()
