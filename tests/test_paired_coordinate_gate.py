import copy
import importlib.util
import unittest
from pathlib import Path


path = Path(__file__).resolve().parents[1] / "experiments/paired_coordinate_gate.py"
spec = importlib.util.spec_from_file_location("paired_coordinate_gate", path)
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)
HASH = "a" * 64


def metrics(states, n50=1):
    discordant = sum(state["status"] == "discordant" for state in states)
    evaluated = sum(state["status"].startswith("concordant_") for state in states)
    mismatches = sum(state["status"] == "concordant_incorrect" for state in states)
    return {
        "nonduplicated_accuracy_coordinate_states": states,
        "diagnostic_mosaic_bases": 0,
        "incompatible_candidate_sequences_bases": 0,
        "accuracy_unresolved_bases": 0,
        "nonduplicated_accuracy_discordant_reference_bases": discordant,
        "nonduplicated_accuracy_evaluated_reference_bases": evaluated,
        "nonduplicated_accuracy_mismatches": mismatches,
        "n50": n50,
        "longest": n50,
    }


def state(position, status, base):
    return {
        "reference": 0, "position_0based": position, "reference_base": "A",
        "observed_bases": [base], "observation_count": 1,
        "origin_resolution": ["unique_origin"], "status": status,
    }


def case(seed=1, scenario="damage_unique", coverage=20, damage=.2, n50=100):
    score = metrics([state(0, "concordant_correct", "A")], n50=n50)
    result = {
        "seed": seed, "scenario": scenario, "target_coverage": coverage,
        "panel": "unique", "damage": damage, "sequencing_error_rate": 0,
        "supplied_profile_scale": 1, "reads": 1,
        "minimum_read_length": 100, "maximum_read_length": 100,
        "simulated_damage_events": int(damage > 0),
        "simulated_sequencing_errors": 0,
        "results": {
            mode: {stage: copy.deepcopy(score) for stage in gate.STAGES}
            for mode in gate.MODES
        },
    }
    result.update({field: HASH for field in gate.CASE_IDENTITY_FIELDS})
    return result


def summary(cases, role="development", assembly_hash=HASH):
    expected = [
        {"seed": item["seed"], "scenario": item["scenario"],
         "target_coverage": item["target_coverage"]}
        for item in cases
    ]
    return {
        "seed_role": role,
        "experiment_manifest": {
            "schema_version": gate.SCHEMA_VERSION,
            "generator_sha256": HASH,
            "evaluator_sha256": HASH,
            "generator_path": "experiments/generator.py",
            "evaluator_path": "experiments/evaluator.py",
            "modes": list(gate.MODES), "stages": list(gate.STAGES),
            "expected_cases": expected,
        },
        "assembly_source_sha256": assembly_hash,
        "source_sha256": {
            "experiments/generator.py": HASH,
            "experiments/evaluator.py": HASH,
        },
        "cases": cases,
    }


def passing_pair(seed=1):
    clean = case(seed, "clean_unique", damage=0)
    damaged = case(seed, "damage_unique", damage=.2)
    before = summary([clean, damaged])
    after = copy.deepcopy(before)
    after["assembly_source_sha256"] = "b" * 64
    for item in after["cases"]:
        for mode in gate.MODES:
            item["results"][mode]["layout"]["n50"] = 150
            item["results"][mode]["consensus"]["n50"] = 150
    return before, after


class PairedCoordinateGateTests(unittest.TestCase):
    def test_reports_new_conflict_instead_of_hiding_previous_error(self):
        before = metrics([state(0, "concordant_incorrect", "C")])
        discordant = state(0, "discordant", "A")
        discordant["observed_bases"] = ["A", "C"]
        discordant["observation_count"] = 2
        result = gate.compare_stage(before, metrics([discordant]))
        self.assertEqual(result["transitions"], {"concordant_incorrect->discordant": 1})
        self.assertIn("new_discordant_coordinate:0:0", result["failures"])

    def test_reports_lost_correct_and_newly_recovered_coordinates(self):
        result = gate.compare_stage(
            metrics([state(0, "concordant_correct", "A")]),
            metrics([state(1, "concordant_correct", "A")]),
        )
        self.assertIn("previously_correct_coordinate_regressed:0:0", result["failures"])
        self.assertEqual(result["correct_coordinate_gain"], 1)

    def test_rejects_duplicate_coordinate_states(self):
        duplicated = metrics([
            state(0, "concordant_correct", "A"),
            state(0, "concordant_correct", "A"),
        ])
        with self.assertRaisesRegex(ValueError, "duplicate coordinate state"):
            gate.compare_stage(duplicated, metrics([]))

    def test_rejects_reference_base_mismatch(self):
        old = state(0, "concordant_correct", "A")
        new = state(0, "concordant_correct", "A")
        new["reference_base"] = "C"
        with self.assertRaisesRegex(ValueError, "reference base differs"):
            gate.compare_stage(metrics([old]), metrics([new]))

    def test_fully_declared_candidate_passes(self):
        before, after = passing_pair()
        result = gate.compare_summaries(before, after)
        self.assertTrue(result["passed"])
        self.assertEqual(result["case_combinations"], 2)
        self.assertEqual(result["distinct_reference_seeds"], [1])

    def test_rejects_no_op_summary(self):
        before, _ = passing_pair()
        after = copy.deepcopy(before)
        after["assembly_source_sha256"] = "b" * 64
        result = gate.compare_summaries(before, after)
        self.assertFalse(result["passed"])
        self.assertIn("candidate_is_no_op", result["global_failures"])

    def test_error_only_control_can_be_safe_no_op(self):
        control = case(scenario="error_only_unique", damage=0)
        before = summary([control], role="development_error_only")
        after = copy.deepcopy(before)
        after["assembly_source_sha256"] = "b" * 64
        self.assertTrue(gate.compare_summaries(before, after)["passed"])

    def test_rejects_exact_only_gain_as_damage_benefit(self):
        before, after = passing_pair()
        damaged = after["cases"][1]
        for stage in gate.STAGES:
            damaged["results"]["damage_aware"][stage]["n50"] = 100
        result = gate.compare_summaries(before, after)
        self.assertIn("no_damage_aware_contiguity_benefit:1", result["global_failures"])

    def test_requires_fifty_percent_clean_gain_in_both_modes(self):
        before, after = passing_pair()
        after["cases"][0]["results"]["exact"]["layout"]["n50"] = 149
        result = gate.compare_summaries(before, after)
        self.assertIn("insufficient_clean_20x_n50_gain:1", result["global_failures"])

    def test_requires_clean_case_for_every_seed(self):
        damaged = case()
        before = summary([damaged])
        after = copy.deepcopy(before)
        after["assembly_source_sha256"] = "b" * 64
        for stage in gate.STAGES:
            after["cases"][0]["results"]["damage_aware"][stage]["n50"] = 150
        result = gate.compare_summaries(before, after)
        self.assertIn("insufficient_clean_20x_n50_gain:1", result["global_failures"])

    def test_rejects_candidate_over_twice_baseline_assembly_time(self):
        before, after = passing_pair()
        for case_before, case_after in zip(before["cases"], after["cases"]):
            for mode in gate.MODES:
                case_before["results"][mode]["assembly_seconds"] = 1.0
                case_after["results"][mode]["assembly_seconds"] = 2.01
        result = gate.compare_summaries(before, after)
        self.assertFalse(result["passed"])
        self.assertTrue(all(
            "assembly_time_exceeded_2x" in row["failures"]
            for row in result["rows"]
        ))

    def test_rejects_missing_or_changed_content_identity(self):
        before, after = passing_pair()
        del after["cases"][0]["input_sha256"]
        with self.assertRaisesRegex(ValueError, "input_sha256"):
            gate.compare_summaries(before, after)
        before, after = passing_pair()
        after["cases"][0]["input_sha256"] = "b" * 64
        with self.assertRaisesRegex(ValueError, "differs in input_sha256"):
            gate.compare_summaries(before, after)

    def test_rejects_missing_or_changed_manifest_identity(self):
        before, after = passing_pair()
        del after["experiment_manifest"]
        with self.assertRaisesRegex(ValueError, "missing experiment manifest"):
            gate.compare_summaries(before, after)
        before, after = passing_pair()
        after["experiment_manifest"]["evaluator_sha256"] = "b" * 64
        after["source_sha256"]["experiments/evaluator.py"] = "b" * 64
        with self.assertRaisesRegex(ValueError, "manifests differ"):
            gate.compare_summaries(before, after)

    def test_rejects_missing_source_manifest_or_identical_assembly(self):
        before, after = passing_pair()
        del after["source_sha256"]
        with self.assertRaisesRegex(ValueError, "source hash manifest"):
            gate.compare_summaries(before, after)
        before, after = passing_pair()
        after["assembly_source_sha256"] = before["assembly_source_sha256"]
        with self.assertRaisesRegex(ValueError, "assembly source hashes are identical"):
            gate.compare_summaries(before, after)

    def test_rejects_case_jointly_omitted_from_both_summaries(self):
        before, after = passing_pair()
        before["cases"].pop()
        after["cases"].pop()
        with self.assertRaisesRegex(ValueError, "cases differ from manifest"):
            gate.compare_summaries(before, after)

    def test_rejects_duplicate_manifest_or_observed_cases(self):
        before, after = passing_pair()
        before["experiment_manifest"]["expected_cases"].append(
            copy.deepcopy(before["experiment_manifest"]["expected_cases"][0])
        )
        after["experiment_manifest"] = copy.deepcopy(before["experiment_manifest"])
        with self.assertRaisesRegex(ValueError, "duplicate expected case"):
            gate.compare_summaries(before, after)
        before, after = passing_pair()
        before["cases"].append(copy.deepcopy(before["cases"][0]))
        with self.assertRaisesRegex(ValueError, "duplicate case"):
            gate.compare_summaries(before, after)


if __name__ == "__main__":
    unittest.main()
