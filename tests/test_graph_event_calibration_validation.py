import importlib.util
import sys
import unittest
from collections import Counter
from pathlib import Path


def _load_experiment():
    path = (
        Path(__file__).parents[1]
        / "experiments"
        / "15_graph_event_calibration_validation.py"
    )
    specification = importlib.util.spec_from_file_location(
        "graph_event_calibration_validation",
        path,
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


experiment = _load_experiment()


class GraphEventCalibrationValidationTests(unittest.TestCase):
    def test_conditions_expand_truth_processes_and_hold_coverage_constant(
        self,
    ) -> None:
        reference_length = 2_000
        total_coverage = 20
        conditions = experiment._conditions(
            42,
            reference_length,
            total_coverage,
        )

        self.assertEqual(
            Counter(condition.truth for condition in conditions),
            {"damage": 3, "sequencing_error": 5, "variation": 8},
        )
        scenarios = {condition.scenario for condition in conditions}
        self.assertIn("independent_error_0.001", scenarios)
        self.assertIn("independent_terminal_error", scenarios)
        self.assertIn("rare_damage_compatible_10x", scenarios)
        self.assertIn("rare_ordinary_10x", scenarios)

        expected_reads = total_coverage * reference_length / experiment.previous.READ_LENGTH
        for condition in conditions:
            if condition.truth != "variation":
                continue
            self.assertLessEqual(abs(len(condition.reads) - expected_reads), 2)


if __name__ == "__main__":
    unittest.main()
