import importlib.util
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]


def _runner():
    spec = importlib.util.spec_from_file_location(
        "overlap_regression", ROOT / "experiments/overlap_regression.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class OverlapRegressionTests(unittest.TestCase):
    def test_runner_compares_real_outputs_and_rejects_changed_hash(self) -> None:
        module = _runner()
        config = json.loads((ROOT / "experiments/overlap_regression.json").read_text())
        config.update(reference_length=400, read_length=80, seeds=[2401], scenarios=["clean", "damage_error"])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            baseline = module.run_regression(config_path, root / "baseline")
            self.assertEqual(len(baseline["records"]), 2)
            for record in baseline["records"]:
                self.assertEqual(len(record["outputs"]), 4)
                self.assertGreater(record["outputs"]["contigs.fasta"]["bases"], 0)
            self.assertEqual(baseline["purpose"], "structural_regression_not_biological_acceptance")
            manifest = root / "baseline/manifest.json"
            argv = ["overlap_regression.py", "--config", str(config_path), "--compare", str(manifest), "--output-dir"]
            with patch.object(sys, "argv", argv + [str(root / "recheck")]), redirect_stdout(io.StringIO()):
                self.assertEqual(module.main(), 0)
            baseline["records"][0]["outputs"]["contigs.fasta"]["sha256"] = "incorrect"
            manifest.write_text(json.dumps(baseline), encoding="utf-8")
            with patch.object(sys, "argv", argv + [str(root / "mismatch")]), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    module.main()
                self.assertEqual(error.exception.code, 1)
            with self.assertRaises(FileExistsError):
                module.run_regression(config_path, root / "baseline")

    def test_unknown_scenario_rejected_before_output_creation(self) -> None:
        module = _runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.json"
            config.write_text(json.dumps({"seeds": [1], "scenarios": ["unknown"]}), encoding="utf-8")
            with self.assertRaises(ValueError):
                module.run_regression(config, root / "out")
            self.assertFalse((root / "out").exists())


if __name__ == "__main__":
    unittest.main()
