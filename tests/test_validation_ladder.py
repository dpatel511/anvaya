import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


def _load_module():
    path = Path(__file__).resolve().parents[1] / "experiments" / "18_validation_ladder.py"
    spec = importlib.util.spec_from_file_location("validation_ladder", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ValidationLadderTests(unittest.TestCase):
    def test_n50_uses_half_of_total_bases(self) -> None:
        module = _load_module()
        self.assertEqual(module._n50([10, 9, 1]), 10)
        self.assertEqual(module._n50([]), 0)

    def test_runner_writes_truth_labelled_outputs(self) -> None:
        module = _load_module()
        config = json.loads(
            (Path(__file__).resolve().parents[1] / "experiments" / "benchmark_ladder.json").read_text()
        )
        config["reference_length"] = 240
        config["read_length"] = 60
        config["seeds"] = [2401]
        config["scenarios"] = ["clean", "damage_error"]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            output = root / "results"
            records = module.run_benchmark(config_path, output)
            self.assertEqual(len(records), 2)
            self.assertTrue((output / "summary.tsv").is_file())
            self.assertTrue((output / "summary.json").is_file())
            self.assertTrue((output / "manifest.json").is_file())
            self.assertEqual({record.scenario for record in records}, {"clean", "damage_error"})
            self.assertGreaterEqual(records[1].damage_events, 0)


if __name__ == "__main__":
    unittest.main()
