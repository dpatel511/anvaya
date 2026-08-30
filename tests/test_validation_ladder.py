import csv
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

    def test_oracle_supplements_only_observed_truth_kmers(self) -> None:
        module = _load_module()
        reads = [
            module.SimulatedRead("truth", "AAC", "AAC", 0, 3),
            module.SimulatedRead("error", "AAT", "AAC", 0, 3),
        ]

        supplements, rescued = module._oracle_supplements(
            reads, ("AACCG",), 3, 2
        )

        self.assertEqual(supplements, ["AAC"])
        self.assertEqual(rescued, 1)

    def test_restricted_oracle_supplements_only_selected_kmers(self) -> None:
        module = _load_module()
        reads = [
            module.SimulatedRead("first", "AAC", "AAC", 0, 3),
            module.SimulatedRead("second", "AAT", "AAT", 0, 3),
        ]

        supplements = module._restricted_oracle_supplements(
            reads, {"AAC"}, 3, 2
        )

        self.assertEqual(supplements, ["AAC"])

    def test_unique_oracle_selects_truth_continuation_at_dead_end(self) -> None:
        module = _load_module()
        reads = [
            module.SimulatedRead("first", "AAAC", "AAAC", 0, 4),
            module.SimulatedRead("second", "AAAC", "AAAC", 0, 4),
            module.SimulatedRead("weak", "AACG", "AACG", 0, 4),
        ]
        source = module.build_bidirected_dbg(
            [read.sequence for read in reads], 3, 2
        )
        graph = module.build_compacted_unitig_graph(source)

        selected = module._oracle_unique_kmers(
            graph, reads, ("AAACG",), 2
        )

        self.assertEqual(selected, {"ACG"})

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
            with (output / "summary.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertEqual(
                {row["scenario"] for row in rows},
                {"clean", "damage_error"},
            )
            self.assertEqual({record.scenario for record in records}, {"clean", "damage_error"})
            self.assertGreaterEqual(records[1].damage_events, 0)
            self.assertGreaterEqual(records[0].oracle_unique_rescued_kmers, 0)
            self.assertEqual(
                records[0].oracle_unique_n50_delta,
                records[0].oracle_unique_n50 - records[0].n50,
            )
            self.assertGreaterEqual(
                records[0].oracle_all_truth_rescued_kmers,
                records[0].oracle_unique_rescued_kmers,
            )


if __name__ == "__main__":
    unittest.main()
