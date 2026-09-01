import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

from anvaya.reads import load_reads


def _load_module():
    path = (
        Path(__file__).resolve().parents[1]
        / "experiments"
        / "20_restore_unpolished_overlap.py"
    )
    spec = importlib.util.spec_from_file_location("restore_unpolished_overlap", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class RestoreUnpolishedOverlapTests(unittest.TestCase):
    def test_restores_reported_bases_without_changing_topology(self) -> None:
        module = _load_module()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            polished = root / "polished.fasta"
            report = root / "corrections.tsv"
            output = root / "unpolished.fasta"
            polished.write_text(
                ">contig_1\nCCGT\n>contig_2\nAAGG\n",
                encoding="utf-8",
            )
            report.write_text(
                "cluster\tposition\toriginal_base\tcorrected_base\n"
                "1\t0\tT\tC\n"
                "2\t3\tA\tG\n",
                encoding="utf-8",
            )

            restored = module.restore_unpolished(polished, report, output)
            sequences = [read.sequence for read in load_reads(output)]

        self.assertEqual(restored, 2)
        self.assertEqual(sequences, ["TCGT", "AAGA"])


if __name__ == "__main__":
    unittest.main()
