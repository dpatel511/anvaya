import csv
import tempfile
import unittest
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.dead_end_attribution import write_dead_end_attribution_report
from anvaya.reads import Read
from anvaya.unitig_graph import build_compacted_unitig_graph


class DeadEndAttributionTests(unittest.TestCase):
    def test_linear_read_ends_are_attributed_to_read_boundaries(self) -> None:
        reads = [Read("read-1", "AACTGGA"), Read("read-2", "AACTGGA")]
        graph = build_bidirected_dbg([read.sequence for read in reads], 3, 2)
        unitigs = build_compacted_unitig_graph(graph)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "dead-ends.tsv"
            summary = write_dead_end_attribution_report(
                unitigs, reads, path, min_count=2
            )
            with path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(summary.dead_ends, 2)
        self.assertEqual(summary.read_boundaries, 2)
        self.assertEqual({row["cause"] for row in rows}, {"read_boundary"})

    def test_filtered_extension_retains_missing_quality_provenance(self) -> None:
        reads = [
            Read("main-1", "AACTG"),
            Read("main-2", "AACTG"),
            Read("extension", "CTGA"),
        ]
        graph = build_bidirected_dbg([read.sequence for read in reads], 3, 2)
        unitigs = build_compacted_unitig_graph(graph)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "dead-ends.tsv"
            summary = write_dead_end_attribution_report(
                unitigs, reads, path, min_count=2
            )
            with path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertGreaterEqual(summary.filtered_unique, 1)
        self.assertGreaterEqual(summary.missing_quality, 1)
        self.assertIn(
            "filtered_unique_missing_quality",
            {row["cause"] for row in rows},
        )

    def test_ambiguous_raw_extension_is_ignored(self) -> None:
        reads = [
            Read("main-1", "AACTG"),
            Read("main-2", "AACTG"),
            Read("ambiguous", "CTGN"),
        ]
        graph = build_bidirected_dbg([read.sequence for read in reads], 3, 2)
        unitigs = build_compacted_unitig_graph(graph)

        with tempfile.TemporaryDirectory() as directory:
            summary = write_dead_end_attribution_report(
                unitigs,
                reads,
                Path(directory) / "dead-ends.tsv",
                min_count=2,
            )

        self.assertEqual(summary.filtered_unique, 0)


if __name__ == "__main__":
    unittest.main()
