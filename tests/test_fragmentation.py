import csv
import tempfile
import unittest
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.fragmentation import write_fragmentation_report
from anvaya.unitig_graph import build_compacted_unitig_graph


class FragmentationReportTests(unittest.TestCase):
    def test_linear_unitig_has_two_dead_ends(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA"], 3)
        compacted = build_compacted_unitig_graph(graph)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fragmentation.tsv"
            summary = write_fragmentation_report(compacted, path)
            with path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(summary.ends, 2)
        self.assertEqual(summary.dead_ends, 2)
        self.assertEqual(summary.unique_continuations, 0)
        self.assertEqual(summary.ambiguous_branches, 0)
        self.assertEqual(summary.isolated_unitigs, 1)
        self.assertEqual(summary.one_sided_unitigs, 0)
        self.assertEqual(summary.connected_unitigs, 0)
        self.assertEqual([row["side"] for row in rows], ["left", "right"])
        self.assertEqual({row["category"] for row in rows}, {"dead_end"})
        self.assertEqual({row["unitig_topology"] for row in rows}, {"isolated"})

    def test_branching_graph_reports_ambiguous_boundaries(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA", "AACTAGA"], 3)
        compacted = build_compacted_unitig_graph(graph)

        with tempfile.TemporaryDirectory() as directory:
            summary = write_fragmentation_report(
                compacted, Path(directory) / "fragmentation.tsv"
            )

        self.assertGreater(summary.ambiguous_branches, 0)
        self.assertEqual(
            summary.ends,
            summary.dead_ends
            + summary.unique_continuations
            + summary.ambiguous_branches,
        )
        self.assertEqual(
            compacted.unitig_count,
            summary.isolated_unitigs
            + summary.one_sided_unitigs
            + summary.connected_unitigs,
        )


if __name__ == "__main__":
    unittest.main()
