import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cli import main
from anvaya.paired_extension import spell_extended_paths
from anvaya.read_threading import (
    audit_read_threads,
    resolve_read_thread_extensions,
)
from anvaya.reads import Read
from anvaya.unitig_graph import build_compacted_unitig_graph


class ReadThreadingTests(unittest.TestCase):
    def test_counts_independent_molecules_crossing_a_branch(self) -> None:
        winner = Read("winner", "AACTGGA", (30,) * 7)
        alternative = Read("alternative", "AACTAGA", (20,) * 7)
        reads = [winner] * 6 + [alternative]
        source = build_bidirected_dbg(
            [read.sequence for read in reads],
            3,
        )
        unitigs = build_compacted_unitig_graph(
            source,
            track_edge_handles=True,
        )

        transitions, summary = audit_read_threads(
            source,
            unitigs,
            ((read,) for read in reads),
            end_window=1,
        )

        self.assertGreaterEqual(summary.junctions, 1)
        self.assertGreaterEqual(summary.supported_transitions, 2)
        self.assertGreaterEqual(summary.resolvable_junctions, 1)
        self.assertGreaterEqual(summary.reciprocal_resolvable_links, 1)
        self.assertEqual(summary.threaded_molecules, 7)
        self.assertTrue(any(item.molecules == 6 for item in transitions))
        self.assertTrue(any(item.mean_quality == 30.0 for item in transitions))
        self.assertTrue(any(item.reciprocal for item in transitions))
        successors = resolve_read_thread_extensions(transitions)
        extended = spell_extended_paths(unitigs, successors)
        self.assertEqual(len(successors), 2)
        self.assertLess(len(extended), unitigs.unitig_count)

    def test_breaks_weaker_link_through_two_copy_unitig(self) -> None:
        winner = Read("winner", "AACTGGA", (30,) * 7)
        alternative = Read("alternative", "AACTAGA", (20,) * 7)
        reads = [winner] * 6 + [alternative]
        source = build_bidirected_dbg([read.sequence for read in reads], 3)
        unitigs = build_compacted_unitig_graph(source, track_edge_handles=True)
        transitions, _ = audit_read_threads(
            source,
            unitigs,
            ((read,) for read in reads),
            end_window=1,
        )
        accepted = [item for item in transitions if item.reciprocal]
        self.assertTrue(accepted)
        first = replace(
            accepted[0],
            source=0,
            target=2,
            molecules=5,
            source_mean_support=10.0,
            target_mean_support=40.0,
        )
        second = replace(
            accepted[0],
            source=2,
            target=4,
            molecules=10,
            source_mean_support=40.0,
            target_mean_support=10.0,
        )

        successors = resolve_read_thread_extensions((first, second))

        self.assertNotIn(0, successors)
        self.assertEqual(successors[2], 4)
        self.assertEqual(successors[5], 3)

    def test_cli_writes_report_without_changing_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            baseline = root / "baseline.fasta"
            audited = root / "audited.fasta"
            extended = root / "extended.fasta"
            report = root / "threads.tsv"
            reads.write_text(
                "".join(
                    [f">winner-{index}\nAACTGGA\n" for index in range(6)]
                    + [">alternative\nAACTAGA\n"]
                ),
                encoding="utf-8",
            )

            common = [
                "assemble",
                "-i",
                str(reads),
                "--k",
                "3",
                "--orientation-aware",
            ]
            self.assertEqual(main(common + ["-o", str(baseline)]), 0)
            self.assertEqual(
                main(
                    common
                    + [
                        "--read-thread-report",
                        str(report),
                        "-o",
                        str(audited),
                    ]
                ),
                0,
            )

            self.assertEqual(baseline.read_bytes(), audited.read_bytes())
            self.assertIn("source_unitig", report.read_text(encoding="utf-8"))
            self.assertEqual(
                main(
                    common
                    + [
                        "--end-window",
                        "1",
                        "--read-thread-extension",
                        "-o",
                        str(extended),
                    ]
                ),
                0,
            )
            self.assertNotEqual(baseline.read_bytes(), extended.read_bytes())


if __name__ == "__main__":
    unittest.main()
