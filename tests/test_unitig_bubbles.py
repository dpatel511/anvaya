import csv
import io
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cli import main
from anvaya.unitig_bubbles import (
    find_unitig_bubbles,
    spell_unitig_path,
    write_unitig_bubble_report,
)
from anvaya.unitig_graph import build_compacted_unitig_graph


class UnitigBubbleTests(unittest.TestCase):
    def setUp(self) -> None:
        graph = build_bidirected_dbg(
            ["GCTTGTTCCGGA"] * 10 + ["GCTTATTCCGGA"] * 3,
            4,
            end_window=4,
        )
        self.compacted = build_compacted_unitig_graph(graph)

    def test_scores_simple_alternative_paths_without_graph_changes(self) -> None:
        degrees_before = bytes(self.compacted.out_degrees)
        links_before = self.compacted.out_links.tobytes()

        bubbles = find_unitig_bubbles(self.compacted)

        self.assertEqual(len(bubbles), 1)
        self.assertEqual(len(bubbles[0].paths), 2)
        strongest = bubbles[0].paths[bubbles[0].strongest_path_index]
        weaker = next(path for path in bubbles[0].paths if path is not strongest)
        self.assertEqual(strongest.mean_edge_support, 10.0)
        self.assertEqual(weaker.mean_edge_support, 3.0)
        self.assertAlmostEqual(weaker.relative_coverage, 0.3)
        self.assertEqual(
            weaker.forward_observations
            + weaker.reverse_observations
            + weaker.palindromic_observations,
            weaker.observations,
        )
        self.assertEqual(
            weaker.terminal_observations + weaker.internal_observations,
            weaker.observations,
        )
        self.assertGreaterEqual(weaker.terminal_fraction, 0.0)
        self.assertLessEqual(weaker.terminal_fraction, 1.0)
        if weaker.strand_balance is not None:
            self.assertGreaterEqual(weaker.strand_balance, 0.0)
            self.assertLessEqual(weaker.strand_balance, 1.0)
        self.assertLess(weaker.similarity_to_strongest, 1.0)
        self.assertEqual(bytes(self.compacted.out_degrees), degrees_before)
        self.assertEqual(self.compacted.out_links.tobytes(), links_before)

    def test_respects_search_bound(self) -> None:
        self.assertEqual(
            find_unitig_bubbles(self.compacted, max_nucleotides=1),
            [],
        )

    def test_spells_reported_path(self) -> None:
        bubble = find_unitig_bubbles(self.compacted)[0]

        for path in bubble.paths:
            sequence = spell_unitig_path(self.compacted, path.handles)
            self.assertEqual(len(sequence), path.nucleotide_length)

    def test_rejects_empty_path(self) -> None:
        with self.assertRaisesRegex(ValueError, "must not be empty"):
            spell_unitig_path(self.compacted, ())

    def test_writes_path_scores(self) -> None:
        bubbles = find_unitig_bubbles(self.compacted)

        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "unitig-bubbles.tsv"
            summary = write_unitig_bubble_report(bubbles, report)
            with report.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(summary.bubbles, 1)
        self.assertEqual(summary.paths, 2)
        self.assertEqual(len(rows), 2)
        self.assertEqual(
            {row["strongest_path"] for row in rows},
            {"true", "false"},
        )
        self.assertIn("relative_coverage", rows[0])
        self.assertIn("relative_internal_coverage", rows[0])
        self.assertIn("terminal_fraction", rows[0])
        self.assertIn("strand_balance", rows[0])
        self.assertIn("similarity_to_strongest", rows[0])

    def test_rejects_invalid_bounds(self) -> None:
        with self.assertRaisesRegex(ValueError, "max_unitigs"):
            find_unitig_bubbles(self.compacted, max_unitigs=0)
        with self.assertRaisesRegex(TypeError, "max_nucleotides"):
            find_unitig_bubbles(self.compacted, max_nucleotides=True)


class UnitigBubbleCliTests(unittest.TestCase):
    def test_cli_report_does_not_change_assembly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            baseline = root / "baseline.fasta"
            reported = root / "reported.fasta"
            report = root / "unitig-bubbles.tsv"
            records = ["GCTTGTTCCGGA"] * 10 + ["GCTTATTCCGGA"] * 3
            reads.write_text(
                "".join(
                    f">read_{index}\n{read}\n"
                    for index, read in enumerate(records)
                ),
                encoding="utf-8",
            )

            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "4",
                        "--orientation-aware", "-o", str(baseline),
                    ]
                )
            stdout = io.StringIO()
            stderr = io.StringIO()
            with redirect_stdout(stdout), redirect_stderr(stderr):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "4",
                        "--orientation-aware",
                        "--unitig-bubble-report", str(report),
                        "-o", str(reported),
                    ]
                )

            self.assertTrue(report.exists())
            self.assertIn("unitig_bubbles_detected=1", stdout.getvalue())
            self.assertIn("Scoring compacted-graph bubbles", stderr.getvalue())
            self.assertEqual(reported.read_bytes(), baseline.read_bytes())

    def test_cli_report_requires_orientation_aware_graph(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            reads.write_text(">read\nAACTGGA\n", encoding="utf-8")

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "3",
                        "--unitig-bubble-report", str(root / "report.tsv"),
                        "-o", str(root / "unitigs.fasta"),
                    ]
                )


if __name__ == "__main__":
    unittest.main()
