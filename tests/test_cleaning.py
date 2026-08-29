import io
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from anvaya.bidirected import (
    build_bidirected_dbg,
    extract_bidirected_unitigs,
    summarize_bidirected_graph,
)
from anvaya.cleaning import remove_damage_aware_tips, remove_weak_tips
from anvaya.cli import main
from anvaya.sequences import canonical_sequence


class TipCleaningTests(unittest.TestCase):
    def test_removes_short_weak_tip(self) -> None:
        main = "AAGCCCAAT"
        tip = "AAGCCTAAA"
        graph = build_bidirected_dbg([main] * 10 + [tip], 5)
        initial_edges = graph.edge_count
        initial_observations = graph.observations

        summary = remove_weak_tips(graph)
        unitigs = extract_bidirected_unitigs(graph)
        graph_summary = summarize_bidirected_graph(graph, unitigs)

        self.assertEqual(summary.tips_removed, 1)
        self.assertGreater(summary.edges_removed, 0)
        self.assertEqual(
            graph_summary.edges,
            initial_edges - summary.edges_removed,
        )
        self.assertEqual(
            graph_summary.observations,
            initial_observations - summary.observations_removed,
        )
        self.assertEqual(len(unitigs), 1)
        self.assertEqual(canonical_sequence(unitigs[0]), canonical_sequence(main))

        repeated = remove_weak_tips(graph)
        self.assertEqual(repeated.tips_removed, 0)
        self.assertEqual(repeated.edges_removed, 0)

    def test_preserves_strong_alternative(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAT"] * 10 + ["AAGCCTAAA"] * 5,
            5,
        )

        summary = remove_weak_tips(graph)

        self.assertEqual(summary.tips_removed, 0)

    def test_preserves_isolated_linear_path(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA"] * 10, 3)

        summary = remove_weak_tips(graph)

        self.assertEqual(summary.tips_removed, 0)
        self.assertEqual(len(extract_bidirected_unitigs(graph)), 1)

    def test_preserves_tip_longer_than_limit(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAT"] * 10 + ["AAGCCTAAA"],
            5,
        )

        summary = remove_weak_tips(graph, max_tip_length=5)

        self.assertEqual(summary.tips_removed, 0)

    def test_damage_aware_cleaning_protects_damage_like_tip(self) -> None:
        graph = build_bidirected_dbg(
            ["ACATACACG"] * 10 + ["ACATACACA"],
            5,
            end_window=3,
        )

        summary = remove_damage_aware_tips(graph)

        self.assertEqual(summary.tips_removed, 0)
        self.assertEqual(summary.tips_protected, 1)
        self.assertEqual(summary.damage_like_protected, 1)

    def test_damage_aware_cleaning_removes_only_error_like_tip(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAT"] * 10 + ["AAGCCTAAA"],
            5,
            end_window=3,
        )

        with (
            patch(
                "anvaya.tip_matching.match_tip_to_backbone",
                return_value=object(),
            ),
            patch(
                "anvaya.tip_classification.classify_tip_match",
                return_value=SimpleNamespace(
                    label="error-like",
                    evidence=SimpleNamespace(tip_strand_balance=0.0),
                ),
            ),
        ):
            summary = remove_damage_aware_tips(graph)

        self.assertEqual(summary.tips_removed, 1)
        self.assertEqual(summary.tips_protected, 0)
        self.assertEqual(summary.rounds, 1)

    def test_damage_aware_cleaning_protects_bidirectional_error_like_tip(
        self,
    ) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAT"] * 10 + ["AAGCCTAAA"],
            5,
            end_window=3,
        )

        with (
            patch(
                "anvaya.tip_matching.match_tip_to_backbone",
                return_value=object(),
            ),
            patch(
                "anvaya.tip_classification.classify_tip_match",
                return_value=SimpleNamespace(
                    label="error-like",
                    evidence=SimpleNamespace(tip_strand_balance=1.0),
                ),
            ),
        ):
            summary = remove_damage_aware_tips(graph)

        self.assertEqual(summary.tips_removed, 0)
        self.assertEqual(summary.tips_protected, 1)
        self.assertEqual(summary.bidirectional_protected, 1)

    def test_damage_aware_cleaning_requires_end_evidence(self) -> None:
        graph = build_bidirected_dbg(["AAGCCCAAT"] * 10 + ["AAGCCTAAA"], 5)

        with self.assertRaisesRegex(ValueError, "requires read-end evidence"):
            remove_damage_aware_tips(graph)


class TipCleaningCliTests(unittest.TestCase):
    def test_tip_cleaning_requires_orientation_aware_graph(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            reads.write_text(">read\nAAGCCCAAT\n", encoding="utf-8")

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(
                    [
                        "assemble",
                        "-i",
                        str(reads),
                        "--k",
                        "5",
                        "--clean-tips",
                        "-o",
                        str(root / "unitigs.fasta"),
                    ]
                )

    def test_cli_reports_tip_cleaning(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            output = root / "unitigs.fasta"
            records = ["AAGCCCAAT"] * 10 + ["AAGCCTAAA"]
            reads.write_text(
                "".join(f">read_{index}\n{read}\n" for index, read in enumerate(records)),
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                exit_code = main(
                    [
                        "assemble",
                        "-i",
                        str(reads),
                        "--k",
                        "5",
                        "--orientation-aware",
                        "--clean-tips",
                        "-o",
                        str(output),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertIn("tip_cleaning=true", stdout.getvalue())
            self.assertIn("tips_removed=1", stdout.getvalue())

    def test_cli_reports_damage_aware_tip_protection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            output = root / "unitigs.fasta"
            records = ["ACATACACG"] * 10 + ["ACATACACA"]
            reads.write_text(
                "".join(
                    f">read_{index}\n{read}\n"
                    for index, read in enumerate(records)
                ),
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                exit_code = main(
                    [
                        "assemble",
                        "-i",
                        str(reads),
                        "--k",
                        "5",
                        "--orientation-aware",
                        "--end-window",
                        "3",
                        "--damage-aware-clean-tips",
                        "-o",
                        str(output),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertIn("damage_aware_tip_cleaning=true", stdout.getvalue())
            self.assertIn("damage_like_tips_protected=1", stdout.getvalue())

    def test_damage_aware_cleaning_requires_end_window(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            reads.write_text(">read\nAAGCCCAAT\n", encoding="utf-8")

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(
                    [
                        "assemble",
                        "-i",
                        str(reads),
                        "--k",
                        "5",
                        "--orientation-aware",
                        "--damage-aware-clean-tips",
                        "-o",
                        str(root / "unitigs.fasta"),
                    ]
                )


if __name__ == "__main__":
    unittest.main()
