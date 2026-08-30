import io
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.cli import main


class CliTests(unittest.TestCase):
    def test_damage_consensus_runs_before_graph_and_writes_report(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads_path = root / "reads.fastq"
            output_path = root / "contigs.fasta"
            report_path = root / "damage-consensus.tsv"
            truth = "CCGATCGTACGTTAGCTACGATCGTACGTA"
            records = [("damaged", "T" + truth[1:])] + [
                (f"truth-{index}", "GGGGG" + truth) for index in range(3)
            ]
            reads_path.write_text(
                "".join(
                    f"@{name}\n{sequence}\n+\n{'D' * len(sequence)}\n"
                    for name, sequence in records
                ),
                encoding="utf-8",
            )

            exit_code = main(
                [
                    "assemble",
                    "--input",
                    str(reads_path),
                    "--k",
                    "7",
                    "--min-count",
                    "2",
                    "--orientation-aware",
                    "--end-window",
                    "5",
                    "--damage-consensus",
                    "--damage-consensus-report",
                    str(report_path),
                    "--output",
                    str(output_path),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue(output_path.is_file())
            self.assertIn("correct", report_path.read_text(encoding="utf-8"))

    def test_dead_end_and_correction_reports_are_report_only(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads_path = root / "reads.fastq"
            baseline_path = root / "baseline.fasta"
            reported_path = root / "reported.fasta"
            dead_end_path = root / "dead-ends.tsv"
            correction_path = root / "corrections.tsv"
            reads_path.write_text(
                "@truth-1\nGGGAACCG\n+\nDDDDDDDD\n"
                "@truth-2\nGGGAACCG\n+\nDDDDDDDD\n"
                "@error\nGGGAATCG\n+\nDDDDD)DD\n",
                encoding="utf-8",
            )

            baseline_exit = main(
                [
                    "assemble",
                    "--input",
                    str(reads_path),
                    "--k",
                    "3",
                    "--min-count",
                    "2",
                    "--orientation-aware",
                    "--output",
                    str(baseline_path),
                ]
            )
            reported_exit = main(
                [
                    "assemble",
                    "--input",
                    str(reads_path),
                    "--k",
                    "3",
                    "--min-count",
                    "2",
                    "--orientation-aware",
                    "--dead-end-report",
                    str(dead_end_path),
                    "--correction-report",
                    str(correction_path),
                    "--output",
                    str(reported_path),
                ]
            )

            self.assertEqual(baseline_exit, 0)
            self.assertEqual(reported_exit, 0)
            self.assertEqual(baseline_path.read_bytes(), reported_path.read_bytes())
            self.assertIn("cause", dead_end_path.read_text(encoding="utf-8"))
            self.assertIn("would_correct", correction_path.read_text(encoding="utf-8"))

    def test_fragmentation_report_is_report_only(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads_path = root / "reads.fasta"
            baseline_path = root / "baseline.fasta"
            reported_path = root / "reported.fasta"
            report_path = root / "fragmentation.tsv"
            reads_path.write_text(
                ">read-1\nAACTGGA\n>read-2\nAACTAGA\n", encoding="utf-8"
            )

            baseline_exit = main(
                [
                    "assemble",
                    "--input",
                    str(reads_path),
                    "--k",
                    "3",
                    "--orientation-aware",
                    "--output",
                    str(baseline_path),
                ]
            )
            reported_exit = main(
                [
                    "assemble",
                    "--input",
                    str(reads_path),
                    "--k",
                    "3",
                    "--orientation-aware",
                    "--fragmentation-report",
                    str(report_path),
                    "--output",
                    str(reported_path),
                ]
            )

            self.assertEqual(baseline_exit, 0)
            self.assertEqual(reported_exit, 0)
            self.assertEqual(baseline_path.read_bytes(), reported_path.read_bytes())
            self.assertIn("category", report_path.read_text(encoding="utf-8"))

    def test_assemble_command_writes_fasta_summary_and_progress(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads_path = root / "reads.fasta"
            output_path = root / "unitigs.fasta"
            reads_path.write_text(">read-1\nATGCA\n", encoding="utf-8")
            stdout = io.StringIO()
            stderr = io.StringIO()

            with redirect_stdout(stdout), redirect_stderr(stderr):
                exit_code = main(
                    [
                        "assemble",
                        "--input",
                        str(reads_path),
                        "--k",
                        "3",
                        "--output",
                        str(output_path),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertEqual(
                output_path.read_text(encoding="utf-8"),
                ">unitig_1 length=5\nATGCA\n",
            )
            self.assertIn("reads=1", stdout.getvalue())
            self.assertIn("unitigs=1", stdout.getvalue())
            self.assertIn("[anvaya] Loading reads", stderr.getvalue())
            self.assertIn("[anvaya] Building directed de Bruijn graph", stderr.getvalue())
            self.assertIn("[anvaya] Extracting unitigs", stderr.getvalue())
            self.assertIn("[anvaya] Writing unitigs", stderr.getvalue())
            self.assertIn("[anvaya] Completed", stderr.getvalue())

    def test_module_entry_point_displays_help(self) -> None:
        result = subprocess.run(
            [sys.executable, "-m", "anvaya", "--help"],
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn("assemble", result.stdout)


if __name__ == "__main__":
    unittest.main()
