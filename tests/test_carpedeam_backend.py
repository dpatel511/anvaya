import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from anvaya.carpedeam_backend import run_carpedeam_safe, scan_sequence_input


def write_profile(prefix: Path) -> None:
    header = [f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b]
    for end, transition in (("5p", "C>T"), ("3p", "G>A")):
        values = ["0.2" if label == transition else "0" for label in header]
        Path(f"{prefix}{end}.prof").write_text(
            "\t".join(header) + "\n" + "\t".join(values) + "\n",
            encoding="utf-8",
        )


def write_fastq(path: Path, lengths=(20, 24)) -> None:
    path.write_text(
        "".join(
            f"@read_{index}\n{'ACGT' * (length // 4)}{'A' * (length % 4)}\n+\n"
            f"{'I' * length}\n"
            for index, length in enumerate(lengths)
        ),
        encoding="utf-8",
    )


class CarpeDeamBackendTests(unittest.TestCase):
    def test_streaming_scan_reports_fastq_lengths(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fq"
            write_fastq(path)
            scan = scan_sequence_input(path)
            self.assertEqual(scan.format, "FASTQ")
            self.assertEqual(scan.reads, 2)
            self.assertEqual(scan.bases, 44)
            self.assertEqual(scan.minimum_length, 20)
            self.assertEqual(scan.maximum_length, 24)

    def test_safe_run_pins_inputs_binary_command_and_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fq"
            executable = root / "carpedeam"
            profile = root / "damage_"
            output = root / "contigs.fasta"
            diagnostics = root / "diagnostics.json"
            temporary = root / "tmp"
            write_fastq(reads)
            write_profile(profile)
            executable.write_bytes(b"pinned fake binary")

            def fake_run(command, **kwargs):
                if command[-1] == "-h":
                    return subprocess.CompletedProcess(
                        command, 0,
                        stdout=(
                            "usage: carpedeam ancient_assemble\n"
                            "--ancient-damage --unsafe --min-contig-len "
                            "--num-iter-reads-only --num-iterations "
                            "--min-merge-seq-id --min-cov-safe\n"
                        ),
                        stderr="",
                    )
                self.assertEqual(command[1], "ancient_assemble")
                self.assertEqual(command[command.index("--unsafe") + 1], "0")
                self.assertEqual(command[command.index("--min-contig-len") + 1], "31")
                self.assertEqual(command[command.index("--num-iter-reads-only") + 1], "5")
                self.assertEqual(command[command.index("--num-iterations") + 1], "10")
                self.assertEqual(command[command.index("--min-merge-seq-id") + 1], "0.99")
                self.assertEqual(command[command.index("--min-cov-safe") + 1], "5")
                Path(command[3]).write_text(
                    ">contig_1\n" + "A" * 40 + "\n>contig_2\n" + "C" * 20 + "\n",
                    encoding="utf-8",
                )
                return subprocess.CompletedProcess(command, 0)

            with patch("anvaya.carpedeam_backend.subprocess.run", side_effect=fake_run):
                result = run_carpedeam_safe(
                    reads, profile, output, temporary, diagnostics,
                    executable=executable,
                )

            self.assertEqual(result["status"], "completed")
            self.assertFalse(result["settings"]["unsafe"])
            self.assertEqual(result["input_scan"]["reads"], 2)
            self.assertEqual(result["output_contigs"], 2)
            self.assertEqual(result["output_n50"], 40)
            self.assertEqual(len(result["input_sha256"]), 64)
            self.assertEqual(len(result["executable_sha256"]), 64)
            self.assertEqual(
                json.loads(diagnostics.read_text(encoding="utf-8")), result,
            )

    def test_rejects_short_input_before_starting_backend(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fq"
            profile = root / "damage_"
            executable = root / "carpedeam"
            write_fastq(reads, lengths=(19,))
            write_profile(profile)
            executable.write_bytes(b"fake")
            with patch("anvaya.carpedeam_backend.subprocess.run") as run:
                with self.assertRaisesRegex(ValueError, "at least 20 bp"):
                    run_carpedeam_safe(
                        reads, profile, root / "out.fa", root / "tmp",
                        root / "diagnostics.json", executable=executable,
                    )
                run.assert_not_called()

    def test_failed_backend_preserves_diagnostics_and_logs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fq"
            profile = root / "damage_"
            executable = root / "carpedeam"
            diagnostics = root / "diagnostics.json"
            write_fastq(reads)
            write_profile(profile)
            executable.write_bytes(b"fake")

            def fake_run(command, **kwargs):
                if command[-1] == "-h":
                    return subprocess.CompletedProcess(
                        command, 0,
                        stdout=(
                            "usage: carpedeam ancient_assemble\n"
                            "--ancient-damage --unsafe --min-contig-len "
                            "--num-iter-reads-only --num-iterations "
                            "--min-merge-seq-id --min-cov-safe\n"
                        ),
                        stderr="",
                    )
                kwargs["stderr"].write("backend failed\n")
                return subprocess.CompletedProcess(command, 7)

            with patch("anvaya.carpedeam_backend.subprocess.run", side_effect=fake_run):
                with self.assertRaisesRegex(ValueError, "exit status 7"):
                    run_carpedeam_safe(
                        reads, profile, root / "out.fa", root / "tmp", diagnostics,
                        executable=executable,
                    )
            saved = json.loads(diagnostics.read_text(encoding="utf-8"))
            self.assertEqual(saved["status"], "failed")
            self.assertEqual(saved["returncode"], 7)
            self.assertIn(
                "backend failed",
                Path(saved["stderr_log"]).read_text(encoding="utf-8"),
            )

    def test_success_with_invalid_fasta_preserves_diagnostics(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fq"
            profile = root / "damage_"
            executable = root / "carpedeam"
            diagnostics = root / "diagnostics.json"
            output = root / "out.fa"
            write_fastq(reads)
            write_profile(profile)
            executable.write_bytes(b"fake")

            def fake_run(command, **kwargs):
                if command[-1] == "-h":
                    return subprocess.CompletedProcess(
                        command, 0,
                        stdout=(
                            "usage: carpedeam ancient_assemble\n"
                            "--ancient-damage --unsafe --min-contig-len "
                            "--num-iter-reads-only --num-iterations "
                            "--min-merge-seq-id --min-cov-safe\n"
                        ),
                        stderr="",
                    )
                output.write_text("not FASTA\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0)

            with patch("anvaya.carpedeam_backend.subprocess.run", side_effect=fake_run):
                with self.assertRaisesRegex(ValueError, "invalid CarpeDeam output FASTA"):
                    run_carpedeam_safe(
                        reads, profile, output, root / "tmp", diagnostics,
                        executable=executable,
                    )
            saved = json.loads(diagnostics.read_text(encoding="utf-8"))
            self.assertEqual(saved["status"], "invalid_output")
            self.assertIn("not FASTA", saved["output_error"])


if __name__ == "__main__":
    unittest.main()
