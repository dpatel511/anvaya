import io
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.cli import main


class PairedCliTests(unittest.TestCase):
    def test_assemble_accepts_paired_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left = root / "left.fastq"
            right = root / "right.fastq"
            output = root / "unitigs.fasta"
            left.write_text("@read-1/1\nATGCA\n+\nIIIII\n", encoding="utf-8")
            right.write_text("@read-1/2\nTGCAT\n+\nIIIII\n", encoding="utf-8")
            stdout = io.StringIO()

            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                exit_code = main(
                    [
                        "assemble",
                        "-1",
                        str(left),
                        "-2",
                        str(right),
                        "--k",
                        "3",
                        "-o",
                        str(output),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertTrue(output.exists())
            self.assertIn("input_files=2", stdout.getvalue())
            self.assertIn("reads=2", stdout.getvalue())

    def test_rejects_missing_right_input(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            left = Path(directory) / "left.fastq"
            left.write_text("@read-1/1\nATGCA\n+\nIIIII\n", encoding="utf-8")

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(["assemble", "-1", str(left), "--k", "3", "-o", "out.fa"])

    def test_rejects_unequal_pair_counts(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left = root / "left.fastq"
            right = root / "right.fastq"
            left.write_text(
                "@read-1/1\nATGCA\n+\nIIIII\n@read-2/1\nATGCA\n+\nIIIII\n",
                encoding="utf-8",
            )
            right.write_text("@read-1/2\nTGCAT\n+\nIIIII\n", encoding="utf-8")

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(
                    [
                        "assemble",
                        "-1",
                        str(left),
                        "-2",
                        str(right),
                        "--k",
                        "3",
                        "-o",
                        str(root / "out.fa"),
                    ]
                )


if __name__ == "__main__":
    unittest.main()
