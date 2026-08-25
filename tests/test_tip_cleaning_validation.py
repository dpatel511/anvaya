import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class TipCleaningValidationRunnerTests(unittest.TestCase):
    def test_dry_run_prints_every_condition_and_command(self) -> None:
        project_root = Path(__file__).resolve().parents[1]
        script = project_root / "experiments" / "03_tip_cleaning_validation.py"

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.fasta"
            second = root / "second.fasta"
            first.write_text(">first\nATGCATGC\n", encoding="utf-8")
            second.write_text(">second\nGCTAGCTA\n", encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--reference",
                    str(first),
                    "--reference",
                    str(second),
                    "--output-dir",
                    str(root / "results"),
                    "--seeds",
                    "7",
                    "8",
                    "--coverages",
                    "5",
                    "--resume",
                    "--dry-run",
                ],
                check=False,
                capture_output=True,
                text=True,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.count("[condition]"), 4)
        self.assertEqual(result.stdout.count("[art] art_illumina"), 4)
        self.assertEqual(result.stdout.count("[baseline] anvaya assemble"), 4)
        self.assertEqual(result.stdout.count("[cleaned] anvaya assemble"), 4)
        self.assertEqual(result.stdout.count("--clean-tips"), 4)
        self.assertEqual(result.stdout.count("[quast] quast.py"), 4)


if __name__ == "__main__":
    unittest.main()
