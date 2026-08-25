import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class BaselineRunnerTests(unittest.TestCase):
    def test_dry_run_prints_reproducible_tool_commands(self) -> None:
        project_root = Path(__file__).resolve().parents[1]
        script = project_root / "experiments" / "02_baseline_accuracy.py"

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.fasta"
            reference.write_text(">reference\nATGCATGC\n", encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--reference",
                    str(reference),
                    "--output-dir",
                    str(root / "results"),
                    "--dry-run",
                ],
                check=False,
                capture_output=True,
                text=True,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("[art] art_illumina", result.stdout)
        self.assertIn("[anvaya] anvaya assemble -1", result.stdout)
        self.assertIn("[megahit] megahit -1", result.stdout)
        self.assertIn("[metaspades] spades.py --meta", result.stdout)
        self.assertIn("[quast] quast.py", result.stdout)
        self.assertIn("--labels Anvaya,MEGAHIT,metaSPAdes", result.stdout)


if __name__ == "__main__":
    unittest.main()
