import argparse
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class BaselineRunnerTests(unittest.TestCase):
    def test_manifest_records_parameters_and_artifact_checksums(self) -> None:
        project_root = Path(__file__).resolve().parents[1]
        script = project_root / "experiments" / "02_baseline_accuracy.py"
        spec = importlib.util.spec_from_file_location("baseline_accuracy", script)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.fasta"
            artifact = root / "contigs.fasta"
            manifest_path = root / "manifest.json"
            reference.write_text(">reference\nATGC\n", encoding="utf-8")
            artifact.write_text(">contig\nATGC\n", encoding="utf-8")
            arguments = argparse.Namespace(
                reference=reference,
                k=21,
                threads=4,
                coverage=20.0,
                read_length=150,
                fragment_mean=300,
                fragment_sd=30,
                seed=42,
            )

            module._write_manifest(
                manifest_path,
                arguments=arguments,
                tools={"anvaya": "/usr/bin/anvaya"},
                commands=[("anvaya", ["anvaya", "assemble"])],
                artifacts=[reference, artifact],
            )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

        self.assertEqual(manifest["schema_version"], 1)
        self.assertEqual(manifest["parameters"]["seed"], 42)
        self.assertEqual(manifest["commands"][0]["name"], "anvaya")
        self.assertEqual(len(manifest["artifacts"]), 2)
        self.assertTrue(
            all(len(item["sha256"]) == 64 for item in manifest["artifacts"].values())
        )

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
