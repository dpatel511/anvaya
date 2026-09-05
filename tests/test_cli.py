import io
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path

import anvaya
from anvaya.cli import build_parser, main


class CliTests(unittest.TestCase):
    def test_module_entry_point_displays_overlap_help(self) -> None:
        result = subprocess.run(
            [sys.executable, "-m", "anvaya", "--help"],
            check=False, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("{overlap-assemble}", result.stdout)
        self.assertNotIn("calibrate-events", result.stdout)
        self.assertNotIn("de Bruijn", result.stdout)

    def test_retired_commands_are_rejected(self) -> None:
        for command in ("assemble", "calibrate-events"):
            with self.subTest(command=command), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    main([command])
                self.assertEqual(error.exception.code, 2)

    def test_graph_options_are_rejected_by_overlap_command(self) -> None:
        for option in ("--orientation-aware", "--clean-tips", "--event-report"):
            with self.subTest(option=option), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    build_parser().parse_args(["overlap-assemble", "-i", "reads.fq", "-o", "out.fa", option])

    def test_rejects_missing_right_input(self) -> None:
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main(["overlap-assemble", "-1", "left.fq", "-o", "out.fa"])

    def test_rejects_unequal_pair_counts(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left, right = root / "left.fa", root / "right.fa"
            left.write_text(">one\nACGT\n>two\nACGT\n", encoding="utf-8")
            right.write_text(">one\nACGT\n", encoding="utf-8")
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(["overlap-assemble", "-1", str(left), "-2", str(right), "-o", str(root / "out.fa")])

    def test_package_exports_overlap_and_standalone_damage_math(self) -> None:
        self.assertTrue(callable(anvaya.assemble_overlap_contigs))
        self.assertTrue(callable(anvaya.fit_candidate_damage_model))
        for name in ("assemble", "assemble_file", "build_dbg", "build_bidirected_dbg", "extract_unitigs", "calibrate_event_report"):
            self.assertFalse(hasattr(anvaya, name), name)


if __name__ == "__main__":
    unittest.main()
