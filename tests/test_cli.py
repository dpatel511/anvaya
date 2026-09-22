import csv
import io
import json
import random
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

import anvaya
from anvaya.cli import build_parser, main


class CliTests(unittest.TestCase):
    def test_low_quality_ry_rescue_is_opt_in_and_requires_progressive_audit(self):
        base = ["overlap-assemble", "-i", "reads.fq", "-o", "out.fa"]
        self.assertFalse(build_parser().parse_args(base).progressive_low_quality_ry_rescue)
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            main(base + ["--progressive-low-quality-ry-rescue"])
        self.assertEqual(error.exception.code, 2)

    def test_support_two_quality_switches(self) -> None:
        base = ["overlap-assemble", "-i", "reads.fq", "-o", "out.fa"]
        for seed in (False, True):
            for recruit in (False, True):
                args = build_parser().parse_args(base + [
                    "--" + ("" if seed else "no-") + "support-two-seed-quality-filter",
                    "--" + ("" if recruit else "no-") + "support-two-recruit-quality-filter",
                ])
                self.assertEqual(args.support_two_seed_quality_filter, seed)
                self.assertEqual(args.support_two_recruit_quality_filter, recruit)
        defaults = build_parser().parse_args(base)
        self.assertFalse(defaults.support_two_seed_quality_filter)
        self.assertTrue(defaults.support_two_recruit_quality_filter)

    def test_module_entry_point_displays_overlap_help(self) -> None:
        result = subprocess.run(
            [sys.executable, "-m", "anvaya", "--help"],
            check=False, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("damage-assemble", result.stdout)
        self.assertIn("carpedeam-assemble", result.stdout)
        self.assertIn("overlap-assemble", result.stdout)
        self.assertNotIn("calibrate-events", result.stdout)
        self.assertNotIn("de Bruijn", result.stdout)

    def test_carpedeam_command_uses_audited_safe_backend(self) -> None:
        result = {
            "elapsed_seconds": 1.25,
            "input_scan": {"reads": 12},
            "output_contigs": 3,
        }
        argv = [
            "carpedeam-assemble",
            "--input", "reads.fq.gz",
            "--profile-prefix", "damage_",
            "--output", "contigs.fa",
            "--temporary-directory", "backend-tmp",
            "--diagnostics", "backend.json",
            "--executable", "/opt/carpedeam",
            "--threads", "4",
        ]
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch("anvaya.cli.run_carpedeam_safe", return_value=result) as run:
            with redirect_stdout(stdout), redirect_stderr(stderr):
                self.assertEqual(main(argv), 0)
        run.assert_called_once_with(
            Path("reads.fq.gz"), Path("damage_"), Path("contigs.fa"),
            Path("backend-tmp"), Path("backend.json"),
            executable="/opt/carpedeam", threads=4, minimum_contig_length=31,
        )
        self.assertIn("reads=12", stdout.getvalue())
        self.assertIn("contigs=3", stdout.getvalue())
        self.assertIn("safe-mode comparator", stderr.getvalue())

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

    def test_damage_assemble_writes_one_layout_and_audits(self) -> None:
        rng = random.Random(95020)
        truth = list("".join(rng.choice("ACGT") for _ in range(170)))
        truth[50] = "C"
        truth = "".join(truth)
        reads = (truth[:90], truth[:90], "T" + truth[51:170])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastq = root / "reads.fq"
            fastq.write_text("".join(
                f"@r{i}\n{sequence}\n+\n{'D' * len(sequence)}\n"
                for i, sequence in enumerate(reads)
            ), encoding="utf-8")
            prefix = root / "profile"
            header = [f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b]
            for end, transition in (("5p", "C>T"), ("3p", "G>A")):
                row = ["0.4" if column == transition else "0" for column in header]
                Path(f"{prefix}{end}.prof").write_text(
                    "\t".join(header) + "\n" + "\t".join(row) + "\n",
                    encoding="utf-8",
                )
            before = root / "before.fa"
            after = root / "after.fa"
            diagnostics = root / "diagnostics.json"
            placements = root / "placements.tsv"
            decisions = root / "decisions.tsv"
            assembled = root / "assembled.fa"
            unresolved = root / "unresolved.fa"
            with redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
                self.assertEqual(main([
                    "damage-assemble", "-i", str(fastq),
                    "--profile-prefix", str(prefix),
                    "--output-before", str(before),
                    "--output-consensus", str(after),
                    "--diagnostics", str(diagnostics),
                    "--placements-report", str(placements),
                    "--consensus-report", str(decisions),
                    "--output-assembled", str(assembled),
                    "--output-unresolved", str(unresolved),
                ]), 0)
            self.assertEqual(before.read_bytes(), after.read_bytes())
            summary = json.loads(diagnostics.read_text(encoding="utf-8"))
            self.assertEqual(summary["input_reads"], 3)
            self.assertEqual(summary["emitted_contigs"], 1)
            self.assertEqual(summary["assembled_contigs"], 1)
            self.assertEqual(summary["unresolved_contigs"], 0)
            self.assertEqual(assembled.read_bytes(), after.read_bytes())
            self.assertEqual(unresolved.read_text(encoding="utf-8"), "")
            self.assertGreater(summary["graph"]["candidate_classifications"]["damage_compatible"], 0)
            self.assertEqual(
                set(summary["stage_seconds"]),
                {
                    "indexing", "candidate_edges", "containment", "phase_audit",
                    "projection", "provenance", "consensus",
                },
            )
            self.assertEqual(len(placements.read_text(encoding="utf-8").splitlines()), 4)

    def test_damage_assemble_filters_layout_before_numbering_audits(self) -> None:
        rng = random.Random(95100)
        reads = (
            "".join(rng.choice("ACGT") for _ in range(31)),
            "".join(rng.choice("ACGT") for _ in range(80)),
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastq = root / "reads.fq"
            fastq.write_text("".join(
                f"@r{i}\n{sequence}\n+\n{'D' * len(sequence)}\n"
                for i, sequence in enumerate(reads)
            ), encoding="utf-8")
            prefix = root / "profile"
            header = [f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b]
            for end in ("5p", "3p"):
                Path(f"{prefix}{end}.prof").write_text(
                    "\t".join(header) + "\n" + "\t".join("0" for _ in header) + "\n",
                    encoding="utf-8",
                )
            before = root / "before.fa"
            after = root / "after.fa"
            diagnostics = root / "diagnostics.json"
            placements = root / "placements.tsv"
            with redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
                self.assertEqual(main([
                    "damage-assemble", "-i", str(fastq),
                    "--profile-prefix", str(prefix),
                    "--output-before", str(before),
                    "--output-consensus", str(after),
                    "--diagnostics", str(diagnostics),
                    "--placements-report", str(placements),
                    "--min-output-length", "50",
                ]), 0)

            expected = f">unitig_1 length=80\n{reads[1]}\n"
            self.assertEqual(before.read_text(encoding="utf-8"), expected)
            self.assertEqual(after.read_text(encoding="utf-8"), expected)
            with placements.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["contig_id"] for row in rows}, {"unitig_1"})
            self.assertEqual({row["read_name"] for row in rows}, {"r1"})
            summary = json.loads(diagnostics.read_text(encoding="utf-8"))
            self.assertEqual(summary["layout_contigs"], 2)
            self.assertEqual(summary["filtered_short_contigs"], 1)
            self.assertEqual(summary["emitted_contigs"], 1)
            self.assertEqual(summary["consensus"]["contigs"], 1)

    def test_damage_assemble_rejects_fasta_and_read_limit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "reads.fa"
            fasta.write_text(">r\n" + "A" * 40 + "\n", encoding="utf-8")
            prefix = root / "profile"
            header = [f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b]
            for end in ("5p", "3p"):
                Path(f"{prefix}{end}.prof").write_text(
                    "\t".join(header) + "\n" + "\t".join("0" for _ in header) + "\n",
                    encoding="utf-8",
                )
            args = ["damage-assemble", "-i", str(fasta), "--profile-prefix", str(prefix),
                    "--output-before", str(root / "before.fa"),
                    "--output-consensus", str(root / "after.fa"),
                    "--diagnostics", str(root / "diagnostics.json"), "--max-reads", "1"]
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(args)
            fastq = root / "reads.fq"
            fastq.write_text(
                "@a\n" + "A" * 40 + "\n+\n" + "D" * 40 + "\n"
                "@b\n" + "C" * 40 + "\n+\n" + "D" * 40 + "\n",
                encoding="utf-8",
            )
            args[2] = str(fastq)
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(args)


if __name__ == "__main__":
    unittest.main()
