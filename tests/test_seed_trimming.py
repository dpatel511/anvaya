"""Synthetic tests for trimming in original seed coordinates (0-based, half-open)."""

import io
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.overlap_adaptive_rescue import project_high_confidence_support_two_rescue
from anvaya.cli import main
from anvaya.overlap_assembly import _Alignment
from anvaya.overlap_progressive import (
    ProgressiveRawCluster, ProgressiveSequencePool, extend_progressive_raw_clusters,
    _seed_trim_bounds,
)
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


class SeedTrimmingTests(unittest.TestCase):
    sequence = "ACGTTGCACTGATCGGACCTAGTA"

    def test_fixed_membership_matches_baseline_and_retains_order(self):
        pool = ProgressiveSequencePool.from_reads([
            Read("seed", self.sequence, (0,) * 6 + (30,) * 12 + (0,) * 6),
            Read("partner", self.sequence[6:18], (30,) * 12),
        ])
        options = dict(anchor_k=3, anchors_per_read=16, minimum_overlap=8)
        baseline, _ = project_high_confidence_support_two_rescue(pool, [], **options)
        for minimum in (12, 13):
            before = []
            after, diagnostics = project_high_confidence_support_two_rescue(
                pool, [], **options, trim_seed_tails=True,
                fixed_membership_before=before, minimum_output_length=minimum,
            )
            self.assertEqual([r.name for r in before], [r.name for r in after])
            self.assertEqual(before, baseline if minimum == 12 else [])
            self.assertEqual([r.sequence for r in after], [self.sequence[6:18]] if minimum == 12 else [])
            self.assertEqual(diagnostics.fixed_membership_short_exclusions, int(minimum == 13))
            self.assertEqual(diagnostics.rescue_contigs,
                diagnostics.contained_by_primary + diagnostics.redundant_with_rescue
                + diagnostics.primary_extension_exclusions + diagnostics.novel_contigs
                + diagnostics.fixed_membership_short_exclusions)

    def test_fixed_membership_does_not_reclassify_primary_extension(self):
        # The untrimmed seed extends primary; trimming removes part of that match.
        pool = ProgressiveSequencePool.from_reads([
            Read("seed", self.sequence, (0,) * 6 + (30,) * 18),
            Read("partner", self.sequence[6:18], (30,) * 12),
        ])
        primary = [Read("primary", self.sequence[:20])]
        before = []
        after, diagnostics = project_high_confidence_support_two_rescue(
            pool, primary, anchor_k=3, anchors_per_read=16, minimum_overlap=8,
            trim_seed_tails=True, fixed_membership_before=before,
        )
        self.assertEqual(before, primary)
        self.assertEqual(after, primary)
        self.assertEqual(diagnostics.primary_extension_exclusions, 1)
        self.assertEqual(diagnostics.novel_contigs, 0)

    def test_fixed_membership_empty_pool_keeps_primary(self):
        primary = [Read("primary", self.sequence)]
        before = []
        after, _ = project_high_confidence_support_two_rescue(
            ProgressiveSequencePool.from_reads([]), primary,
            trim_seed_tails=True, fixed_membership_before=before,
        )
        self.assertEqual(before, primary)
        self.assertEqual(after, primary)

    @unittest.skipUnless(shutil.which("bash") and Path("/usr/bin/time").exists(), "WSL/Linux comparison script")
    def test_fixed_comparison_script_on_synthetic_fixture(self):
        rng = random.Random(804)
        sequence = "".join(rng.choice("ACGT") for _ in range(140))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastq = root / "synthetic.fq"
            fastq.write_text(
                f"@seed\n{sequence[:100]}\n+\n!{'?' * 99}\n"
                f"@partner\n{sequence[40:]}\n+\n{'?' * 100}\n", encoding="ascii",
            )
            import os
            environment = dict(os.environ, ANVAYA_PYTHON=sys.executable)
            result = subprocess.run(
                ["bash", "experiments/support_two_quality_ablation.sh", str(fastq), str(root / "out"), "fixed"],
                cwd=Path(__file__).resolve().parents[1], env=environment,
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            before = (root / "out/fixed/before-contigs.fasta").read_text().splitlines()
            after = (root / "out/fixed/support-two-contigs.fasta").read_text().splitlines()
            self.assertEqual([line.split()[0] for line in before if line.startswith(">")],
                             [line.split()[0] for line in after if line.startswith(">")])
            self.assertEqual("".join(line for line in before if not line.startswith(">")), sequence)
            self.assertEqual("".join(line for line in after if not line.startswith(">")), sequence[1:])

    def test_bounds_match_per_base_support_oracle(self):
        rng = random.Random(4908)
        for _ in range(200):
            length = rng.randrange(10, 50)
            qualities = tuple(rng.choice((0, 19, 20, 30)) for _ in range(length))
            seed = Read("seed", "A" * length, qualities)
            offset = rng.choice((0, 0, 5))
            contig_length = offset + length + rng.choice((0, 0, 5))
            members = [
                _Alignment(index, "A" * rng.randrange(2, length), rng.randrange(-5, contig_length), 2)
                for index in range(1, 5)
            ]
            molecules = [0, 0, 2, 3, 4]  # same-molecule evidence must not protect a tail
            coverage = [False] * length
            for member in members:
                if molecules[member.read_index] != 0:
                    for position in range(length):
                        if member.offset <= position + offset < member.offset + len(member.sequence):
                            coverage[position] = True
            left, right = 0, contig_length
            if any(coverage):
                while left < right and 0 <= left - offset < length:
                    position = left - offset
                    if coverage[position] or qualities[position] >= 20:
                        break
                    left += 1
                while right > left and 0 <= right - offset - 1 < length:
                    position = right - offset - 1
                    if coverage[position] or qualities[position] >= 20:
                        break
                    right -= 1
            self.assertEqual(_seed_trim_bounds(seed, contig_length, offset, members, molecules, 0, 20), (left, right))

    def test_cli_projects_trimmed_sequence_and_diagnostics(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastq = root / "synthetic.fq"
            fastq.write_text(
                f"@seed\n{self.sequence}\n+\n{'!' * 6}{'?' * 12}{'!' * 6}\n"
                f"@partner\n{self.sequence[6:18]}\n+\n{'?' * 12}\n", encoding="ascii",
            )
            for trim in (False, True):
                output = root / f"rescue-{trim}.fa"
                args = [
                    "overlap-assemble", "-i", str(fastq), "-o", str(root / "primary.fa"),
                    "--anchor-k", "3", "--anchors-per-read", "16", "--min-overlap", "8",
                    "--min-cluster-size", "5", "--min-output-length", "1", "--max-contig-iterations", "0",
                    "--progressive-raw-phase-audit", "--selective-support-rescue-audit",
                    "--ranked-extension", "--reciprocal-best-extension", "--damage-aware-ranking",
                    "--high-confidence-support-two-rescue-audit",
                    "--high-confidence-support-two-rescue-projection", str(output),
                ]
                if trim:
                    args.append("--support-two-trim-seed-tails")
                    args.extend(["--support-two-fixed-membership-before", str(root / "before.fa")])
                log = io.StringIO()
                with redirect_stdout(log), redirect_stderr(io.StringIO()):
                    self.assertEqual(main(args), 0)
                self.assertEqual(output.read_text().splitlines()[1], self.sequence[6:18] if trim else self.sequence)
                self.assertIn(f"overlap_support_two_rescue_seed_trimmed_contigs={int(trim)}", log.getvalue())
                if trim:
                    self.assertEqual((root / "before.fa").read_text().splitlines()[1], self.sequence)

    def test_cli_rejects_invalid_trim_options_before_opening_input(self):
        base = ["overlap-assemble", "-i", "does-not-exist.fq", "-o", "out.fa", "--support-two-trim-seed-tails"]
        for extra, message in (([], "seed trimming requires"), (["--support-two-seed-quality-filter"], "seed rejection and trimming")):
            error = io.StringIO()
            with redirect_stderr(error), self.assertRaises(SystemExit) as result:
                main(base + extra)
            self.assertEqual(result.exception.code, 2)
            self.assertIn(message, error.getvalue())

    def extend(self, seed, partner, offset, **kwargs):
        pool = ProgressiveSequencePool.from_reads([seed, partner])
        cluster = ProgressiveRawCluster(0, (0, 1), (_Alignment(1, partner.sequence, offset, 2),))
        updated, diagnostics = extend_progressive_raw_clusters(
            pool, (cluster,), anchor_k=3, anchors_per_read=16, minimum_overlap=8,
            minimum_consensus_support=2, minimum_correction_support=2,
            reciprocal_best_extension=False, trim_seed_min_base_quality=20, **kwargs,
        )
        self.assertIs(updated.records[0].raw, seed)
        self.assertEqual(pool.records[0].current, seed)
        return updated.records[0], diagnostics

    def test_both_tails_trim_to_supported_interval_and_preserve_raw_coordinates(self):
        for qualities in (None, (0,) * 6 + (30,) * 12 + (19,) * 6):
            with self.subTest(qualities=qualities):
                seed = Read("seed", self.sequence, qualities)
                record, diagnostics = self.extend(seed, Read("partner", self.sequence[6:18], (30,) * 12), 6)
                self.assertEqual(record.current.sequence, self.sequence[6:18])
                self.assertEqual(record.seed_interval, (6, 18))
                self.assertEqual(record.seed_offset, -6)
                self.assertEqual(diagnostics.seed_trimmed_contigs, 1)
                self.assertEqual(diagnostics.seed_trimmed_left_bases, 6)
                self.assertEqual(diagnostics.seed_trimmed_right_bases, 6)

    def test_trim_stops_at_high_quality_base_and_leaves_internal_low_quality(self):
        seed = Read("seed", self.sequence, (0, 20, 0, 0, 0, 0) + (30,) * 18)
        record, diagnostics = self.extend(seed, Read("partner", self.sequence[6:18], (30,) * 12), 6)
        self.assertEqual(record.current.sequence, self.sequence[1:])
        self.assertEqual(record.seed_interval, (1, 24))
        self.assertEqual(diagnostics.seed_trimmed_left_bases, 1)

    def test_supported_low_quality_seed_bases_are_retained(self):
        seed = Read("seed", self.sequence, (30,) * 6 + (0,) * 12 + (30,) * 6)
        record, diagnostics = self.extend(seed, Read("partner", self.sequence[6:18], (30,) * 12), 6)
        self.assertEqual(record.current.sequence, self.sequence)
        self.assertEqual(record.seed_interval, (0, 24))
        self.assertEqual(diagnostics.seed_trimmed_contigs, 0)

    def test_extension_and_reverse_orientation_preserve_seed_mapping(self):
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                seed_sequence = self.sequence[:18]
                partner_sequence = self.sequence[6:]
                qualities = (0,) * 6 + (30,) * 12
                offset = 6
                if reverse:
                    seed_sequence = reverse_complement(seed_sequence)
                    partner_sequence = reverse_complement(partner_sequence)
                    qualities = qualities[::-1]
                    offset = -6
                record, diagnostics = self.extend(
                    Read("seed", seed_sequence, qualities),
                    Read("partner", partner_sequence, (30,) * 18), offset,
                )
                expected = reverse_complement(self.sequence[6:]) if reverse else self.sequence[6:]
                self.assertEqual(record.current.sequence, expected)
                self.assertEqual(record.seed_interval, (0, 12) if reverse else (6, 18))
                self.assertEqual(record.seed_offset, 6 if reverse else -6)
                self.assertEqual(diagnostics.seed_trimmed_left_bases + diagnostics.seed_trimmed_right_bases, 6)

    def test_projection_trimming_and_minimum_output_length(self):
        pool = ProgressiveSequencePool.from_reads([
            Read("seed", self.sequence, (0,) * 6 + (30,) * 12 + (0,) * 6),
            Read("partner", self.sequence[6:18], (30,) * 12),
        ])
        for minimum in (12, 13):
            with self.subTest(minimum=minimum):
                projected, diagnostics = project_high_confidence_support_two_rescue(
                    pool, [], anchor_k=3, anchors_per_read=16, minimum_overlap=8,
                    trim_seed_tails=True, minimum_output_length=minimum,
                )
                self.assertEqual([read.sequence for read in projected], [self.sequence[6:18]] if minimum == 12 else [])
                self.assertEqual(diagnostics.seed_trimmed_contigs, 1)
                self.assertEqual(diagnostics.trimmed_short_contigs, int(minimum == 13))

    def test_rejection_and_trimming_cannot_be_combined(self):
        with self.assertRaisesRegex(ValueError, "seed rejection and trimming"):
            project_high_confidence_support_two_rescue(
                ProgressiveSequencePool.from_reads([]), [],
                filter_seed_quality=True, trim_seed_tails=True,
            )
