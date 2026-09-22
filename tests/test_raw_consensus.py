"""Raw-fragment consensus tests; all coordinates are 0-based and half-open."""

import math
import io
import random
import tempfile
import unittest
from dataclasses import replace
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.raw_consensus import (
    DamageProfile, RawPlacement, observation_likelihoods, project_raw_consensus,
)
from anvaya.overlap_progressive import ProgressiveSequencePool
from anvaya.overlap_progressive import (
    ProgressiveRawCluster, _raw_placement, extend_progressive_raw_clusters,
    iterate_progressive_raw_extension,
)
from anvaya.overlap_assembly import _Alignment
from anvaya.cli import main
from anvaya.sequences import reverse_complement
from anvaya.reads import Read


class RawConsensusTests(unittest.TestCase):
    def test_palindromic_read_keeps_explicit_reverse_orientation(self):
        read = Read("palindrome", "ACGT", (10, 20, 30, 40))
        placement = _raw_placement(_Alignment(0, "ACGT", 2, 2, True), read)
        self.assertTrue(placement.reverse)
        with self.assertRaisesRegex(ValueError, "ambiguous"):
            _raw_placement(_Alignment(0, "ACGT", 2, 2), read)

    def test_reverse_placement_survives_extension_and_trimming(self):
        sequence = "ACGTTGCACTGATCGGACCTAGTA"
        reads = [Read("seed", sequence[:18], (0,) * 6 + (30,) * 12),
                 Read("partner", reverse_complement(sequence[6:]), (30,) * 18)]
        pool = ProgressiveSequencePool.from_reads(reads)
        cluster = ProgressiveRawCluster(0, (0, 1), (_Alignment(1, sequence[6:], 6, 2, True),))
        updated, _ = extend_progressive_raw_clusters(
            pool, (cluster,), anchor_k=3, anchors_per_read=16, minimum_overlap=8,
            minimum_consensus_support=2, minimum_correction_support=2,
            reciprocal_best_extension=False, trim_seed_min_base_quality=20, track_raw_placements=True,
        )
        record = updated.records[0]
        self.assertEqual(record.current.sequence, sequence[6:])
        self.assertEqual(record.raw_placements, (RawPlacement(0, -6, False, 6, 18), RawPlacement(1, 0, True, 0, 18)))
        self.assertIs(record.raw, reads[0])

    def test_later_left_extension_shifts_existing_placements_and_global_ids(self):
        rng = random.Random(319)
        sequence = "".join(rng.choice("ACGT") for _ in range(80))
        reads = [Read("seed", sequence[20:], (30,) * 60)] + [Read(str(i), sequence[:60], (30,) * 60) for i in range(3)]
        pool = ProgressiveSequencePool.from_reads(reads)
        center = replace(pool.records[0].corrected(reads[0]), raw_placements=(RawPlacement(0, 0, False, 0, 60),))
        pool = pool.replace_record(center)
        updated, diagnostics = iterate_progressive_raw_extension(
            pool, anchor_k=5, anchors_per_read=30, minimum_overlap=20,
            minimum_consensus_support=3, maximum_iterations=1,
        )
        self.assertEqual(diagnostics.extended_centers, 1)
        record = updated.records[0]
        self.assertEqual(record.current.sequence, sequence)
        self.assertEqual(record.raw_placements[0].offset, 20)
        self.assertEqual({p.read_index for p in record.raw_placements}, {0, 1, 2, 3})
        for placement in record.raw_placements:
            raw = updated.records[placement.read_index].raw.sequence
            for position in range(placement.read_start, placement.read_stop):
                contig_position = placement.offset + (len(raw) - 1 - position if placement.reverse else position)
                observed = reverse_complement(raw[position]) if placement.reverse else raw[position]
                self.assertEqual(record.current.sequence[contig_position], observed)

    def test_high_depth_and_ambiguous_mapping(self):
        pool = self.pool("C" * 1000)
        output, diagnostics = project_raw_consensus(pool, DamageProfile((), ()))
        self.assertEqual(output[0].sequence, "C")
        self.assertEqual(diagnostics.changed_bases, 1)
        pool = self.pool("TCC", positions=[0, 5, 5])
        center = pool.records[0]
        pool = pool.replace_record(replace(center, raw_placements=center.raw_placements + (RawPlacement(1, -4, False, 5, 6),)))
        output, diagnostics = project_raw_consensus(pool, DamageProfile((0.4,), ()))
        self.assertEqual(output[0].sequence, "T")
        self.assertEqual(diagnostics.ambiguous_placement_molecules, 1)

    def test_unobserved_allele_is_not_invented(self):
        output, diagnostics = project_raw_consensus(self.pool("TTTTT"), DamageProfile((1.0,), ()))
        self.assertEqual(output[0].sequence, "T")
        self.assertEqual(diagnostics.changed_bases, 0)

    def test_uncertain_quality_weighted_majority_does_not_force_a_change(self):
        pool = self.pool("TCC", qualities=[30, 20, 20], positions=[5, 5, 5])
        output, diagnostics = project_raw_consensus(pool, DamageProfile((), ()))
        self.assertEqual(output[0].sequence, "T")
        self.assertEqual(diagnostics.posterior_rejections, 1)

    def test_conflicting_duplicate_molecule_is_excluded(self):
        pool = self.pool("TCCCC", molecules=[0, 0, 1, 2, 3], positions=[5] * 5)
        output, diagnostics = project_raw_consensus(pool, DamageProfile((), ()))
        self.assertEqual(output[0].sequence, "C")
        self.assertEqual(diagnostics.conflicting_molecule_observations, 1)

    def test_cli_fixed_layout_damage_and_zero_damage_control(self):
        rng = random.Random(388)
        reference = "".join(rng.choice("ACGT") for _ in range(140))
        reference = reference[:40] + "C" + reference[41:]
        damaged = "T" + reference[41:]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prefix = root / "profile"
            header = [f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b]
            for end, transition in (("5p", "C>T"), ("3p", "G>A")):
                Path(f"{prefix}{end}.prof").write_text("\t".join(header) + "\n" + "\t".join("0.4" if c == transition else "0" for c in header) + "\n")
            fastq = root / "synthetic.fq"
            fragments = [damaged] * 5 + [reference[:100]] * 2
            fastq.write_text("".join(f"@{i}\n{seq}\n+\n{'?' * len(seq)}\n" for i, seq in enumerate(fragments)))
            args = ["overlap-assemble", "-i", str(fastq), "-o", str(root / "primary.fa"),
                "--anchor-k", "5", "--anchors-per-read", "40", "--min-overlap", "20",
                "--min-cluster-size", "5", "--max-rounds", "1", "--max-contig-iterations", "0",
                "--ranked-extension", "--reciprocal-best-extension", "--damage-aware-ranking",
                "--progressive-raw-phase-audit", "--progressive-raw-phase-projection", str(root / "before.fa"),
                "--raw-consensus-profile-prefix", str(prefix), "--raw-consensus-projection", str(root / "after.fa"),
                "--raw-consensus-quality-control", str(root / "quality.fa"),
                "--raw-consensus-report", str(root / "decisions.tsv"),
                "--raw-consensus-placements", str(root / "placements.tsv"),
                "--raw-consensus-mixture-report", str(root / "mixture.tsv"),
                "--raw-consensus-linkage-report", str(root / "linkage.tsv")]
            with redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
                self.assertEqual(main(args), 0)
            sequences = {name: "".join(line for line in (root / name).read_text().splitlines() if not line.startswith(">"))
                         for name in ("before.fa", "after.fa", "quality.fa")}
            self.assertEqual(sequences["before.fa"], reference[:40] + "T" + reference[41:])
            self.assertEqual(sequences["quality.fa"], sequences["before.fa"])
            self.assertEqual(sequences["after.fa"], reference)
            self.assertIn("unitig_1\t40\tT\tC", (root / "decisions.tsv").read_text())
            self.assertEqual(len((root / "placements.tsv").read_text().splitlines()), 8)
            self.assertIn("log_bf_mixture_vs_single", (root / "mixture.tsv").read_text())
            self.assertIn("neighbour_position_0based", (root / "linkage.tsv").read_text())
            unguarded_outputs = {name: (root / name).read_bytes() for name in ("before.fa", "after.fa", "quality.fa")}
            log_output = io.StringIO()
            with redirect_stderr(io.StringIO()), redirect_stdout(log_output):
                self.assertEqual(main(args + ["--raw-consensus-linked-allele-guard"]), 0)
            self.assertIn("raw_consensus_linked_allele_guard=true", log_output.getvalue())
            self.assertIn("raw_consensus_linked_allele_rejections=0", log_output.getvalue())
            self.assertEqual(unguarded_outputs, {name: (root / name).read_bytes() for name in unguarded_outputs})

    def test_likelihood_is_normalized_and_reduces_to_quality_error(self):
        for profile in (DamageProfile((), ()), DamageProfile((0.4,), (0.2,))):
            for reverse in (False, True):
                for position in (0, 4, 9):
                    vectors = [observation_likelihoods(base, 30, position, 10, reverse, profile) for base in "ACGT"]
                    for true_index in range(4):
                        self.assertAlmostEqual(sum(vector[true_index] for vector in vectors), 1)
        values = observation_likelihoods("C", 30, 0, 10, False, DamageProfile((), ()))
        self.assertAlmostEqual(values[1], 0.999)
        self.assertAlmostEqual(values[0], 0.001 / 3)

    def test_damage_direction_and_reverse_complement(self):
        profile = DamageProfile((0.4,), (0.2,))
        values = observation_likelihoods("T", 30, 0, 10, False, profile)
        self.assertAlmostEqual(values[1], 0.4 * 0.999 + 0.6 * 0.001 / 3)
        reverse = observation_likelihoods("T", 30, 0, 10, True, profile)
        self.assertEqual(reverse, tuple(reversed(values)))
        interior = observation_likelihoods("T", 30, 4, 10, False, profile)
        self.assertAlmostEqual(interior[1], 0.001 / 3)

    def pool(self, bases, *, qualities=None, molecules=None, positions=None):
        qualities = qualities or [30] * len(bases)
        positions = positions or [0] * len(bases)
        reads = [Read(str(i), "A" * p + base + "A" * 10,
                      None if q is None else (q,) * (p + 11))
                 for i, (base, p, q) in enumerate(zip(bases, positions, qualities))]
        pool = ProgressiveSequencePool.from_reads(reads, molecules)
        center = pool.records[0].corrected(Read("contig", "T"))
        center = replace(center, raw_placements=tuple(
            RawPlacement(i, -p, False, p, p + 1) for i, p in enumerate(positions)
        ))
        return pool.replace_record(center)

    def test_terminal_damage_is_corrected_using_internal_raw_support(self):
        pool = self.pool("TTTCC", positions=[0, 0, 0, 5, 5])
        damaged, diagnostics = project_raw_consensus(pool, DamageProfile((0.4,), ()))
        quality_only, _ = project_raw_consensus(pool, DamageProfile((), ()))
        self.assertEqual(damaged[0].sequence, "C")
        self.assertEqual(quality_only[0].sequence, "T")
        self.assertEqual(diagnostics.changed_bases, 1)
        self.assertEqual(pool.records[0].current.sequence, "T")

    def test_two_internal_alleles_are_flagged_not_collapsed(self):
        pool = self.pool("TTCCCCCC", positions=[5] * 8)
        output, diagnostics = project_raw_consensus(pool, DamageProfile((0.4,), ()))
        self.assertEqual(output[0].sequence, "T")
        self.assertEqual(diagnostics.allele_conflict_sites, 1)

    def test_internal_mixture_guard_across_abundances_and_profile_strengths(self):
        for minor_count in (2, 3, 5):
            for major_count in (5, 10, 20):
                for rate in (0, 0.1, 0.4, 0.9):
                    with self.subTest(minor=minor_count, major=major_count, rate=rate):
                        pool = self.pool("T" * minor_count + "C" * major_count,
                                         positions=[5] * (minor_count + major_count))
                        output, diagnostics = project_raw_consensus(pool, DamageProfile((rate,), (rate,)))
                        self.assertEqual(output[0].sequence, "T")
                        self.assertEqual(diagnostics.allele_conflict_sites, 1)

    def test_wrong_profile_can_hide_terminal_strain_evidence(self):
        # A true T-majority/C-minority mixture can have exactly the observations
        # of damaged C molecules. This is a known limitation, not strain safety.
        pool = self.pool("TTTCC", positions=[0, 0, 0, 5, 5])
        zero, zero_diagnostics = project_raw_consensus(pool, DamageProfile((), ()))
        wrong, wrong_diagnostics = project_raw_consensus(pool, DamageProfile((0.4,), ()))
        self.assertEqual(zero[0].sequence, "T")
        self.assertEqual(zero_diagnostics.allele_conflict_sites, 1)
        self.assertEqual(wrong[0].sequence, "C")
        self.assertEqual(wrong_diagnostics.allele_conflict_sites, 0)

    def test_extreme_profile_cannot_replace_required_observed_support(self):
        for rate in (0, 0.4, 0.9, 1):
            with self.subTest(rate=rate):
                pool = self.pool("TTTC", positions=[0, 0, 0, 5])
                output, _ = project_raw_consensus(pool, DamageProfile((rate,), ()))
                self.assertEqual(output[0].sequence, "T")

    def test_duplicate_molecules_cannot_manufacture_support(self):
        pool = self.pool("TCCCC", molecules=[0, 1, 1, 1, 1], positions=[0, 5, 5, 5, 5])
        output, diagnostics = project_raw_consensus(pool, DamageProfile((0.4,), ()))
        self.assertEqual(output[0].sequence, "T")
        self.assertEqual(diagnostics.insufficient_support_sites, 1)

    def test_missing_and_uninformative_quality_are_not_invented(self):
        pool = self.pool("TCCCC", qualities=[30, None, None, 0, 1])
        output, diagnostics = project_raw_consensus(pool, DamageProfile((), ()))
        self.assertEqual(output[0].sequence, "T")
        self.assertEqual(diagnostics.skipped_quality_observations, 4)

    def test_untracked_layout_is_rejected(self):
        pool = ProgressiveSequencePool.from_reads([Read("seed", "ACGT")])
        pool = pool.replace_record(pool.records[0].corrected(Read("contig", "ACGT")))
        with self.assertRaisesRegex(ValueError, "raw placements"):
            project_raw_consensus(pool, DamageProfile((), ()))

    def test_profile_validation(self):
        for values in ((math.nan,), (-0.1,), (1.1,)):
            with self.assertRaises(ValueError):
                DamageProfile(values, ())
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "profile"
            header = "\t".join(f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b)
            for end, transition in (("5p", "C>T"), ("3p", "G>A")):
                row = ["0.4" if column == transition else "0" for column in header.split("\t")]
                Path(f"{prefix}{end}.prof").write_text(header + "\n" + "\t".join(row) + "\n")
            self.assertEqual(DamageProfile.from_prefix(prefix), DamageProfile((0.4,), (0.4,)))
            Path(f"{prefix}5p.prof").write_text(header + "\n0\n")
            with self.assertRaises(ValueError):
                DamageProfile.from_prefix(prefix)
