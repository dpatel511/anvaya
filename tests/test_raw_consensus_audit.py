"""Synthetic alignment fixtures only; no actual FASTQ datasets."""

import csv
import importlib.util
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None

SCRIPT = Path(__file__).resolve().parents[1] / "experiments/raw_consensus_audit.py"
if pysam is not None:
    spec = importlib.util.spec_from_file_location("raw_consensus_audit", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)


@unittest.skipIf(pysam is None, "audit tests require optional pysam")
class RawConsensusAuditTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.paths = []
        for name, sequence in (("before", "TCGA"), ("quality", "CCGA"),
                               ("damage", "CCTA"), ("reference", "CCGT")):
            path = self.root / f"{name}.fa"
            path.write_text(f">{'ref' if name == 'reference' else 'contig'}\n{sequence}\n")
            pysam.faidx(str(path))
            self.paths.append(path)
        self.sam = self.root / "before.sam"
        self.write_sam()

    def write_sam(self, flag=0, mapq=60, cigar="4M", sequence="TCGA", extra=""):
        self.sam.write_text("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:ref\tLN:4\n"
                            f"contig\t{flag}\tref\t1\t{mapq}\t{cigar}\t*\t0\t0\t{sequence}\t*{extra}\n")

    def run_audit(self):
        return module.audit(*self.paths, self.sam, self.root / "audit")

    def test_cli_fixed_denominator_and_damage_increment(self):
        args = [sys.executable, str(SCRIPT)]
        for name, path in zip(("before", "quality", "damage", "reference"), self.paths):
            args += [f"--{name}", str(path)]
        args += ["--alignments", str(self.sam), "--output-dir", str(self.root / "audit")]
        result = subprocess.run(args, check=True, capture_output=True, text=True)
        summary = json.loads(result.stdout)
        self.assertEqual(summary["comparisons"]["before_to_quality"]["net_reference_agreement_gain"], 1)
        self.assertEqual(summary["comparisons"]["quality_to_damage"]["net_reference_agreement_gain"], -1)
        self.assertEqual(summary["comparisons"]["before_to_damage"]["changed_sites"], 2)
        with (self.root / "audit/sites.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({int(row["position_0based"]) for row in rows}, {0, 2})
        self.assertEqual(len(rows), 4)

    def test_reverse_coordinates_and_complement(self):
        self.write_sam(flag=16, sequence="TCGA")  # TCGA is reverse-complement palindromic.
        with pysam.AlignmentFile(str(self.sam)) as sam, pysam.FastaFile(str(self.paths[3])) as ref:
            status, columns, _, strand = module.reference_columns(list(sam), "TCGA", ref)
        self.assertEqual((status, strand), ("aligned", "-"))
        self.assertEqual(columns, {3: (0, "G"), 2: (1, "G"), 1: (2, "C"), 0: (3, "A")})

    def test_soft_clips_and_insertions_are_unresolved(self):
        self.write_sam(cigar="1S1M1I1M")
        summary = self.run_audit()
        self.assertEqual(summary["comparisons"]["before_to_damage"]["classifications"], {"unaligned_position": 2})

    def test_deletion_preserves_reference_coordinates(self):
        self.write_sam(cigar="1M1D2M1S")
        with pysam.AlignmentFile(str(self.sam)) as sam, pysam.FastaFile(str(self.paths[3])) as ref:
            _, columns, _, _ = module.reference_columns(list(sam), "TCGA", ref)
        self.assertEqual(columns, {0: (0, "C"), 1: (2, "G"), 2: (3, "T")})

    def test_unknown_mapq_is_not_high_confidence(self):
        self.write_sam(mapq=255)
        summary = self.run_audit()
        self.assertEqual(summary["comparisons"]["before_to_damage"]["unresolved_sites"], 2)

    def test_secondary_and_supplementary_are_unresolved(self):
        for flag in (256, 2048):
            with self.subTest(flag=flag):
                self.write_sam(flag=flag)
                with pysam.AlignmentFile(str(self.sam)) as sam, pysam.FastaFile(str(self.paths[3])) as ref:
                    self.assertEqual(module.reference_columns(list(sam), "TCGA", ref)[0], "ambiguous_alignment")

    def test_alternative_alignment_tag_is_unresolved(self):
        self.write_sam(extra="\tXA:Z:ref,+1,4M,1;")
        summary = self.run_audit()
        self.assertEqual(summary["baseline_alignment_status"], {"ambiguous_alignment": 1})

    def test_wrong_baseline_sequence_fails(self):
        self.write_sam(sequence="AAAA")
        with self.assertRaisesRegex(ValueError, "sequence does not match"):
            self.run_audit()

    def test_missing_alignment_is_accounted_for(self):
        self.sam.write_text("@HD\tVN:1.6\n@SQ\tSN:ref\tLN:4\n")
        summary = self.run_audit()
        self.assertEqual(summary["contigs_absent_from_alignment"], 1)
        self.assertEqual(summary["comparisons"]["before_to_damage"]["unresolved_sites"], 2)

    def test_membership_mismatch_fails(self):
        self.paths[1].write_text(">wrong_id\nCCGA\n")
        pysam.faidx(str(self.paths[1]))
        with self.assertRaisesRegex(ValueError, "membership"):
            self.run_audit()

    def test_reference_ambiguity_and_two_wrong_bases(self):
        self.assertEqual(module.classify("T", "C", "N"), "ambiguous_base")
        self.assertEqual(module.classify("T", "C", "A"), "neither_matches_reference")

    def test_output_directory_cannot_be_overwritten(self):
        self.run_audit()
        with self.assertRaises(FileExistsError):
            self.run_audit()

    def alternative_fixture(self, *, flag=256, score=8, cigar="4M", limit=20, primary_mapq=0):
        self.write_sam(mapq=primary_mapq, extra="\tAS:i:8")
        with pysam.AlignmentFile(str(self.sam)) as sam:
            primary = next(sam)
            secondary = pysam.AlignedSegment.fromstring(
                f"contig\t{flag}\tref\t1\t0\t{cigar}\t*\t0\t0\t*\t*\tAS:i:{score}", sam.header)
        with pysam.FastaFile(str(self.paths[3])) as ref:
            return module.alternative_sites([primary, secondary], "TCGA", ref, {0, 2}, limit)

    def test_secondary_without_sequence_can_support_site_agreement(self):
        result = self.alternative_fixture()
        self.assertEqual(result["all_reported"][0], ("reported_reference_agreement", "C", 2))
        self.assertEqual(result["near_best_95pct"][0], ("reported_reference_agreement", "C", 2))

    def test_equal_score_reverse_alternative_disagrees_at_site(self):
        result = self.alternative_fixture(flag=272)
        self.assertEqual(result["all_reported"][0][0], "reference_base_disagreement")
        self.assertEqual(result["near_best_95pct"][0][0], "reference_base_disagreement")
        self.assertEqual(result["all_reported"][2], ("reported_reference_agreement", "G", 2))

    def test_weaker_hit_only_excluded_in_sensitivity_view(self):
        result = self.alternative_fixture(flag=272, score=4)
        self.assertEqual(result["all_reported"][0][0], "reference_base_disagreement")
        self.assertEqual(result["near_best_95pct"][0], ("reported_reference_agreement", "C", 1))

    def test_alternative_clipping_is_not_silently_ignored(self):
        result = self.alternative_fixture(cigar="1S3M")
        self.assertEqual(result["all_reported"][0][0], "alternative_does_not_cover_site")

    def test_reported_secondary_limit_stays_unresolved(self):
        result = self.alternative_fixture(limit=1)
        self.assertEqual(result["all_reported"][0][0], "secondary_limit_reached")

    def test_supplementary_stays_unresolved_in_both_views(self):
        result = self.alternative_fixture(flag=2048)
        for view in result.values():
            self.assertEqual(view[0][0], "split_or_unexpanded_alternatives")

    def test_low_mapq_singleton_is_not_rescued(self):
        self.write_sam(mapq=0, extra="\tAS:i:8")
        with pysam.AlignmentFile(str(self.sam)) as sam, pysam.FastaFile(str(self.paths[3])) as ref:
            result = module.alternative_sites(list(sam), "TCGA", ref, {0}, 20)
        self.assertEqual(result["all_reported"][0][0], "low_or_unknown_mapq")

    def test_alternative_summary_keeps_strict_counts_and_site_totals(self):
        self.write_sam(extra="\tAS:i:8")
        with self.sam.open("a") as handle:
            handle.write("contig\t256\tref\t1\t0\t4M\t*\t0\t0\t*\t*\tAS:i:8\n")
        summary = self.run_audit()
        self.assertEqual(summary["comparisons"]["before_to_damage"]["resolved_sites"], 0)
        for policy in summary["alternative_comparisons"].values():
            self.assertEqual(policy["before_to_damage"]["changed_sites"], 2)
            self.assertEqual(policy["before_to_damage"]["resolved_sites"], 2)
            self.assertEqual(policy["before_to_damage"]["net_reference_agreement_gain"], 0)
        with (self.root / "audit/alignments.tsv").open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 2)

    @unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools") and shutil.which("bash"),
                         "wrapper smoke test requires minimap2, samtools and bash")
    def test_shell_wrapper_with_synthetic_contig(self):
        rng = random.Random(60823)
        reference = "".join(rng.choice("ACGT") for _ in range(2000))
        truth = reference[500:650]
        before = list(truth)
        before[30] = "A" if truth[30] != "A" else "C"
        damage = list(truth)
        damage[60] = "A" if truth[60] != "A" else "C"
        input_dir = self.root / "input"
        input_dir.mkdir()
        for name, sequence in (("before-contigs", "".join(before)), ("quality-only", truth),
                               ("damage-consensus", "".join(damage))):
            (input_dir / f"{name}.fasta").write_text(f">contig\n{sequence}\n")
        reference_path = self.root / "synthetic-reference.fa"
        reference_path.write_text(f">ref\n{reference}\n")
        env = dict(os.environ, ANVAYA_AUDIT_PYTHON=sys.executable)
        result = subprocess.run(["bash", str(SCRIPT.with_suffix(".sh")), str(input_dir),
                                 str(reference_path), str(self.root / "wrapper")],
                                cwd=SCRIPT.parents[1], env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        summary = json.loads((self.root / "wrapper/audit/summary.json").read_text())
        self.assertEqual(summary["comparisons"]["before_to_quality"]["net_reference_agreement_gain"], 1)
        self.assertEqual(summary["comparisons"]["quality_to_damage"]["net_reference_agreement_gain"], -1)
