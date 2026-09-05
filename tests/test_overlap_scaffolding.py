import io
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.cli import main
from anvaya.overlap_scaffolding import scaffold_progressive_contigs
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


class ProgressivePairedScaffoldingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.first = "ACGTTGCACTGATCGGACCTAGTA"
        self.second = "TTACCGATGCTAGGCTAACGTCGA"
        self.left = self.first[12:20]
        self.right = reverse_complement(self.second[2:10])

    def test_projects_a_reciprocal_supported_pair_link(self) -> None:
        contigs = [Read("first", self.first), Read("second", self.second)]
        left_reads = [Read(f"pair-{index}/1", self.left) for index in range(3)]
        right_reads = [Read(f"pair-{index}/2", self.right) for index in range(3)]

        projected, diagnostics = scaffold_progressive_contigs(
            contigs,
            left_reads,
            right_reads,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_mapping_identity=1.0,
            minimum_mapping_ry_identity=1.0,
            minimum_link_support=3,
            fragment_mean=30,
            fragment_sd=2,
        )

        self.assertEqual(len(projected), 1)
        self.assertEqual(len(projected[0].sequence), 56)
        self.assertIn("N" * 8, projected[0].sequence)
        self.assertEqual(diagnostics.uniquely_mapped_pairs, 3)
        self.assertEqual(diagnostics.reciprocal_links, 1)
        self.assertEqual(diagnostics.linear_paths, 1)
        self.assertEqual(diagnostics.joined_contigs, 2)
        self.assertEqual(diagnostics.inserted_gap_bases, 8)

    def test_does_not_use_ambiguously_mapped_mates(self) -> None:
        contigs = [
            Read("first-copy-1", self.first),
            Read("first-copy-2", self.first),
            Read("second", self.second),
        ]
        left_reads = [Read(f"pair-{index}/1", self.left) for index in range(3)]
        right_reads = [Read(f"pair-{index}/2", self.right) for index in range(3)]

        projected, diagnostics = scaffold_progressive_contigs(
            contigs,
            left_reads,
            right_reads,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_mapping_identity=1.0,
            minimum_mapping_ry_identity=1.0,
            minimum_link_support=3,
            fragment_mean=30,
            fragment_sd=2,
        )

        self.assertEqual(projected, contigs)
        self.assertEqual(diagnostics.uniquely_mapped_pairs, 0)
        self.assertEqual(diagnostics.ambiguous_reads, 3)
        self.assertEqual(diagnostics.reciprocal_links, 0)

    def test_abstains_when_fragment_distribution_cannot_be_calibrated(self) -> None:
        contigs = [Read("first", self.first), Read("second", self.second)]
        left_reads = [Read("pair/1", self.left)]
        right_reads = [Read("pair/2", self.right)]

        projected, diagnostics = scaffold_progressive_contigs(
            contigs,
            left_reads,
            right_reads,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_mapping_identity=1.0,
            minimum_mapping_ry_identity=1.0,
            minimum_calibration_pairs=2,
        )

        self.assertEqual(projected, contigs)
        self.assertEqual(diagnostics.calibration_pairs, 0)
        self.assertEqual(diagnostics.fragment_mean, 0.0)
        self.assertEqual(diagnostics.reciprocal_links, 0)

    def test_rejects_fragment_calibration_shorter_than_typical_mates(self) -> None:
        contig = "ACGTTGCACTGATCGGACCTAGTATTCGAGCACTGA"
        short_left = Read("short/1", contig[2:10])
        short_right = Read("short/2", reverse_complement(contig[6:14]))
        long_left = Read("long/1", "AAAAAAAAAAAAAAAAAAAA")
        long_right = Read("long/2", "CCCCCCCCCCCCCCCCCCCC")

        projected, diagnostics = scaffold_progressive_contigs(
            [Read("contig", contig)],
            [short_left, long_left],
            [short_right, long_right],
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_mapping_identity=1.0,
            minimum_mapping_ry_identity=1.0,
            minimum_calibration_pairs=1,
        )

        self.assertEqual(projected, [Read("contig", contig)])
        self.assertEqual(diagnostics.calibration_pairs, 1)
        self.assertEqual(diagnostics.fragment_mean, 12.0)
        self.assertFalse(diagnostics.calibration_valid)
        self.assertEqual(diagnostics.reciprocal_links, 0)

    def test_cli_accepts_paired_input_and_writes_separate_projection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left_path = root / "left.fastq"
            right_path = root / "right.fastq"
            primary_path = root / "contigs.fasta"
            progressive_path = root / "progressive.fasta"
            scaffold_path = root / "scaffolds.fasta"
            left_path.write_text(
                "".join(
                    f"@pair-{index}/1\n{self.left}\n+\nIIIIIIII\n"
                    for index in range(3)
                ),
                encoding="utf-8",
            )
            right_path.write_text(
                "".join(
                    f"@pair-{index}/2\n{self.right}\n+\nIIIIIIII\n"
                    for index in range(3)
                ),
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                result = main(
                    [
                        "overlap-assemble",
                        "-1", str(left_path),
                        "-2", str(right_path),
                        "--anchor-k", "3",
                        "--anchors-per-read", "8",
                        "--min-anchor-matches", "1",
                        "--min-overlap", "6",
                        "--min-cluster-size", "3",
                        "--ranked-extension",
                        "--reciprocal-best-extension",
                        "--damage-aware-ranking",
                        "--max-contig-iterations", "0",
                        "--progressive-raw-phase-audit",
                        "--progressive-raw-phase-projection",
                        str(progressive_path),
                        "--paired-progressive-scaffold-audit",
                        "--paired-progressive-scaffold-projection",
                        str(scaffold_path),
                        "--paired-fragment-mean", "30",
                        "--paired-fragment-sd", "2",
                        "--min-output-length", "1",
                        "-o", str(primary_path),
                    ]
                )

            self.assertEqual(result, 0)
            self.assertTrue(primary_path.exists())
            self.assertTrue(progressive_path.exists())
            self.assertTrue(scaffold_path.exists())
            self.assertIn(
                "overlap_paired_scaffold_audit=true", stdout.getvalue()
            )
            self.assertIn(
                "overlap_paired_scaffold_reciprocal_links=1",
                stdout.getvalue(),
            )


if __name__ == "__main__":
    unittest.main()
