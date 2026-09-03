"""Tests for the progressive overlap sequence lifecycle."""

import tempfile
import unittest
from pathlib import Path

from anvaya.cli import main
from anvaya.overlap_progressive import (
    ProgressiveSequencePool,
    SequenceState,
    discover_progressive_raw_clusters,
    extend_progressive_raw_clusters,
    iterate_progressive_raw_extension,
)
from anvaya.reads import Read


class ProgressiveSequencePoolTests(unittest.TestCase):
    def test_keeps_raw_evidence_after_correction_and_consumption(self) -> None:
        raw = Read("raw-0", "ACCT", (30, 30, 30, 30))
        pool = ProgressiveSequencePool.from_reads([raw], [7])
        corrected = pool.records[0].corrected(Read("center-0", "ATCT"))
        pool = pool.replace_record(corrected.consumed())

        self.assertEqual(pool.raw_evidence, (raw,))
        self.assertEqual(pool.records[0].state, SequenceState.CONSUMED)

    def test_separates_raw_and_derived_representatives(self) -> None:
        reads = [Read("raw-0", "ACCT"), Read("raw-1", "CCTA")]
        pool = ProgressiveSequencePool.from_reads(reads, [10, 11])
        derived = pool.records[0].extended(
            Read("contig-0", "ACCTA"),
            frozenset({11}),
        )
        pool = pool.replace_record(derived)

        self.assertEqual([record.index for record in pool.active_raw], [1])
        self.assertEqual([record.index for record in pool.active_derived], [0])
        self.assertEqual(derived.contributing_molecules, frozenset({10, 11}))
        self.assertEqual(derived.generation, 1)

    def test_rejects_mutating_a_consumed_representative(self) -> None:
        record = ProgressiveSequencePool.from_reads([Read("raw", "ACCT")]).records[0]
        record = record.consumed()

        with self.assertRaises(ValueError):
            record.corrected(Read("corrected", "ATCT"))
        with self.assertRaises(ValueError):
            record.extended(Read("extended", "ACCTA"), frozenset({1}))

    def test_validates_molecule_alignment(self) -> None:
        with self.assertRaises(ValueError):
            ProgressiveSequencePool.from_reads([Read("raw", "ACCT")], [])

    def test_raw_clustering_uses_longest_center_and_stable_input_indices(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        reads = [
            Read("short-right", genome[8:20]),
            Read("center", genome[2:18]),
            Read("short-left", genome[0:12]),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)

        clusters, diagnostics = discover_progressive_raw_clusters(
            pool,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_cluster_size=3,
        )

        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0].center_index, 1)
        self.assertEqual(clusters[0].member_indices, (0, 1, 2))
        self.assertEqual(diagnostics.clustered_reads, 3)
        self.assertEqual(diagnostics.retained_unclustered_reads, 0)

    def test_undersized_cluster_does_not_hide_a_later_valid_center(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        reads = [
            Read("long-decoy", "TTTT" + genome[4:18]),
            Read("center", genome[2:18]),
            Read("left", genome[0:14]),
            Read("right", genome[6:20]),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)

        clusters, diagnostics = discover_progressive_raw_clusters(
            pool,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=10,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_cluster_size=3,
        )

        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0].center_index, 1)
        self.assertEqual(clusters[0].member_indices, (1, 2, 3))
        self.assertEqual(diagnostics.retained_unclustered_reads, 1)

    def test_clustering_excludes_derived_representatives(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        pool = ProgressiveSequencePool.from_reads(
            [Read("raw-0", genome[:14]), Read("raw-1", genome[6:])]
        )
        pool = pool.replace_record(
            pool.records[0].extended(Read("contig", genome), frozenset({1}))
        )

        clusters, diagnostics = discover_progressive_raw_clusters(
            pool,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_cluster_size=2,
        )

        self.assertEqual(clusters, ())
        self.assertEqual(diagnostics.candidate_centers, 1)
        self.assertEqual(diagnostics.retained_unclustered_reads, 1)

    def test_raw_phase_extends_center_and_consumes_only_cluster_members(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        reads = [
            Read("center", genome[2:18]),
            Read("left", genome[0:14]),
            Read("right", genome[6:20]),
            Read("unclustered", "TTTTTTTTTTTT"),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        clusters, _ = discover_progressive_raw_clusters(
            pool,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_cluster_size=3,
        )

        updated, diagnostics = extend_progressive_raw_clusters(
            pool,
            clusters,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_consensus_support=2,
            reciprocal_best_extension=False,
        )

        self.assertEqual(updated.records[0].current.sequence, genome)
        self.assertEqual(updated.records[0].state, SequenceState.EXTENDED_CONTIG)
        self.assertEqual(updated.records[1].state, SequenceState.CONSUMED)
        self.assertEqual(updated.records[2].state, SequenceState.CONSUMED)
        self.assertEqual(updated.records[3].state, SequenceState.RAW)
        self.assertEqual(updated.raw_evidence, tuple(reads))
        self.assertEqual(diagnostics.extended_centers, 1)
        self.assertEqual(diagnostics.added_bases, 4)
        self.assertEqual(diagnostics.consumed_reads, 2)

    def test_raw_phase_reindexes_extended_center_before_consuming_members(self) -> None:
        genome = "ACGTTGCACTGATCGGACCTTAGGCACTGAAAACGTT"
        reads = [
            Read("center", genome[8:24]),
            Read("left", genome[4:20]),
            Read("right", genome[12:28]),
            Read("later-right-1", genome[20:36]),
            Read("later-right-2", genome[20:36]),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        clusters, _ = discover_progressive_raw_clusters(
            pool,
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_cluster_size=3,
        )

        updated, diagnostics = extend_progressive_raw_clusters(
            pool,
            clusters,
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_consensus_support=2,
            reciprocal_best_extension=False,
            maximum_rounds=3,
        )

        self.assertEqual(updated.records[0].current.sequence, genome[4:36])
        self.assertEqual(updated.records[3].state, SequenceState.CONSUMED)
        self.assertEqual(updated.records[4].state, SequenceState.CONSUMED)
        self.assertEqual(diagnostics.extension_rounds, 2)
        self.assertEqual(diagnostics.recruited_reads, 2)

    def test_cli_progressive_projection_preserves_primary_output(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads_path = root / "reads.fasta"
            baseline_path = root / "baseline.fasta"
            audited_path = root / "audited.fasta"
            projection_path = root / "progressive.fasta"
            reads_path.write_text(
                "".join(
                    (
                        f">center\n{genome[2:18]}\n",
                        f">left\n{genome[0:14]}\n",
                        f">right\n{genome[6:20]}\n",
                    )
                ),
                encoding="utf-8",
            )
            common = [
                "overlap-assemble",
                "-i", str(reads_path),
                "--anchor-k", "3",
                "--anchors-per-read", "8",
                "--min-anchor-matches", "1",
                "--min-overlap", "6",
                "--min-cluster-size", "3",
                "--max-contig-iterations", "0",
                "--ranked-extension",
                "--reciprocal-best-extension",
                "--damage-aware-ranking",
                "--min-output-length", "1",
            ]

            self.assertEqual(main(common + ["-o", str(baseline_path)]), 0)
            self.assertEqual(
                main(
                    common
                    + [
                        "-o", str(audited_path),
                        "--progressive-raw-phase-audit",
                        "--progressive-raw-phase-projection",
                        str(projection_path),
                    ]
                ),
                0,
            )

            self.assertEqual(audited_path.read_bytes(), baseline_path.read_bytes())
            self.assertTrue(projection_path.exists())

    def test_persistent_center_recruits_raw_reads_after_initial_extension(self) -> None:
        genome = "ACGTTGCACTGATCGGACCTTAGGCACTGAAAAC"
        reads = [
            Read("center", genome[4:20]),
            Read("left", genome[0:16]),
            Read("right", genome[8:24]),
            Read("later-right-1", genome[16:32]),
            Read("later-right-2", genome[18:34]),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        clusters, _ = discover_progressive_raw_clusters(
            pool,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_cluster_size=3,
        )
        pool, _ = extend_progressive_raw_clusters(
            pool,
            clusters,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_consensus_support=2,
            reciprocal_best_extension=False,
            maximum_rounds=1,
        )

        pool, diagnostics = iterate_progressive_raw_extension(
            pool,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_consensus_support=2,
            maximum_iterations=3,
        )

        self.assertEqual(pool.records[0].current.sequence, genome[:32])
        self.assertEqual(pool.records[3].state, SequenceState.CONSUMED)
        self.assertEqual(pool.records[4].state, SequenceState.CONSUMED)
        self.assertEqual(diagnostics.extended_centers, 1)
        self.assertGreater(diagnostics.added_bases, 0)
        self.assertTrue(diagnostics.converged)


if __name__ == "__main__":
    unittest.main()
