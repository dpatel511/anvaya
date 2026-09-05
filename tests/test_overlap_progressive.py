"""Tests for the progressive overlap sequence lifecycle."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from anvaya.cli import main
from anvaya.overlap_assembly import _Alignment
from anvaya.overlap_progressive import (
    _RejectedExtensionDiagnostics,
    _audit_rejected_extension_sides,
    _evidence_priority_key,
    _filter_strict_boundary_evidence,
    ProgressiveSequencePool,
    SequenceState,
    discover_progressive_raw_clusters,
    extend_progressive_raw_clusters,
    iterate_progressive_raw_extension,
)
from anvaya.overlap_progressive_links import (
    audit_raw_supported_progressive_links,
)
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


class ProgressiveSequencePoolTests(unittest.TestCase):
    def test_strict_boundary_filter_requires_unanimous_q20_support(self) -> None:
        target = "ACGTCCGT"
        reads = [
            Read(f"support-{index}", "CCGTTT", (30,) * 6)
            for index in range(3)
        ]
        candidates = [
            _Alignment(index, read.sequence, 4, 1)
            for index, read in enumerate(reads)
        ]

        accepted, sides, conflicts, quality = _filter_strict_boundary_evidence(
            target,
            candidates,
            candidates[:1],
            reads,
            minimum_support=3,
            minimum_base_quality=20,
        )

        self.assertEqual(accepted, candidates[:1])
        self.assertEqual((sides, conflicts, quality), (1, 0, 0))

        conflicting_reads = reads + [Read("conflict", "CCGTGG", (30,) * 6)]
        conflicting_candidates = candidates + [
            _Alignment(3, "CCGTGG", 4, 1)
        ]
        accepted, _, conflicts, _ = _filter_strict_boundary_evidence(
            target,
            conflicting_candidates,
            candidates[:1],
            conflicting_reads,
            minimum_support=3,
            minimum_base_quality=20,
        )

        self.assertEqual(accepted, [])
        self.assertEqual(conflicts, 1)

    def test_failure_audit_classifies_support_conflict_quality_and_damage(self) -> None:
        target = "ACGTCCGT"
        reads = [
            Read("support-1", "CCGTTT", (30,) * 6),
            Read("support-2", "CCGTTT", (30,) * 6),
            Read("conflict", "TCGTGG", (10,) * 6),
        ]
        candidates = [
            _Alignment(index, read.sequence, 4, 1)
            for index, read in enumerate(reads)
        ]
        diagnostics = _RejectedExtensionDiagnostics()

        _audit_rejected_extension_sides(
            target,
            candidates,
            reads,
            diagnostics,
            minimum_base_quality=20,
            damage_end_window=5,
        )

        self.assertEqual(diagnostics.right_sides, 1)
        self.assertEqual(diagnostics.support_2, 1)
        self.assertEqual(diagnostics.conflicting_sides, 1)
        self.assertEqual(diagnostics.extension_1_5, 1)
        self.assertEqual(diagnostics.high_quality_boundary_observations, 2)
        self.assertEqual(diagnostics.low_quality_boundary_observations, 1)
        self.assertEqual(diagnostics.damage_compatible_mismatches, 1)

    def test_evidence_priority_prefers_raw_support_before_length(self) -> None:
        pool = ProgressiveSequencePool.from_reads(
            [Read("long", "A" * 20), Read("supported", "C" * 12)]
        )
        long = pool.records[0].corrected(Read("long", "A" * 20))
        supported = pool.records[1].corrected(Read("supported", "C" * 12))
        candidates = [
            _Alignment(2, "GG" + supported.current.sequence, -2, 1),
            _Alignment(3, "TT" + supported.current.sequence, -2, 1),
        ]

        self.assertGreater(
            _evidence_priority_key(supported, candidates),
            _evidence_priority_key(long, []),
        )

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

    def test_raw_cluster_support_counts_independent_molecules(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        reads = [
            Read("pair/1", genome[2:18]),
            Read("pair/2", genome[2:18]),
            Read("independent", genome[4:20]),
        ]
        pool = ProgressiveSequencePool.from_reads(reads, [0, 0, 1])

        clusters, diagnostics = discover_progressive_raw_clusters(
            pool,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=8,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_cluster_size=3,
        )

        self.assertEqual(clusters, ())
        self.assertEqual(diagnostics.clustered_reads, 0)

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
            support_three_projection_path = root / "support-three.fasta"
            selective_rescue_path = root / "selective-rescue.fasta"
            support_two_rescue_path = root / "support-two-rescue.fasta"
            support_two_master_path = root / "support-two-master.fasta"
            selective_master_path = root / "selective-master.fasta"
            raw_selective_master_path = root / "raw-selective-master.fasta"
            link_projection_path = root / "progressive-links.fasta"
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
                        "--support-three-progressive-extension-audit",
                        "--support-three-progressive-projection",
                        str(support_three_projection_path),
                        "--selective-support-rescue-audit",
                        "--selective-support-rescue-projection",
                        str(selective_rescue_path),
                        "--high-confidence-support-two-rescue-audit",
                        "--high-confidence-support-two-rescue-projection",
                        str(support_two_rescue_path),
                        "--support-two-master-graph-audit",
                        "--support-two-master-graph-projection",
                        str(support_two_master_path),
                        "--selective-rescue-master-graph-audit",
                        "--selective-rescue-master-graph-projection",
                        str(selective_master_path),
                        "--raw-confirmed-selective-master-graph-audit",
                        "--raw-confirmed-selective-master-graph-projection",
                        str(raw_selective_master_path),
                        "--raw-supported-progressive-link-audit",
                        "--raw-confirmed-progressive-links",
                        "--raw-supported-progressive-link-projection",
                        str(link_projection_path),
                    ]
                ),
                0,
            )

            self.assertEqual(audited_path.read_bytes(), baseline_path.read_bytes())
            self.assertTrue(projection_path.exists())
            self.assertTrue(support_three_projection_path.exists())
            self.assertTrue(selective_rescue_path.exists())
            self.assertTrue(support_two_rescue_path.exists())
            self.assertTrue(support_two_master_path.exists())
            self.assertTrue(selective_master_path.exists())
            self.assertTrue(raw_selective_master_path.exists())
            self.assertTrue(link_projection_path.exists())

    def test_progressive_link_requires_two_unique_spanning_molecules(self) -> None:
        left = "ACGTTGCACTGA"
        right = "CTGATCGGACCT"
        merged = "ACGTTGCACTGATCGGACCT"
        reads = [
            Read("left", left),
            Read("right", right),
            Read("bridge-1", "CACTGATCGG"),
            Read("bridge-2", "CACTGATCGG"),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(
            pool.records[0].corrected(Read("left-contig", left))
        )
        pool = pool.replace_record(
            pool.records[1].corrected(Read("right-contig", right))
        )

        projection, diagnostics = audit_raw_supported_progressive_links(
            pool,
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=1,
            minimum_overlap=4,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_read_support=2,
        )

        self.assertEqual(len(projection), 1)
        self.assertEqual(
            min(projection[0].sequence, reverse_complement(projection[0].sequence)),
            min(merged, reverse_complement(merged)),
        )
        self.assertEqual(diagnostics.supported_dovetails, 1)
        self.assertEqual(diagnostics.linear_paths, 1)
        self.assertEqual(diagnostics.merged_contigs, 2)

    def test_progressive_link_abstains_without_enough_molecules(self) -> None:
        left = "ACGTTGCACTGA"
        right = "CTGATCGGACCT"
        reads = [
            Read("left", left),
            Read("right", right),
            Read("bridge", "CACTGATCGG"),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(
            pool.records[0].corrected(Read("left-contig", left))
        )
        pool = pool.replace_record(
            pool.records[1].corrected(Read("right-contig", right))
        )

        projection, diagnostics = audit_raw_supported_progressive_links(
            pool,
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=1,
            minimum_overlap=4,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_read_support=2,
        )

        self.assertEqual([read.sequence for read in projection], [left, right])
        self.assertEqual(diagnostics.support_at_least_1, 1)
        self.assertEqual(diagnostics.supported_dovetails, 0)

    def test_progressive_link_rejects_molecules_supporting_competing_edges(self) -> None:
        sequences = [
            "ACGTTGCACTGA",
            "CTGATCGGACCT",
            "CTGATTCCGGAA",
        ]
        pool = ProgressiveSequencePool.from_reads(
            [Read(f"raw-{index}", sequence) for index, sequence in enumerate(sequences)]
        )
        for index, sequence in enumerate(sequences):
            pool = pool.replace_record(
                pool.records[index].corrected(Read(f"contig-{index}", sequence))
            )

        with patch(
            "anvaya.overlap_progressive_links._spanning_molecules",
            return_value={100, 101},
        ):
            projection, diagnostics = audit_raw_supported_progressive_links(
                pool,
                anchor_k=3,
                anchors_per_read=16,
                minimum_anchor_matches=1,
                minimum_overlap=4,
                minimum_identity=1.0,
                minimum_ry_identity=1.0,
                minimum_read_support=2,
            )

        self.assertEqual([read.sequence for read in projection], sequences)
        self.assertEqual(diagnostics.ambiguous_bridge_molecules, 2)
        self.assertEqual(diagnostics.supported_dovetails, 0)

    def test_progressive_link_accepts_raw_confirmed_near_exact_overlap(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        first = genome[:14]
        damaged_first = first[:10] + "A" + first[11:]
        second = genome[6:]
        reads = [
            Read(f"raw-{index}", genome, [30] * len(genome))
            for index in range(3)
        ] + [
            Read("first-raw", damaged_first, [30] * len(damaged_first)),
            Read("second-raw", second, [30] * len(second)),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(
            pool.records[3].corrected(Read("first-contig", damaged_first))
        )
        pool = pool.replace_record(
            pool.records[4].corrected(Read("second-contig", second))
        )

        projection, diagnostics = audit_raw_supported_progressive_links(
            pool,
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            allow_near_exact=True,
            near_exact_minimum_identity=0.80,
            damage_end_window=0,
        )

        self.assertEqual([contig.sequence for contig in projection], [genome])
        self.assertEqual(diagnostics.near_exact_accepted, 1)
        self.assertEqual(diagnostics.corrected_overlap_bases, 1)

    def test_progressive_link_rejects_raw_supported_strain_conflict(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        first = genome[:14]
        alternate = first[:10] + "A" + first[11:]
        second = genome[6:]
        alternate_genome = genome[:10] + "A" + genome[11:]
        reads = [
            Read(f"primary-{index}", genome, [30] * len(genome))
            for index in range(2)
        ] + [
            Read(f"alternate-{index}", alternate_genome, [30] * len(genome))
            for index in range(2)
        ] + [
            Read("first-raw", alternate, [30] * len(alternate)),
            Read("second-raw", second, [30] * len(second)),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(
            pool.records[4].corrected(Read("first-contig", alternate))
        )
        pool = pool.replace_record(
            pool.records[5].corrected(Read("second-contig", second))
        )

        projection, diagnostics = audit_raw_supported_progressive_links(
            pool,
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            allow_near_exact=True,
            near_exact_minimum_identity=0.80,
            damage_end_window=0,
        )

        self.assertCountEqual(
            [contig.sequence for contig in projection],
            [alternate, second],
        )
        self.assertGreaterEqual(diagnostics.strain_conflicts, 1)
        self.assertEqual(diagnostics.near_exact_accepted, 0)

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
            evidence_priority=True,
        )

        self.assertEqual(pool.records[0].current.sequence, genome[:32])
        self.assertEqual(pool.records[3].state, SequenceState.CONSUMED)
        self.assertEqual(pool.records[4].state, SequenceState.CONSUMED)
        self.assertEqual(diagnostics.extended_centers, 1)
        self.assertGreater(diagnostics.added_bases, 0)
        self.assertTrue(diagnostics.converged)
        self.assertGreater(diagnostics.evidence_priority_sweeps, 0)

    def test_support_three_shadow_extends_with_unanimous_q20_reads(self) -> None:
        genome = "ACGTTGCACTGATCGGACCTTAGGCACTGAAAAC"
        sequences = [
            genome[4:20],
            genome[0:16],
            genome[8:24],
            genome[18:36],
            genome[18:36],
            genome[18:36],
        ]
        reads = [
            Read(f"read-{index}", sequence, (30,) * len(sequence))
            for index, sequence in enumerate(sequences)
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(
            pool.records[0].extended(
                Read("primary", genome[:24]),
                frozenset({1, 2}),
            )
        )
        pool = pool.replace_record(pool.records[1].consumed())
        pool = pool.replace_record(pool.records[2].consumed())

        shadow, diagnostics = iterate_progressive_raw_extension(
            pool,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=1.0,
            minimum_ry_identity=1.0,
            minimum_consensus_support=3,
            maximum_iterations=2,
            require_strict_boundary_evidence=True,
        )

        self.assertEqual(shadow.records[0].current.sequence, genome[:36])
        self.assertEqual(diagnostics.extended_centers, 1)
        self.assertGreaterEqual(diagnostics.strict_boundary_accepted_sides, 1)
        self.assertGreater(diagnostics.added_bases, 0)


if __name__ == "__main__":
    unittest.main()
