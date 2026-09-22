import random
import unittest
from unittest.mock import Mock, patch

from anvaya.overlap_adaptive_rescue import (
    _retain_primary_anchored_edges,
    audit_adaptive_support_rescue,
    project_high_confidence_support_two_rescue,
    project_selective_support_rescue,
)
from anvaya.overlap_assembly import _MasterOverlapEdge
from anvaya.overlap_progressive import ProgressiveSequencePool
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


class AdaptiveSupportRescueTests(unittest.TestCase):
    def test_support_two_quality_controls_are_independent(self) -> None:
        rng = random.Random(971)
        sequence = "".join(rng.choice("ACGT") for _ in range(140))
        for seed_filter in (False, True):
            for recruit_filter in (False, True):
                for bad_seed in (False, True):
                    with self.subTest(seed=seed_filter, recruit=recruit_filter, bad_seed=bad_seed):
                        pool = ProgressiveSequencePool.from_reads([
                            Read("seed", sequence[:60], (0 if bad_seed else 30,) * 40 + (30,) * 20),
                            Read("partner", sequence[40:100], (30,) * 60),
                            Read("recruit", sequence[80:], (30,) * 20 + (0,) * 40),
                        ])
                        projected, diagnostics = project_high_confidence_support_two_rescue(
                            pool, [], anchor_k=5, anchors_per_read=50, minimum_overlap=20,
                            filter_seed_quality=seed_filter,
                            filter_recruit_quality=recruit_filter,
                        )
                        rejected_seed = seed_filter and bad_seed
                        expected = [] if rejected_seed else [sequence[:100] if recruit_filter else sequence]
                        self.assertEqual([read.sequence for read in projected], expected)
                        self.assertEqual(diagnostics.seed_quality_rejections, int(rejected_seed))
                        self.assertEqual(diagnostics.seed_quality_rejected_bases, 60 if rejected_seed else 0)
                        rejected_recruit = recruit_filter and not rejected_seed
                        self.assertEqual(diagnostics.recruit_quality_rejections, int(rejected_recruit))
                        self.assertEqual(diagnostics.recruit_quality_rejected_overhang_bases, 40 if rejected_recruit else 0)

    def test_support_two_rejects_unsupported_seed_without_quality(self) -> None:
        sequence = "ACGTTGCACTGATCGGACCTAGTA"
        for qualities in (None, (0,) * 6 + (30,) * 12):
            with self.subTest(qualities=qualities):
                pool = ProgressiveSequencePool.from_reads([
                    Read("seed", sequence[:18], qualities),
                    Read("partner", sequence[6:], (30,) * 18),
                ])
                projected, diagnostics = project_high_confidence_support_two_rescue(
                    pool, [], anchor_k=3, anchors_per_read=16, minimum_overlap=8,
                    filter_seed_quality=True,
                )
                self.assertEqual(projected, [])
                self.assertEqual(diagnostics.boundary_quality_rejections, 1)

    def test_support_two_checks_later_recruit_overhang_quality(self) -> None:
        rng = random.Random(971)
        sequence = "".join(rng.choice("ACGT") for _ in range(140))
        for quality, reverse in ((q, r) for q in (0, 19, 20, 30, None) for r in (False, True)):
            with self.subTest(quality=quality, reverse=reverse):
                qualities = None if quality is None else (30,) * 20 + (quality,) * 40
                recruit = sequence[80:]
                if reverse:
                    recruit = reverse_complement(recruit)
                    qualities = None if qualities is None else qualities[::-1]
                pool = ProgressiveSequencePool.from_reads([
                    Read("seed", sequence[:60], (30,) * 60),
                    Read("partner", sequence[40:100], (30,) * 60),
                    Read("recruit", recruit, qualities),
                ])
                projected, _ = project_high_confidence_support_two_rescue(
                    pool, [], anchor_k=5, anchors_per_read=50, minimum_overlap=20,
                )
                expected = sequence if quality is not None and quality >= 20 else sequence[:100]
                self.assertEqual([read.sequence for read in projected], [expected])

    def test_support_two_allows_low_quality_seed_bases_inside_overlap(self) -> None:
        sequence = "ACGTTGCACTGATCGGACCTAGTA"
        pool = ProgressiveSequencePool.from_reads([
            Read("seed", sequence[:18], (30,) * 6 + (0,) * 12),
            Read("partner", sequence[6:], (30,) * 18),
        ])
        projected, _ = project_high_confidence_support_two_rescue(
            pool, [], anchor_k=3, anchors_per_read=16, minimum_overlap=8,
        )
        self.assertEqual([read.sequence for read in projected], [sequence])

    def test_support_two_rescue_accepts_high_quality_exact_overlap(self) -> None:
        primary = Read("primary", "GGGGAAAACCCCTTTT")
        orphan = "ACGTTGCACTGATCGGACCTAGTA"
        reads = [
            Read("primary", primary.sequence, (30,) * len(primary.sequence)),
            Read("orphan-left", orphan[:18], (30,) * 18),
            Read("orphan-right", orphan[6:], (30,) * 18),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(pool.records[0].corrected(primary))

        projected, diagnostics = project_high_confidence_support_two_rescue(
            pool,
            [primary],
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=2,
            minimum_overlap=8,
        )

        self.assertIn(orphan, [item.sequence for item in projected])
        self.assertEqual(diagnostics.admitted_clusters, 1)
        self.assertEqual(diagnostics.novel_contigs, 1)

    def test_support_two_rescue_rejects_low_quality_extension(self) -> None:
        primary = Read("primary", "GGGGAAAACCCCTTTT")
        orphan = "ACGTTGCACTGATCGGACCTAGTA"
        reads = [
            Read("primary", primary.sequence, (30,) * len(primary.sequence)),
            Read("orphan-left", orphan[:18], (30,) * 18),
            Read("orphan-right", orphan[6:], (30,) * 12 + (10,) * 6),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(pool.records[0].corrected(primary))

        projected, diagnostics = project_high_confidence_support_two_rescue(
            pool,
            [primary],
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=2,
            minimum_overlap=8,
        )

        self.assertEqual([item.sequence for item in projected], [primary.sequence])
        self.assertEqual(diagnostics.boundary_quality_rejections, 1)

    def test_selective_rescue_keeps_unattached_support_three_cluster(self) -> None:
        primary = "ACGTTGCACTGATCGGACCT"
        orphan = "TTTTCCCCAAAAGGGGTTAA"
        reads = [
            Read("primary", primary),
            Read("orphan-center", orphan[1:19]),
            Read("orphan-left", orphan[:18]),
            Read("orphan-right", orphan[2:]),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(
            pool.records[0].corrected(Read("primary", primary))
        )

        projected, diagnostics = project_selective_support_rescue(
            pool,
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=1,
            minimum_overlap=8,
        )

        self.assertEqual(
            [item.sequence for item in projected],
            [primary, orphan[1:19]],
        )
        self.assertEqual(diagnostics.novel_contigs, 1)
        self.assertEqual(diagnostics.novel_bases, len(orphan[1:19]))
        self.assertEqual(diagnostics.insufficient_molecule_support, 0)

    def test_keeps_only_rescue_chains_anchored_to_primary(self) -> None:
        edges = {
            ((0, False), (1, False)): _MasterOverlapEdge(
                (0, False), (1, False), 5, 10
            ),
            ((1, False), (2, False)): _MasterOverlapEdge(
                (1, False), (2, False), 5, 10
            ),
            ((3, False), (4, False)): _MasterOverlapEdge(
                (3, False), (4, False), 5, 10
            ),
        }

        retained, anchored, rejected = _retain_primary_anchored_edges(edges, 1)

        self.assertEqual(
            set(retained),
            {((0, False), (1, False)), ((1, False), (2, False))},
        )
        self.assertEqual(anchored, 1)
        self.assertEqual(rejected, 1)

    def test_spells_rescue_chain_only_when_anchored_to_primary(self) -> None:
        genome = "ACGTTGCACTGATCGGACCTAGTATTCGAG"
        primary = genome[:14]
        rescue_pool = ProgressiveSequencePool.from_reads(
            [Read("rescue-1", genome[6:22]), Read("rescue-2", genome[14:30])]
        )
        for index in range(2):
            rescue_pool = rescue_pool.replace_record(
                rescue_pool.records[index].corrected(
                    Read(
                        f"rescue-{index + 1}",
                        rescue_pool.records[index].raw.sequence,
                    )
                )
            )
        pool = ProgressiveSequencePool.from_reads(
            [Read("primary", primary), Read("eligible", "TTTTAAAACCCC")]
        )
        pool = pool.replace_record(
            pool.records[0].corrected(Read("primary", primary))
        )

        with (
            patch(
                "anvaya.overlap_adaptive_rescue.discover_progressive_raw_clusters",
                return_value=([], Mock(clusters=2, clustered_reads=6)),
            ),
            patch(
                "anvaya.overlap_adaptive_rescue.extend_progressive_raw_clusters",
                return_value=(rescue_pool, None),
            ),
        ):
            projected, diagnostics = audit_adaptive_support_rescue(
                pool,
                anchor_k=3,
                anchors_per_read=16,
                minimum_anchor_matches=1,
                minimum_overlap=8,
                allow_rescue_chains=True,
            )

        self.assertEqual([item.sequence for item in projected], [genome[:30]])
        self.assertGreaterEqual(diagnostics.candidate_rescue_links, 1)
        self.assertEqual(diagnostics.anchored_components, 1)

    def test_promotes_only_rescue_contig_attached_to_primary(self) -> None:
        genome = "ACGTTGCACTGATCGGACCTAGTATTCGAGCA"
        reads = [
            Read("primary-center", genome[2:20]),
            Read("primary-left", genome[:18]),
            Read("primary-right", genome[4:22]),
            Read("rescue-center", genome[12:30]),
            Read("rescue-left", genome[10:28]),
            Read("rescue-right", genome[14:32]),
            Read("orphan-center", "TTTTCCCCAAAAGGGGTT"),
            Read("orphan-left", "AATTTTCCCCAAAAGGGG"),
            Read("orphan-right", "TTCCCCAAAAGGGGTTAA"),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        primary = pool.records[0].extended(
            Read("primary", genome[:22]), frozenset({0, 1, 2})
        )
        pool = pool.replace_record(primary)
        pool = pool.replace_record(pool.records[1].consumed())
        pool = pool.replace_record(pool.records[2].consumed())

        projected, diagnostics = audit_adaptive_support_rescue(
            pool,
            minimum_rescue_support=3,
            anchor_k=3,
            anchors_per_read=16,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=1,
            minimum_overlap=8,
        )

        self.assertTrue(
            any(
                item.sequence.startswith(genome[:22])
                and len(item.sequence) > 22
                for item in projected
            )
        )
        self.assertEqual(diagnostics.promoted_rescue_contigs, 1)
        self.assertEqual(diagnostics.rejected_unattached_contigs, 1)
        self.assertEqual(diagnostics.reciprocal_edges, 1)

    def test_does_not_promote_an_unattached_support_three_cluster(self) -> None:
        primary = "ACGTTGCACTGATCGGACCT"
        orphan = "TTTTCCCCAAAAGGGGTTAA"
        reads = [
            Read("primary", primary),
            Read("orphan-center", orphan[1:19]),
            Read("orphan-left", orphan[:18]),
            Read("orphan-right", orphan[2:]),
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(
            pool.records[0].corrected(Read("primary", primary))
        )

        projected, diagnostics = audit_adaptive_support_rescue(
            pool,
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=1,
            minimum_overlap=8,
        )

        self.assertEqual([item.sequence for item in projected], [primary])
        self.assertEqual(diagnostics.promoted_rescue_contigs, 0)
        self.assertEqual(diagnostics.rejected_unattached_contigs, 1)

    def test_accepts_raw_confirmed_near_exact_rescue_link(self) -> None:
        genome = "ACGTTGCACTGATCGGACCTAGTATTCGAGCA"
        damaged_primary = genome[:15] + "A" + genome[16:22]
        reads = [
            Read("primary", damaged_primary),
            Read("rescue-center", genome[12:30]),
            Read("rescue-left", genome[10:28]),
            Read("rescue-right", genome[14:32]),
        ] + [
            Read(f"evidence-{index}", genome, [30] * len(genome))
            for index in range(3)
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(
            pool.records[0].corrected(Read("primary", damaged_primary))
        )
        for index in range(4, 7):
            pool = pool.replace_record(pool.records[index].consumed())

        projected, diagnostics = audit_adaptive_support_rescue(
            pool,
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=1,
            minimum_overlap=8,
            allow_near_exact=True,
            near_exact_minimum_identity=0.80,
            damage_end_window=0,
        )

        self.assertEqual([item.sequence for item in projected], [genome[:30]])
        self.assertEqual(diagnostics.near_exact_accepted, 1)
        self.assertEqual(diagnostics.corrected_overlap_bases, 1)

    def test_rejects_near_exact_rescue_link_with_strain_conflict(self) -> None:
        genome = "ACGTTGCACTGATCGGACCTAGTATTCGAGCA"
        alternate = genome[:15] + "A" + genome[16:]
        reads = [
            Read("primary", alternate[:22]),
            Read("rescue-center", genome[12:30]),
            Read("rescue-left", genome[10:28]),
            Read("rescue-right", genome[14:32]),
        ] + [
            Read(f"primary-{index}", genome, [30] * len(genome))
            for index in range(2)
        ] + [
            Read(f"alternate-{index}", alternate, [30] * len(alternate))
            for index in range(2)
        ]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(
            pool.records[0].corrected(Read("primary", alternate[:22]))
        )
        for index in range(4, 8):
            pool = pool.replace_record(pool.records[index].consumed())

        projected, diagnostics = audit_adaptive_support_rescue(
            pool,
            anchor_k=3,
            anchors_per_read=16,
            minimum_anchor_matches=1,
            minimum_overlap=8,
            allow_near_exact=True,
            near_exact_minimum_identity=0.80,
            damage_end_window=0,
        )

        self.assertEqual([item.sequence for item in projected], [alternate[:22]])
        self.assertGreaterEqual(diagnostics.strain_conflicts, 1)
        self.assertEqual(diagnostics.near_exact_accepted, 0)


if __name__ == "__main__":
    unittest.main()
