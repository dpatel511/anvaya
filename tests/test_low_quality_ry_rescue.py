"""Synthetic recruitment controls; all positions and offsets are 0-based."""

import random
import unittest

from anvaya.overlap_assembly import _CandidateDiagnostics, _candidate_alignments
from anvaya.overlap_index import _anchor_index
from anvaya.overlap_progressive import ProgressiveSequencePool, iterate_progressive_raw_extension
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


class LowQualityRYRescueTests(unittest.TestCase):
    def setUp(self):
        rng = random.Random(81203)
        self.truth = "".join(rng.choice("ACGT") for _ in range(80))

    def partner(self, positions=(15,), quality=15, reverse=False):
        sequence = list(self.truth[30:])
        qualities = [35] * len(sequence)
        for position in positions:
            sequence[position] = "A" if sequence[position] in "CT" else "C"
            qualities[position] = quality
        sequence = "".join(sequence)
        if reverse:
            sequence = reverse_complement(sequence)
            qualities.reverse()
        return Read("partner", sequence, tuple(qualities))

    def search(self, read, enabled=True, **overrides):
        options = dict(anchor_k=5, anchors_per_read=60, maximum_anchor_occurrences=100,
                       minimum_anchor_matches=1, minimum_overlap=30, minimum_identity=.9,
                       minimum_ry_identity=.99, position_bits=6, target_window=60)
        options.update(overrides)
        diagnostics = _CandidateDiagnostics()
        result = _candidate_alignments(self.truth[:60], [read], [0],
            _anchor_index([read], 5, 0, 60, 100), set(), diagnostics=diagnostics,
            low_quality_ry_rescue=enabled, **options)
        return result, diagnostics

    def test_single_low_quality_mismatch_rescued_both_strands(self):
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                read = self.partner(reverse=reverse)
                self.assertEqual(self.search(read, enabled=False)[0], [])
                result, diagnostics = self.search(read)
                self.assertEqual([(a.offset, a.reverse) for a in result], [(30, reverse)])
                self.assertEqual(diagnostics.ry_rescued_alignments, 1)

    def test_high_quality_strain_difference_and_two_errors_rejected(self):
        for read in (self.partner(quality=16), self.partner(quality=35),
                     self.partner(positions=(10, 15)), Read("missing", self.partner().sequence)):
            self.assertEqual(self.search(read)[0], [])

    def test_quality_at_wrong_raw_coordinate_does_not_rescue(self):
        read = self.partner(quality=35, reverse=True)
        qualities = list(read.qualities)
        qualities[15] = 10  # Correct raw coordinate is 49 - 15 = 34.
        self.assertEqual(self.search(Read(read.name, read.sequence, tuple(qualities)))[0], [])

    def test_other_gates_and_ambiguous_bases_preserved(self):
        read = self.partner()
        for options in (dict(minimum_overlap=31), dict(minimum_identity=1),
                        dict(minimum_anchor_matches=100)):
            self.assertEqual(self.search(read, **options)[0], [])
        sequence = list(read.sequence)
        sequence[10] = "N"
        self.assertEqual(self.search(Read("ambiguous", "".join(sequence), read.qualities))[0], [])

    def test_low_quality_biological_difference_is_indistinguishable(self):
        # The gate cannot identify the cause of a low-Q observed substitution.
        self.assertEqual(len(self.search(self.partner(quality=10))[0]), 1)

    def test_rescue_enables_supported_truth_extension_with_and_without_cache(self):
        for priority in (False, True):
            reads = [Read("center", self.truth[:60], (35,) * 60)]
            for index in range(5):
                partner = self.partner(positions=() if index < 3 else (10 + index,), quality=10)
                reads.append(Read(f"read-{index}", partner.sequence, partner.qualities))
            pool = ProgressiveSequencePool.from_reads(reads)
            pool = pool.replace_record(pool.records[0].corrected(Read("center", self.truth[:60])))
            options = dict(anchor_k=5, anchors_per_read=60, minimum_anchor_matches=1,
                           maximum_iterations=1, evidence_priority=priority)
            baseline, before = iterate_progressive_raw_extension(pool, **options)
            rescued, after = iterate_progressive_raw_extension(pool, low_quality_ry_rescue=True, **options)
            self.assertEqual(before.added_bases, 0)
            self.assertEqual(baseline.records[0].current.sequence, self.truth[:60])
            self.assertEqual(after.added_bases, 20)
            self.assertEqual(rescued.records[0].current.sequence, self.truth)
            self.assertEqual(dict(after.recruit_search)["ry_rescued_alignments"], 2)
