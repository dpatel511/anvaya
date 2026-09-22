"""Synthetic checks of candidate-stage tracing and raw round outcomes."""

import random
import unittest

from anvaya.overlap_assembly import _Alignment, _CandidateDiagnostics, _candidate_alignments
from anvaya.overlap_index import _anchor_index
from anvaya.overlap_progressive import ProgressiveSequencePool, ProgressiveRawCluster, extend_progressive_raw_clusters
from anvaya.reads import Read


class ExtensionDiagnosticsTests(unittest.TestCase):
    def search(self, mutation=False, trace=None):
        rng = random.Random(81203)
        truth = "".join(rng.choice("ACGT") for _ in range(80))
        target = truth[:60]
        partner = list(truth[30:])
        if mutation:
            partner[15] = "A" if partner[15] in "CT" else "C"
        reads = [Read("partner", "".join(partner))]
        return _candidate_alignments(target, reads, [0], _anchor_index(reads, 5, 0, 60, 100), set(),
            anchor_k=5, anchors_per_read=60, maximum_anchor_occurrences=100, minimum_anchor_matches=1,
            minimum_overlap=30, minimum_identity=.9, minimum_ry_identity=.99,
            position_bits=(50).bit_length(), target_window=60, stage_trace=trace)

    def test_trace_is_output_neutral_and_tracks_true_offset(self):
        trace = {}
        self.assertEqual(self.search(trace=trace), self.search())
        self.assertEqual(trace[(0, False, 30)], "selected")

    def test_ry_rejection_distinguished_from_missing_seed(self):
        trace = {}
        self.search(mutation=True, trace=trace)
        self.assertEqual(trace[(0, False, 30)], "ry_identity_rejected")

    def test_candidate_index_diagnostics_separate_missing_and_capped_anchors(self):
        diagnostics = _CandidateDiagnostics()
        _candidate_alignments(
            "ACGTACGTACGT", [], [], {}, set(), anchor_k=5,
            anchors_per_read=4, maximum_anchor_occurrences=1,
            minimum_anchor_matches=1, minimum_overlap=5, minimum_identity=1,
            minimum_ry_identity=1, position_bits=1, target_window=12,
            diagnostics=diagnostics,
        )
        self.assertGreater(diagnostics.target_anchors, 0)
        self.assertEqual(diagnostics.missing_index_anchors, diagnostics.target_anchors)

        reads = [Read("a", "ACGTACGTACGT"), Read("b", "ACGTACGTACGT")]
        diagnostics = _CandidateDiagnostics()
        _candidate_alignments(
            reads[0].sequence, reads, [0, 1],
            _anchor_index(reads, 5, 0, 4, 100), set(), anchor_k=5,
            anchors_per_read=4, maximum_anchor_occurrences=1,
            minimum_anchor_matches=1, minimum_overlap=5, minimum_identity=1,
            minimum_ry_identity=1, position_bits=4, target_window=12,
            diagnostics=diagnostics,
        )
        self.assertGreater(diagnostics.occurrence_capped_anchors, 0)
        self.assertEqual(diagnostics.indexed_anchor_hits, 0)

    def test_contained_reads_are_not_extending_candidates(self):
        rng = random.Random(81204)
        sequence = "".join(rng.choice("ACGT") for _ in range(60))
        reads = [Read("seed", sequence), Read("contained", sequence[10:50])]
        pool = ProgressiveSequencePool.from_reads(reads)
        cluster = ProgressiveRawCluster(0, (0, 1), (_Alignment(1, reads[1].sequence, 10, 2, False),))
        updated, diagnostics = extend_progressive_raw_clusters(pool, (cluster,), maximum_rounds=1)
        self.assertEqual(dict(diagnostics.round_outcomes), {"no_dovetail": 1})
        self.assertEqual(updated.records[0].current.sequence, sequence)
        self.assertEqual(diagnostics.added_bases, 0)
