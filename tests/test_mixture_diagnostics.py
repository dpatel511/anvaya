"""Small synthetic fixtures; no actual FASTQ datasets or fitted thresholds."""

import csv
import io
import math
import unittest
from dataclasses import replace
from itertools import combinations

from anvaya.mixture_diagnostics import (
    FRACTIONS, diagnostic_writers, mixture_evidence, write_diagnostics,
)
from anvaya.raw_consensus import DamageProfile, RawPlacement, observation_likelihoods, project_raw_consensus
from anvaya.overlap_progressive import ProgressiveSequencePool
from anvaya.reads import Read


def observation(base, q=30, position=10, length=31, reverse=False, rate=0):
    vector = observation_likelihoods(base, q, position, length, reverse, DamageProfile((rate,), (rate,)))
    allele = "ACGT".index(base)
    return (3 - allele if reverse else allele, q, vector)


class MixtureDiagnosticsTests(unittest.TestCase):
    def test_discrete_evidence_matches_direct_enumeration(self):
        observations = [observation("C"), observation("T"), observation("T", position=0, rate=0.4)]
        vectors = [o[2] for o in observations]
        single = sum(math.prod(v[a] for v in vectors) for a in range(4)) / 4
        mixture = sum(math.prod(f * v[a] + (1-f) * v[b] for v in vectors)
                      for a, b in combinations(range(4), 2) for f in FRACTIONS) / (6 * len(FRACTIONS))
        self.assertAlmostEqual(mixture_evidence(observations)["log_bf_mixture_vs_single"], math.log(mixture / single))

    def test_no_observations_or_one_observation_does_not_favour_mixture(self):
        for observations in ([], [observation("T", position=0, rate=0.4)]):
            self.assertAlmostEqual(mixture_evidence(observations)["log_bf_mixture_vs_single"], 0)

    def test_internal_mixture_differs_from_single_allele(self):
        mixed = [observation("T")] * 3 + [observation("C")] * 3
        pure = [observation("C")] * 6
        self.assertGreater(mixture_evidence(mixed)["log_bf_mixture_vs_single"], 0)
        self.assertEqual(mixture_evidence(mixed)["best_pair"], "CT")
        self.assertLess(mixture_evidence(pure)["log_bf_mixture_vs_single"], 0)

    def test_damage_explanation_reduces_apparent_mixture_evidence(self):
        def observations(rate):
            return [observation("T", position=0, rate=rate)] * 3 + [observation("C")] * 2
        self.assertLess(mixture_evidence(observations(.4))["log_bf_mixture_vs_single"],
                        mixture_evidence(observations(0))["log_bf_mixture_vs_single"])

    def test_reverse_complement_and_order_invariance(self):
        forward = [observation("C"), observation("T", position=0, rate=.4), observation("C")]
        reverse = [(3-a, q, tuple(reversed(v))) for a, q, v in forward]
        self.assertAlmostEqual(mixture_evidence(forward)["log_bf_mixture_vs_single"],
                               mixture_evidence(reverse)["log_bf_mixture_vs_single"])
        self.assertEqual(mixture_evidence(forward), mixture_evidence(list(reversed(forward))))

    def test_high_depth_is_finite(self):
        value = mixture_evidence([observation("C")] * 5000 + [observation("T")] * 5000)
        self.assertTrue(math.isfinite(value["log_bf_mixture_vs_single"]))

    def reports(self, columns):
        mixture, linkage = io.StringIO(), io.StringIO()
        writers = diagnostic_writers(mixture, linkage)
        write_diagnostics(columns, "contig", "T" * len(columns), "C" * len(columns), *writers)
        return (list(csv.DictReader(io.StringIO(mixture.getvalue()), delimiter="\t")),
                list(csv.DictReader(io.StringIO(linkage.getvalue()), delimiter="\t")))

    def test_low_support_is_reported_not_promoted(self):
        rows, links = self.reports([{0: observation("T"), 1: observation("C")}])
        self.assertEqual(rows[0]["support_status"], "insufficient_support")
        self.assertEqual(links, [])

    def test_shared_molecule_linkage_and_focal_exclusion(self):
        first = {i: observation("T" if i < 2 else "C") for i in range(4)}
        second = {i: observation("A" if i < 2 else "G") for i in range(4)}
        rows, links = self.reports([first, second])
        self.assertEqual([r["linked_neighbour_sites"] for r in rows], ["1", "1"])
        focal = [r for r in links if r["focal_position_0based"] == "0"]
        self.assertEqual({(r["focal_observed_allele"], r["neighbour_robust_allele"], r["molecules"])
                          for r in focal}, {("T", "A", "2"), ("C", "G", "2")})
        self.assertTrue(all(r["focal_position_0based"] != r["neighbour_position_0based"] for r in links))

    def test_unlinked_molecules_do_not_manufacture_linkage(self):
        first = {i: observation("T" if i < 2 else "C") for i in range(4)}
        second = {i+4: observation("A" if i < 2 else "G") for i in range(4)}
        _, links = self.reports([first, second])
        self.assertEqual(links, [])

    def test_terminal_ambiguous_neighbour_not_used_to_define_groups(self):
        first = {i: observation("T" if i < 2 else "C") for i in range(4)}
        second = {i: observation("T" if i < 2 else "C", position=0, rate=.4) for i in range(4)}
        _, links = self.reports([first, second])
        self.assertFalse(any(r["focal_position_0based"] == "0" for r in links))

    def test_excluded_molecule_not_counted(self):
        rows, _ = self.reports([{0: observation("T"), 1: observation("C"), 2: None}])
        self.assertEqual(rows[0]["molecules"], "2")

    def test_optional_reports_leave_consensus_and_diagnostics_unchanged(self):
        reads = [Read(str(i), sequence, (30,) * 12) for i, sequence in enumerate(
            ["TAAAAAAAAAAA"] * 3 + ["AAAAACAAAAAA"] * 2)]
        pool = ProgressiveSequencePool.from_reads(reads)
        center = replace(pool.records[0].corrected(Read("contig", "T")), raw_placements=tuple(
            RawPlacement(i, 0 if i < 3 else -5, False, 0 if i < 3 else 5, 1 if i < 3 else 6)
            for i in range(5)))
        pool = pool.replace_record(center)
        profile = DamageProfile((.4,), ())
        expected, expected_diagnostics = project_raw_consensus(pool, profile)
        mixture, linkage = io.StringIO(), io.StringIO()
        actual, diagnostics = project_raw_consensus(pool, profile, mixture_report=mixture, linkage_report=linkage)
        self.assertEqual(actual, expected)
        self.assertEqual(diagnostics, expected_diagnostics)
        rows = list(csv.DictReader(io.StringIO(mixture.getvalue()), delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["before"], "T")
        self.assertEqual(rows[0]["after"], "C")
