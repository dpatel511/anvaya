"""Verify synthetic benchmark geometry and metric accounting independently."""

import importlib.util
import unittest
from pathlib import Path

from anvaya.sequences import reverse_complement

path = Path(__file__).resolve().parents[1] / "experiments/linked_guard_validation.py"
spec = importlib.util.spec_from_file_location("linked_guard_validation", path)
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


class LinkedGuardValidationTests(unittest.TestCase):
    def test_damage_free_reads_match_source_on_both_strands(self):
        pool, _, truth, covered, competitors = benchmark.fixture(83001, 24, 0, 0, "matched", 6, 93)
        self.assertEqual(competitors, 0)
        placements = pool.records[0].raw_placements
        self.assertEqual({p.reverse for p in placements}, {True, False})
        union = set()
        for placement in placements:
            raw = pool.records[placement.read_index].raw.sequence
            oriented = reverse_complement(raw) if placement.reverse else raw
            self.assertEqual(oriented, truth[placement.offset:placement.offset + len(raw)])
            union.update(range(placement.offset, placement.offset + len(raw)))
        self.assertEqual(union, covered)

    def test_competitor_coordinates_are_known(self):
        pool, _, truth, _, competitors = benchmark.fixture(83002, 24, 1, 0, "matched", 35, 93)
        self.assertEqual(competitors, 24)
        other = list(truth)
        other[40], other[75] = "C", "T"
        other = "".join(other)
        for placement in pool.records[0].raw_placements:
            raw = pool.records[placement.read_index].raw.sequence
            oriented = reverse_complement(raw) if placement.reverse else raw
            self.assertEqual(oriented, other[placement.offset:placement.offset + len(raw)])

    def test_reproducibility_and_profile_change_preserves_observations(self):
        args = (83003, 12, .5, .2, "matched", 6, 30)
        first = benchmark.fixture(*args)
        self.assertEqual(first, benchmark.fixture(*args))
        over = benchmark.fixture(*args[:4], "overestimated", *args[5:])
        self.assertEqual(first[0], over[0])
        self.assertNotEqual(first[1], over[1])

    def test_metric_partition_and_error_delta(self):
        for seed in range(83004, 83014):
            pool, profile, truth, covered, _ = benchmark.fixture(seed, 6, .5, .4, "overestimated", 6, 20)
            metrics = benchmark.evaluate(pool, profile, truth, covered, [40])
            self.assertEqual(metrics["blocked_changes"], metrics["prevented_errors"] +
                             metrics["blocked_beneficial_corrections"] + metrics["both_wrong"])
            self.assertEqual(metrics["unguarded_errors"] - metrics["guarded_errors"],
                             metrics["prevented_errors"] - metrics["blocked_beneficial_corrections"])
            self.assertLessEqual(metrics["guarded_source_alleles"], metrics["source_variant_sites"])
