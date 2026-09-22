"""Synthetic guard controls, including the observed two-group pattern."""

import csv
import io
import unittest
from dataclasses import replace

from anvaya.linked_allele_guard import blocking_neighbour, variable_neighbours
from anvaya.overlap_progressive import ProgressiveSequencePool
from anvaya.raw_consensus import DamageProfile, RawPlacement, observation_likelihoods, project_raw_consensus
from anvaya.reads import Read


def obs(base, terminal=False):
    return ("ACGT".index(base), 38, observation_likelihoods(
        base, 38, 0 if terminal else 10, 30, False, DamageProfile((.4,), ())))


class LinkedAlleleGuardTests(unittest.TestCase):
    def columns(self):
        focal = {0: obs("T"), 1: obs("T", True), **{i: obs("C") for i in range(2, 6)}}
        neighbour = {i: obs("C" if i < 2 else "T") for i in range(6)}
        return focal, neighbour

    def test_two_groups_block_with_neighbour_at_zero(self):
        focal, neighbour = self.columns()
        self.assertEqual(blocking_neighbour(1, focal, 3, 1, variable_neighbours([neighbour, focal])), 0)

    def test_focal_site_cannot_supply_its_own_witness(self):
        focal, neighbour = self.columns()
        self.assertIsNone(blocking_neighbour(0, focal, 3, 1, variable_neighbours([neighbour])))

    def test_discordant_combination_prevents_blocking(self):
        focal, neighbour = self.columns()
        neighbour[2] = obs("C")
        self.assertIsNone(blocking_neighbour(0, focal, 3, 1, variable_neighbours([focal, neighbour])))

    def test_missing_neighbour_observation_prevents_cherry_picking(self):
        focal, neighbour = self.columns()
        del neighbour[5]
        self.assertIsNone(blocking_neighbour(0, focal, 3, 1, variable_neighbours([focal, neighbour])))

    def test_singleton_error_group_is_not_protected(self):
        focal, neighbour = self.columns()
        del focal[1]
        self.assertIsNone(blocking_neighbour(0, focal, 3, 1, variable_neighbours([focal, neighbour])))

    def test_damage_only_focal_support_does_not_block(self):
        focal, neighbour = self.columns()
        focal[0] = obs("T", True)
        self.assertIsNone(blocking_neighbour(0, focal, 3, 1, variable_neighbours([focal, neighbour])))

    def test_unrelated_molecules_do_not_link(self):
        focal, neighbour = self.columns()
        neighbour = {i+10: o for i, o in neighbour.items()}
        self.assertIsNone(blocking_neighbour(0, focal, 3, 1, variable_neighbours([focal, neighbour])))

    def test_excluded_duplicate_cannot_complete_a_group(self):
        focal, neighbour = self.columns()
        focal[1] = None
        self.assertIsNone(blocking_neighbour(0, focal, 3, 1, variable_neighbours([focal, neighbour])))

    def pool(self, damage_only=False):
        # Target positions 0 and 6. Raw intervals and offsets are 0-based/half-open.
        reads, placements = [], []
        for i in range(6):
            p = 1 if i == 1 or (damage_only and i == 0) else 20
            focal, neighbour = ("T", "C") if i < 2 else ("C", "T")
            sequence = "A" * p + focal + "A" * 5 + neighbour + "A" * 10
            reads.append(Read(str(i), sequence, (38,) * len(sequence)))
            placements.append(RawPlacement(i, -p, False, p, p+7))
        pool = ProgressiveSequencePool.from_reads(reads)
        center = replace(pool.records[0].corrected(Read("contig", "TAAAAAC")), raw_placements=tuple(placements))
        return pool.replace_record(center)

    def test_guard_changes_only_proposed_substitution_and_reports_witness(self):
        pool = self.pool()
        profile = DamageProfile((.4,) * 5, ())
        baseline, diagnostics = project_raw_consensus(pool, profile)
        report = io.StringIO()
        guarded, guarded_diagnostics = project_raw_consensus(pool, profile, linked_allele_guard=True, report=report)
        self.assertEqual(baseline[0].sequence, "CAAAAAC")
        self.assertEqual(guarded[0].sequence, "TAAAAAC")
        self.assertEqual(guarded_diagnostics.linked_allele_rejections, 1)
        self.assertEqual(guarded_diagnostics.changed_bases, diagnostics.changed_bases - 1)
        rows = list(csv.DictReader(io.StringIO(report.getvalue()), delimiter="\t"))
        focal = next(row for row in rows if row["position_0based"] == "0")
        self.assertEqual(focal["reason"], "linked_allele_conflict")
        self.assertEqual(focal["guard_neighbour_0based"], "6")
        self.assertEqual(pool.records[0].current.sequence, "TAAAAAC")

    def test_genuine_terminal_damage_correction_is_retained(self):
        pool = self.pool(damage_only=True)
        profile = DamageProfile((.4,) * 5, ())
        baseline, _ = project_raw_consensus(pool, profile)
        guarded, diagnostics = project_raw_consensus(pool, profile, linked_allele_guard=True)
        self.assertEqual(guarded, baseline)
        self.assertEqual(guarded[0].sequence[0], "C")
        self.assertEqual(diagnostics.linked_allele_rejections, 0)
