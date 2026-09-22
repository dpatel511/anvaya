import importlib.util
import random
import sys
import unittest
from pathlib import Path

from anvaya.reads import Read
from anvaya.sequences import reverse_complement
from anvaya.overlap_assembly import _MasterOverlapEdge
from anvaya.overlap_graph import _project_master_edges
from anvaya.overlap_progressive_links import _add_bidirected_edge
from anvaya.damage_string_graph import (
    _exact_containment_candidates,
    _physical_edge_signature,
)
from anvaya.raw_consensus import DamageProfile

path = Path(__file__).resolve().parents[1] / "experiments/string_layout.py"
spec = importlib.util.spec_from_file_location("string_layout", path)
layout = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = layout
spec.loader.exec_module(layout)


class StringLayoutTests(unittest.TestCase):
    def setUp(self):
        self.damage = DamageProfile((0.4,) * 5, (0.4,) * 5)

    def test_conflicting_offsets_remove_both_orientations(self):
        reads = [Read("a", "A" * 80), Read("b", "A" * 80)]
        edges, ambiguous = {}, set()
        first = _MasterOverlapEdge((0, False), (1, False), 30, 50)
        second = _MasterOverlapEdge((0, False), (1, False), 31, 49)
        self.assertTrue(_add_bidirected_edge(edges, reads, first, ambiguous))
        self.assertFalse(_add_bidirected_edge(edges, reads, second, ambiguous))
        self.assertEqual(edges, {})
        self.assertEqual(len(ambiguous), 1)
        self.assertFalse(_add_bidirected_edge(edges, reads, first, ambiguous))

    def test_physical_edge_signature_handles_unequal_reverse_offsets(self):
        reads = [Read("source", "A" * 100), Read("target", "A" * 80)]
        forward = _MasterOverlapEdge(
            (0, False), (1, False), 40, 60, ((50, "C"),),
        )
        edges = {}
        _add_bidirected_edge(edges, reads, forward)
        reverse = edges[((1, True), (0, True))]
        self.assertEqual(reverse.shift, 20)
        self.assertEqual(reverse.corrections, ((69, "G"),))
        self.assertEqual(
            _physical_edge_signature(forward, reads),
            _physical_edge_signature(reverse, reads),
        )
        distinct = _MasterOverlapEdge(
            (0, False), (1, False), 41, 59, ((50, "C"),),
        )
        self.assertNotEqual(
            _physical_edge_signature(forward, reads),
            _physical_edge_signature(distinct, reads),
        )

    def test_indexed_containment_matches_exhaustive_placements(self):
        rng = random.Random(198101)
        parents = [
            "".join(rng.choice("ACGT") for _ in range(length))
            for length in (45, 63, 80, 101)
        ]
        sequences = list(dict.fromkeys((
            *parents,
            parents[2][7:38],
            reverse_complement(parents[3][19:59]),
            "A" * 31,
            "A" * 70,
            "ACGT" * 7,
            "TT" + "ACGT" * 7 + "GG" + "ACGT" * 7,
        )))

        expected = {}
        for child_index, child in enumerate(sequences):
            placements = set()
            orientations = [(False, child)]
            reverse_child = reverse_complement(child)
            if reverse_child != child:
                orientations.append((True, reverse_child))
            for parent_index, parent in enumerate(sequences):
                if len(child) >= len(parent):
                    continue
                for reverse, oriented_child in orientations:
                    offset = parent.find(oriented_child)
                    while offset >= 0:
                        placements.add((parent_index, reverse, offset))
                        offset = parent.find(oriented_child, offset + 1)
            if placements:
                expected[child_index] = placements

        observed = _exact_containment_candidates(sequences)
        self.assertEqual(dict(observed), expected)
        self.assertEqual(list(observed), sorted(observed))

    def test_reciprocal_corrections_have_canonical_order(self):
        reads = [Read("a", "A" * 80), Read("b", "A" * 80)]
        edges, ambiguous = {}, set()
        edge = _MasterOverlapEdge(
            (0, False), (1, False), 30, 50, ((31, "C"), (70, "G"))
        )
        self.assertTrue(_add_bidirected_edge(edges, reads, edge, ambiguous))
        reverse = edges[((1, True), (0, True))]
        self.assertEqual(reverse.corrections, tuple(sorted(reverse.corrections)))
        self.assertTrue(_add_bidirected_edge(edges, reads, reverse, ambiguous))
        self.assertFalse(ambiguous)

    def test_equal_spelled_paths_merge_all_physical_provenance(self):
        left = Read("left", "ACGT" * 20)
        right = Read("right", left.sequence[40:] + "TGCA" * 10)
        reads = [left, right, left, right]
        edges = {}
        for source, target in ((0, 1), (2, 3)):
            _add_bidirected_edge(
                edges,
                reads,
                _MasterOverlapEdge((source, False), (target, False), 40, 40),
            )
        layouts = {}
        projection, diagnostics = _project_master_edges(
            reads, edges, name_prefix="path", path_layouts=layouts
        )
        self.assertEqual(len(projection), 1)
        self.assertEqual(diagnostics.merged_contigs, 4)
        self.assertEqual({node for node, _, _ in layouts["path_1"]}, {0, 1, 2, 3})

    def test_cycle_is_reported_and_not_spelled(self):
        reads = [Read(str(index), base * 80) for index, base in enumerate("ACG")]
        edges = {}
        for source, target in ((0, 1), (1, 2), (2, 0)):
            _add_bidirected_edge(
                edges,
                reads,
                _MasterOverlapEdge((source, False), (target, False), 40, 40),
            )
        projection, diagnostics = _project_master_edges(
            reads, edges, name_prefix="cycle"
        )
        self.assertEqual([read.sequence for read in projection], [r.sequence for r in reads])
        self.assertGreater(diagnostics.cyclic_components, 0)

    def test_duplicate_molecule_does_not_inflate_provenance(self):
        truth = "ACGTTGCA" * 20
        reads = [Read("a", truth[:80]), Read("b", truth[40:120])]
        pool, _ = layout.assemble(reads, [7, 7])
        self.assertEqual(pool.active_derived[0].contributing_molecules, frozenset({7}))

    def test_palindromic_raw_orientation_remains_ambiguous(self):
        palindrome = Read("palindrome", "ACGT" * 20, tuple(range(80)))
        pool, _ = layout.assemble([palindrome])
        placements = pool.active_derived[0].raw_placements
        self.assertEqual({placement.reverse for placement in placements}, {False, True})
        self.assertEqual({placement.read_index for placement in placements}, {0})

    def test_clean_chain_and_reverse_duplicate_provenance(self):
        rng = random.Random(94001)
        truth = "".join(rng.choice("ACGT") for _ in range(200))
        reads = [Read(f"r{i}", truth[start:start + 80]) for i, start in enumerate(range(0, 121, 20))]
        reads.append(Read("duplicate_reverse", reverse_complement(reads[2].sequence)))
        pool, diagnostics = layout.assemble(reads)
        self.assertEqual(len(pool.active_derived), 1)
        center = pool.active_derived[0]
        self.assertIn(center.current.sequence, (truth, reverse_complement(truth)))
        self.assertEqual({p.read_index for p in center.raw_placements}, set(range(len(reads))))
        for p in center.raw_placements:
            raw = reads[p.read_index].sequence
            oriented = reverse_complement(raw) if p.reverse else raw
            self.assertEqual(center.current.sequence[p.offset:p.offset + 80], oriented)
        self.assertGreater(diagnostics["transitive_edges_removed"], 0)

    def test_branch_does_not_spell_mosaic(self):
        rng = random.Random(94002)
        trunk = "".join(rng.choice("ACGT") for _ in range(80))
        branches = [trunk[40:] + base * 40 for base in "AC"]
        reads = [Read("trunk", trunk)] + [Read(f"b{i}", b) for i, b in enumerate(branches)]
        pool, diagnostics = layout.assemble(reads)
        self.assertGreater(diagnostics["ambiguous_ends"], 0)
        self.assertGreater(
            diagnostics["input_physical_edges"], diagnostics["reciprocal_edges"]
        )
        self.assertGreater(diagnostics["branch_rejected_physical_edges"], 0)
        self.assertTrue(all(len(r.current.sequence) == 80 for r in pool.active_derived))

    def test_terminal_damage_reconnects_chain(self):
        rng = random.Random(94003)
        truth = list("".join(rng.choice("ACGT") for _ in range(120)))
        truth[40] = "C"
        truth = "".join(truth)
        damaged_right = "T" + truth[41:120]
        reads = [
            Read("left", truth[:80], (35,) * 80),
            Read("left_support", truth[:80], (35,) * 80),
            Read("right", damaged_right, (35,) * 80),
        ]
        exact, _ = layout.assemble(reads)
        aware, diagnostics = layout.assemble(reads, damage_profile=self.damage)
        self.assertEqual(len(exact.active_derived), 2)
        self.assertEqual(len(aware.active_derived), 1)
        self.assertEqual(aware.active_derived[0].current.sequence, truth)
        self.assertGreater(diagnostics["candidate_classifications"]["damage_compatible"], 0)

    def test_reverse_complemented_terminal_damage_reconnects_chain(self):
        rng = random.Random(94006)
        truth = list("".join(rng.choice("ACGT") for _ in range(120)))
        truth[40] = "C"
        truth = "".join(truth)
        reads = [
            Read("left", reverse_complement(truth[:80]), (35,) * 80),
            Read("left_support", reverse_complement(truth[:80]), (35,) * 80),
            Read("right", reverse_complement("T" + truth[41:120]), (35,) * 80),
        ]
        pool, diagnostics = layout.assemble(reads, damage_profile=self.damage)
        self.assertEqual(len(pool.active_derived), 1)
        self.assertIn(pool.active_derived[0].current.sequence,
                      (truth, reverse_complement(truth)))
        self.assertGreater(diagnostics["candidate_classifications"]["damage_compatible"], 0)

    def test_zero_and_asymmetric_profiles_obey_raw_end_direction(self):
        rng = random.Random(94007)
        truth = list("".join(rng.choice("ACGT") for _ in range(120)))
        truth[40] = "C"
        truth = "".join(truth)
        forward = [
            Read("left", truth[:80], (35,) * 80),
            Read("support", truth[:80], (35,) * 80),
            Read("right", "T" + truth[41:120], (35,) * 80),
        ]
        zero = DamageProfile((0.0,) * 5, (0.0,) * 5)
        self.assertEqual(len(layout.assemble(forward, damage_profile=zero)[0].active_derived), 2)
        five_prime_only = DamageProfile((0.4,) * 5, (0.0,) * 5)
        self.assertEqual(len(layout.assemble(forward, damage_profile=five_prime_only)[0].active_derived), 1)
        reversed_reads = [
            Read(r.name, reverse_complement(r.sequence), tuple(reversed(r.qualities)))
            for r in forward
        ]
        self.assertEqual(
            len(layout.assemble(reversed_reads, damage_profile=five_prime_only)[0].active_derived),
            2,
        )

    def test_variable_length_chain_preserves_coordinates(self):
        rng = random.Random(95001)
        truth = "".join(rng.choice("ACGT") for _ in range(300))
        intervals = ((0, 100), (60, 180), (130, 260), (220, 300))
        reads = [Read(f"r{i}", truth[start:stop]) for i, (start, stop) in enumerate(intervals)]
        pool, _ = layout.assemble(reads)
        self.assertEqual(len(pool.active_derived), 1)
        contig = pool.active_derived[0]
        self.assertIn(contig.current.sequence, (truth, reverse_complement(truth)))
        for placement in contig.raw_placements:
            raw = reads[placement.read_index].sequence
            oriented = reverse_complement(raw) if placement.reverse else raw
            self.assertEqual(
                contig.current.sequence[placement.offset:placement.offset + len(raw)],
                oriented,
            )

    def test_unique_exact_containment_preserves_raw_molecule(self):
        rng = random.Random(95002)
        parent = "".join(rng.choice("ACGT") for _ in range(200))
        reads = [Read("parent", parent), Read("child", parent[40:71])]
        pool, diagnostics = layout.assemble(reads, [11, 12])
        self.assertEqual(len(pool.active_derived), 1)
        contig = pool.active_derived[0]
        self.assertEqual({p.read_index for p in contig.raw_placements}, {0, 1})
        self.assertEqual(contig.contributing_molecules, frozenset({11, 12}))
        self.assertEqual(diagnostics["contained_nodes"], 1)
        for placement in contig.raw_placements:
            raw = reads[placement.read_index].sequence
            oriented = reverse_complement(raw) if placement.reverse else raw
            self.assertEqual(
                contig.current.sequence[placement.offset:placement.offset + len(raw)],
                oriented,
            )

    def test_reverse_exact_containment_preserves_orientation(self):
        rng = random.Random(95005)
        parent = "".join(rng.choice("ACGT") for _ in range(200))
        child = reverse_complement(parent[40:71])
        reads = [Read("parent", parent), Read("child", child)]
        pool, diagnostics = layout.assemble(reads, [21, 22])
        self.assertEqual(len(pool.active_derived), 1)
        contig = pool.active_derived[0]
        self.assertEqual(contig.contributing_molecules, frozenset({21, 22}))
        self.assertEqual(diagnostics["contained_nodes"], 1)
        for placement in contig.raw_placements:
            raw = reads[placement.read_index].sequence
            oriented = reverse_complement(raw) if placement.reverse else raw
            self.assertEqual(
                contig.current.sequence[placement.offset:placement.offset + len(raw)],
                oriented,
            )

    def test_ambiguous_exact_containment_is_not_retired(self):
        rng = random.Random(95003)
        child = "".join(rng.choice("ACGT") for _ in range(40))
        left = "".join(rng.choice("ACGT") for _ in range(30))
        right = "".join(rng.choice("ACGT") for _ in range(30))
        reads = [
            Read("parent_a", left + child + "A" * 30),
            Read("parent_b", right + child + "C" * 30),
            Read("child", child),
        ]
        pool, diagnostics = layout.assemble(reads)
        self.assertEqual(len(pool.active_derived), 3)
        self.assertEqual(diagnostics["contained_nodes"], 0)
        self.assertEqual(diagnostics["ambiguous_containments"], 1)

    def test_unequal_parent_containment_is_not_retired(self):
        rng = random.Random(95101)
        child = "".join(rng.choice("ACGT") for _ in range(40))
        reads = [
            Read("short_parent", "".join(rng.choice("ACGT") for _ in range(30)) + child + "A" * 30),
            Read("long_parent", "".join(rng.choice("ACGT") for _ in range(50)) + child + "C" * 30),
            Read("child", child),
        ]
        pool, diagnostics = layout.assemble(reads)
        self.assertEqual(len(pool.active_derived), 3)
        self.assertEqual(diagnostics["contained_nodes"], 0)
        self.assertEqual(diagnostics["ambiguous_containments"], 1)

    def test_repeated_parent_placement_is_not_retired(self):
        rng = random.Random(95102)
        child = "".join(rng.choice("ACGT") for _ in range(31))
        parent = "".join(rng.choice("ACGT") for _ in range(20)) + child
        parent += "".join(rng.choice("ACGT") for _ in range(20)) + child
        parent += "".join(rng.choice("ACGT") for _ in range(20))
        pool, diagnostics = layout.assemble([
            Read("parent", parent), Read("child", child),
        ])
        self.assertEqual(len(pool.active_derived), 2)
        self.assertEqual(diagnostics["contained_nodes"], 0)
        self.assertEqual(diagnostics["ambiguous_containments"], 1)

    def test_ambiguous_repeat_containment_cannot_be_used_as_dovetail(self):
        motif = "ACGATTCG"
        parent = motif * 20
        child = motif * 5
        reads = [Read("parent", parent, (38,) * len(parent)),
                 Read("child", child, (38,) * len(child))]

        pool, diagnostics = layout.assemble(reads, [10, 11])

        self.assertEqual(len(pool.active_derived), 2)
        self.assertEqual(diagnostics["contained_nodes"], 0)
        self.assertEqual(diagnostics["ambiguous_containments"], 1)
        self.assertEqual(
            {record.contributing_molecules for record in pool.active_derived},
            {frozenset({10}), frozenset({11})},
        )
        for record in pool.active_derived:
            for placement in record.raw_placements:
                raw = reads[placement.read_index].sequence
                oriented = reverse_complement(raw) if placement.reverse else raw
                self.assertEqual(
                    record.current.sequence[placement.offset:placement.offset + len(raw)],
                    oriented,
                )

    def test_variable_length_layout_is_deterministic(self):
        rng = random.Random(95006)
        truth = "".join(rng.choice("ACGT") for _ in range(300))
        reads = [
            Read(f"r{i}", truth[start:stop])
            for i, (start, stop) in enumerate(((0, 100), (60, 180), (130, 260), (220, 300)))
        ]
        first, first_diagnostics = layout.assemble(reads)
        second, second_diagnostics = layout.assemble(reads)
        timings = {"stale": 1.0}
        timed, timed_diagnostics = layout.assemble(reads, stage_timings=timings)
        self.assertEqual(first, second)
        self.assertEqual(first, timed)
        self.assertEqual(first_diagnostics, second_diagnostics)
        self.assertEqual(first_diagnostics, timed_diagnostics)
        self.assertEqual(
            set(timings),
            {
                "indexing", "candidate_edges", "containment", "phase_audit",
                "projection", "provenance",
            },
        )
        self.assertTrue(all(seconds >= 0 for seconds in timings.values()))

    def test_variable_length_terminal_damage_reconnects(self):
        rng = random.Random(95004)
        truth = list("".join(rng.choice("ACGT") for _ in range(170)))
        truth[50] = "C"
        truth = "".join(truth)
        reads = [
            Read("left", truth[:90], (35,) * 90),
            Read("left_support", truth[:90], (35,) * 90),
            Read("right", "T" + truth[51:170], (35,) * 120),
        ]
        pool, _ = layout.assemble(reads, damage_profile=self.damage)
        self.assertEqual(len(pool.active_derived), 1)
        self.assertIn(pool.active_derived[0].current.sequence,
                      (truth, reverse_complement(truth)))

    def test_ordinary_terminal_mismatch_is_rejected(self):
        rng = random.Random(94008)
        truth = "".join(rng.choice("ACGT") for _ in range(120))
        changed = ("C" if truth[40] != "C" else "A") + truth[41:120]
        reads = [
            Read("left", truth[:80], (35,) * 80),
            Read("support", truth[:80], (35,) * 80),
            Read("right", changed, (35,) * 80),
        ]
        pool, _ = layout.assemble(reads, damage_profile=self.damage)
        self.assertEqual(len(pool.active_derived), 2)

    def test_internal_transition_is_rejected(self):
        rng = random.Random(94004)
        truth = list("".join(rng.choice("ACGT") for _ in range(120)))
        truth[50] = "C"
        truth = "".join(truth)
        changed = list(truth[40:120])
        changed[10] = "T"
        reads = [
            Read("left", truth[:80], (35,) * 80),
            Read("left_support", truth[:80], (35,) * 80),
            Read("right", "".join(changed), (35,) * 80),
        ]
        pool, diagnostics = layout.assemble(reads, damage_profile=self.damage)
        self.assertEqual(len(pool.active_derived), 2)
        classifications = diagnostics["candidate_classifications"]
        self.assertGreater(classifications.get("non_directional_or_internal", 0), 0)

    def test_damage_edge_requires_quality_and_distinct_molecules(self):
        rng = random.Random(94005)
        truth = list("".join(rng.choice("ACGT") for _ in range(120)))
        truth[40] = "C"
        truth = "".join(truth)
        damaged = "T" + truth[41:120]
        missing = [Read("left", truth[:80]), Read("left_support", truth[:80]),
                   Read("right", damaged)]
        pool, diagnostics = layout.assemble(missing, damage_profile=self.damage)
        self.assertEqual(len(pool.active_derived), 2)
        self.assertGreater(diagnostics["candidate_classifications"]["missing_quality"], 0)
        quality = [Read(r.name, r.sequence, (35,) * 80) for r in missing]
        pool, diagnostics = layout.assemble(quality, [7, 8, 7], self.damage)
        self.assertEqual(len(pool.active_derived), 2)
        self.assertGreater(diagnostics["candidate_classifications"]["insufficient_molecules"], 0)
        low_quality = list(quality)
        low_quality[2] = Read("right", damaged, (5,) + (35,) * 79)
        pool, _ = layout.assemble(low_quality, damage_profile=self.damage)
        self.assertEqual(len(pool.active_derived), 2)

    def test_incompatible_transitive_spelling_is_not_reduced(self):
        reads = [Read("a", "A" * 80), Read("b", "A" * 80), Read("c", "A" * 80)]
        edges = {}
        for edge in (
            _MasterOverlapEdge((0, False), (2, False), 40, 40, ((30, "A"),)),
            _MasterOverlapEdge((0, False), (1, False), 20, 60, ((30, "C"),)),
            _MasterOverlapEdge((1, False), (2, False), 20, 60),
        ):
            _add_bidirected_edge(edges, reads, edge)
        _, diagnostics = _project_master_edges(
            reads, edges, name_prefix="guard", verify_transitive_alleles=True
        )
        self.assertEqual(diagnostics.transitive_edges_removed, 0)

    def test_transitive_path_with_an_extra_correction_is_not_reduced(self):
        reads = [Read("a", "A" * 80), Read("b", "A" * 80), Read("c", "A" * 80)]
        edges = {}
        for edge in (
            _MasterOverlapEdge((0, False), (2, False), 40, 40),
            _MasterOverlapEdge((0, False), (1, False), 20, 60, ((50, "C"),)),
            _MasterOverlapEdge((1, False), (2, False), 20, 60),
        ):
            _add_bidirected_edge(edges, reads, edge)
        _, diagnostics = _project_master_edges(
            reads, edges, name_prefix="guard", verify_transitive_alleles=True
        )
        self.assertEqual(diagnostics.transitive_edges_removed, 0)
