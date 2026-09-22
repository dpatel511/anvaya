"""Pre-layout phase evidence over complete synthetic overlap components."""

import unittest

from anvaya.overlap_assembly import _MasterOverlapEdge
from anvaya.phase_graph import (
    graph_component_coordinates,
    phase_graph_components,
)
from anvaya.raw_consensus import DamageProfile
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


def edge(source, target, shift):
    return _MasterOverlapEdge(source, target, shift, 1, ())


def edge_map(*items):
    return {(item.source, item.target): item for item in items}


class GraphCoordinateTests(unittest.TestCase):
    def test_consistent_branch_retains_parallel_node_coordinates(self):
        sequences = ["AAAA", "CCCC", "GGGG"]
        edges = edge_map(
            edge((0, False), (2, False), 12),
            edge((1, False), (2, False), 12),
        )

        component, = graph_component_coordinates(sequences, edges)

        self.assertTrue(component.consistent)
        self.assertEqual(component.physical_nodes, (0, 1, 2))
        self.assertEqual(
            tuple((item.node_index, item.reverse, item.offset)
                  for item in component.placements),
            ((0, False, 0), (1, False, 0), (2, False, 12)),
        )
        self.assertEqual(component.span, 16)

    def test_bidirected_mirror_component_is_reported_once(self):
        sequences = ["AAAA", "CCCC"]
        edges = edge_map(
            edge((0, False), (1, False), 3),
            edge((1, True), (0, True), 3),
        )

        components = graph_component_coordinates(sequences, edges)

        self.assertEqual(len(components), 1)
        self.assertEqual(components[0].physical_nodes, (0, 1))

    def test_contradictory_cycle_is_unresolved_without_arbitrary_coordinates(self):
        sequences = ["AAAA", "CCCC", "GGGG"]
        edges = edge_map(
            edge((0, False), (1, False), 10),
            edge((0, False), (2, False), 20),
            edge((1, False), (2, False), 15),
        )

        component, = graph_component_coordinates(sequences, edges)

        self.assertFalse(component.consistent)
        self.assertTrue(component.conflicts)
        proposals = {
            (conflict.assigned_offset, conflict.proposed_offset)
            for conflict in component.conflicts
        }
        self.assertTrue(any(left != right for left, right in proposals))

    def test_both_orientations_of_one_physical_node_are_unresolved(self):
        sequences = ["AAAA", "CCCC"]
        edges = edge_map(
            edge((0, False), (1, False), 4),
            edge((1, False), (0, True), 4),
        )

        component, = graph_component_coordinates(sequences, edges)

        self.assertFalse(component.consistent)
        self.assertEqual(component.orientation_conflicts, (0,))


class PreLayoutPhaseTests(unittest.TestCase):
    def test_linked_alternatives_survive_before_path_selection_on_both_strands(self):
        first = "AAAAAAAAAGAA"
        second = "AACAAAAAATAA"
        second_canonical = reverse_complement(second)
        reads = [
            Read("first-0", first, (38,) * 12),
            Read("first-1", first, (38,) * 12),
            Read("second-0", second, (38,) * 12),
            Read("second-1", second, (38,) * 12),
        ]
        sequences = [first, second_canonical, "GGGG"]
        groups = {
            first: ((0, False), (1, False)),
            second_canonical: ((2, True), (3, True)),
            "GGGG": (),
        }
        edges = edge_map(
            edge((0, False), (2, False), 12),
            edge((1, True), (2, False), 12),
        )

        result, = phase_graph_components(
            reads, [10, 11, 20, 21], sequences, groups, edges,
            DamageProfile((), ()),
        )

        self.assertTrue(result.coordinates.consistent)
        self.assertIsNotNone(result.phase)
        self.assertEqual(result.phase.phase.variable_positions_0based, (2, 9))
        self.assertEqual(result.phase.phase.blocks, ((2, 9),))
        self.assertEqual(
            result.phase.phase.links[0].allele_pair_molecules,
            ((0, 2, (10, 11)), (1, 3, (20, 21))),
        )

    def test_unknown_phase_conserved_bridge_does_not_link_separate_markers(self):
        reads = [
            Read(str(index), base, (38,))
            for index, base in enumerate("AACCGGTT")
        ] + [Read("bridge", "ACGT", (38,) * 4)]
        sequences = ["A", "C", "G", "T", "ACGT"]
        groups = {
            "A": ((0, False), (1, False)),
            "C": ((2, False), (3, False)),
            "G": ((4, False), (5, False)),
            "T": ((6, False), (7, False)),
            "ACGT": ((8, False),),
        }
        edges = edge_map(
            edge((0, False), (4, False), 3),
            edge((1, False), (4, False), 3),
            edge((4, False), (2, False), 5),
            edge((4, False), (3, False), 5),
        )

        result, = phase_graph_components(
            reads, list(range(9)), sequences, groups, edges,
            DamageProfile((), ()),
        )

        self.assertEqual(result.phase.phase.variable_positions_0based, (0, 8))
        self.assertEqual(result.phase.phase.blocks, ((0,), (8,)))
        self.assertEqual(result.phase.phase.unlinked_pairs_0based, ((0, 8),))
        self.assertEqual(result.phase.phase.links, ())

    def test_terminal_damage_does_not_create_a_robust_phase_marker(self):
        reads = [
            Read("c0", "C", (38,)), Read("c1", "C", (38,)),
            Read("t0", "T", (38,)), Read("t1", "T", (38,)),
        ]
        sequences = ["C", "T", "GG"]
        groups = {
            "C": ((0, False), (1, False)),
            "T": ((2, False), (3, False)),
            "GG": (),
        }
        edges = edge_map(
            edge((0, False), (2, False), 1),
            edge((1, False), (2, False), 1),
        )

        result, = phase_graph_components(
            reads, list(range(4)), sequences, groups, edges,
            DamageProfile((0.9,), ()),
        )

        self.assertEqual(result.phase.phase.variable_positions_0based, ())

    def test_single_ordinary_error_does_not_create_a_phase_marker(self):
        reads = [
            Read("a0", "A", (38,)), Read("a1", "A", (38,)),
            Read("a2", "A", (38,)), Read("error", "C", (38,)),
        ]
        sequences = ["A", "C", "GG"]
        groups = {
            "A": ((0, False), (1, False), (2, False)),
            "C": ((3, False),),
            "GG": (),
        }
        edges = edge_map(
            edge((0, False), (2, False), 1),
            edge((1, False), (2, False), 1),
        )

        result, = phase_graph_components(
            reads, list(range(4)), sequences, groups, edges,
            DamageProfile((), ()),
        )

        self.assertEqual(result.phase.phase.variable_positions_0based, ())

    def test_coordinate_conflict_blocks_phase_projection(self):
        reads = [Read("a", "A", (38,)), Read("c", "C", (38,))]
        sequences = ["A", "C", "GG"]
        groups = {"A": ((0, False),), "C": ((1, False),), "GG": ()}
        edges = edge_map(
            edge((0, False), (1, False), 3),
            edge((0, False), (2, False), 6),
            edge((1, False), (2, False), 4),
        )

        result, = phase_graph_components(
            reads, [0, 1], sequences, groups, edges, DamageProfile((), ()),
        )

        self.assertFalse(result.coordinates.consistent)
        self.assertIsNone(result.phase)

    def test_palindromic_node_has_ambiguous_strand_and_contributes_no_observation(self):
        reads = [Read("palindrome", "ACGT", (38,) * 4)]
        sequences = ["ACGT"]
        groups = {"ACGT": ((0, False),)}

        result, = phase_graph_components(
            reads, [7], sequences, groups, {}, DamageProfile((), ()),
        )

        self.assertEqual(result.phase.diagnostics.ambiguous_molecules, 1)
        self.assertEqual(result.phase.diagnostics.distinct_observed_molecules, 0)


if __name__ == "__main__":
    unittest.main()
