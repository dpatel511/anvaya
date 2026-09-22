"""Independent small-oracle checks for the Phase A exact string graph."""

import random
import unittest

from anvaya.exact_string_graph import (
    ExactEdge,
    _candidate_edges,
    _containments,
    _transitive_reduction,
    assemble_damage_aware_unitigs,
    assemble_exact_unitigs,
    evidence_pool,
)
from anvaya.paired_reads import merge_overlapping_pairs
from anvaya.raw_consensus import DamageProfile, project_raw_consensus
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


def exhaustive_edges(sequences, minimum_overlap):
    edges = {}
    for source_index, source in enumerate(sequences):
        for source_reverse in (False, True):
            oriented_source = reverse_complement(source) if source_reverse else source
            for target_index, target in enumerate(sequences):
                if target_index == source_index:
                    continue
                for target_reverse in (False, True):
                    oriented_target = reverse_complement(target) if target_reverse else target
                    for overlap in range(
                        minimum_overlap,
                        min(len(oriented_source), len(oriented_target)),
                    ):
                        if oriented_source[-overlap:] != oriented_target[:overlap]:
                            continue
                        edge = ExactEdge(
                            (source_index, source_reverse),
                            (target_index, target_reverse),
                            len(oriented_source) - overlap,
                            overlap,
                        )
                        key = (edge.source, edge.target)
                        if key not in edges or overlap > edges[key].overlap:
                            edges[key] = edge
    return edges


def exhaustive_containments(sequences):
    placements = {}
    for child_index, child in enumerate(sequences):
        found = set()
        oriented = [(False, child)]
        reverse = reverse_complement(child)
        if reverse != child:
            oriented.append((True, reverse))
        for parent_index, parent in enumerate(sequences):
            if len(child) >= len(parent):
                continue
            for is_reverse, oriented_child in oriented:
                start = 0
                while (offset := parent.find(oriented_child, start)) >= 0:
                    found.add((parent_index, is_reverse, offset))
                    start = offset + 1
        if found:
            placements[child_index] = found
    return placements


def exhaustive_reduction(edges):
    outgoing = {}
    for edge in edges.values():
        outgoing.setdefault(edge.source, []).append(edge)

    def alternative(direct):
        def visit(node, shift, used):
            if node == direct.target and shift == direct.shift:
                return True
            for edge in outgoing.get(node, ()):
                key = (edge.source, edge.target)
                total = shift + edge.shift
                if edge == direct or key in used or total > direct.shift:
                    continue
                if visit(edge.target, total, used | {key}):
                    return True
            return False
        return visit(direct.source, 0, set())

    return {
        key: edge for key, edge in edges.items()
        if not alternative(edge)
    }


class ExactStringGraphTests(unittest.TestCase):
    def test_indexed_dovetails_match_exhaustive_oracle(self):
        reads = (
            "AACCGGTTAACC",
            "TTAACCGGATTA",
            "CCGGATTATGCA",
            "GGTTAACCGGATTA",
            "TGCAAAAACCCC",
        )
        expected = exhaustive_edges(reads, 4)
        observed = _candidate_edges(reads, 4)
        self.assertEqual(observed, expected)
        for edge in observed.values():
            self.assertLess(edge.overlap, len(reads[edge.source[0]]))
            self.assertLess(edge.overlap, len(reads[edge.target[0]]))

    def test_indexed_geometry_matches_random_exhaustive_oracles(self):
        rng = random.Random(731)
        for _ in range(12):
            reference = "".join(rng.choice("ACGT") for _ in range(80))
            sequences = []
            for _ in range(18):
                start = rng.randrange(65)
                stop = min(80, start + rng.randrange(8, 18))
                sequence = reference[start:stop]
                sequences.append(min(sequence, reverse_complement(sequence)))
            sequences = tuple(dict.fromkeys(sequences))
            with self.subTest(reads=len(sequences)):
                indexed = _candidate_edges(sequences, 5)
                self.assertEqual(indexed, exhaustive_edges(sequences, 5))
                self.assertEqual(_containments(sequences), exhaustive_containments(sequences))
                self.assertEqual(
                    _transitive_reduction(indexed, sequences),
                    exhaustive_reduction(indexed),
                )

    def test_containment_is_separate_and_all_evidence_placements_survive(self):
        reads = [
            Read("parent_a", "AACCGGTTACGA"),
            Read("parent_b", "TTTAACCGGTTACGAGG"),
            Read("child", "CCGGTT"),
            Read("child_duplicate", "CCGGTT"),
            Read("reverse_child", reverse_complement("CCGGTT")),
        ]
        graph = assemble_exact_unitigs(reads, minimum_overlap=4)
        self.assertEqual(graph.contained_groups, 2)
        self.assertEqual(len(graph.topology_sequences), 1)
        by_raw = {index: [] for index in range(len(reads))}
        for placement in graph.evidence_placements:
            by_raw[placement.raw_index].append(placement)
        self.assertTrue(all(by_raw.values()))
        self.assertEqual(len(by_raw[2]), 2)
        self.assertEqual(len(by_raw[3]), 2)
        self.assertEqual(len(by_raw[4]), 2)
        self.assertEqual(graph.ambiguous_evidence_reads, 3)
        self.assertFalse(graph.candidate_edges)

    def test_multiple_maximal_placements_remain_ambiguous(self):
        graph = assemble_exact_unitigs([
            Read("left", "AAAACCCGTTGG"),
            Read("right", "TTTAACCCGTTAA"),
            Read("child", "AACCCGTT"),
        ], minimum_overlap=5)
        child = [
            placement for placement in graph.evidence_placements
            if placement.raw_index == 2
        ]
        self.assertEqual(len(child), 2)
        self.assertEqual(graph.ambiguous_evidence_reads, 1)
        for placement in graph.evidence_placements:
            raw = ["AAAACCCGTTGG", "TTTAACCCGTTAA", "AACCCGTT"][placement.raw_index]
            if placement.reverse:
                raw = reverse_complement(raw)
            topology = graph.topology_sequences[placement.topology_node]
            self.assertEqual(
                topology[placement.offset:placement.offset + len(raw)], raw
            )

    def test_transitive_edge_is_removed_and_unitig_is_deterministic(self):
        reads = [
            Read("a", "AAAACCCC"),
            Read("b", "AACCCCGG"),
            Read("c", "CCCCGGTT"),
        ]
        first = assemble_exact_unitigs(reads, minimum_overlap=4)
        permuted = assemble_exact_unitigs(list(reversed(reads)), minimum_overlap=4)
        self.assertLess(len(first.reduced_edges), len(first.candidate_edges))
        reduced = {(edge.source, edge.target): edge for edge in first.reduced_edges}
        for unitig in first.unitigs:
            offsets = [offset for _, _, offset in unitig.layout]
            self.assertEqual(offsets, sorted(offsets))
            for source, target in zip(unitig.layout, unitig.layout[1:]):
                edge = reduced[((source[0], source[1]), (target[0], target[1]))]
                self.assertEqual(target[2] - source[2], edge.shift)
        self.assertEqual(
            sorted(unitig.sequence for unitig in first.unitigs),
            sorted(unitig.sequence for unitig in permuted.unitigs),
        )
        self.assertEqual(
            sorted(unitig.sequence for unitig in first.unitigs),
            ["AAAACCCCGGTT"],
        )

    def test_branch_cycle_and_isolated_node_are_deterministic(self):
        branch_reads = [
            Read("stem", "AAAACGTC"),
            Read("first", "CGTCGGAA"),
            Read("second", "CGTCATAT"),
        ]
        branch = assemble_exact_unitigs(branch_reads, minimum_overlap=4)
        reversed_branch = assemble_exact_unitigs([
            Read(read.name, reverse_complement(read.sequence))
            for read in reversed(branch_reads)
        ], minimum_overlap=4)
        self.assertEqual(len(branch.unitigs), 2)
        self.assertEqual({len(unitig.layout) for unitig in branch.unitigs}, {2})
        self.assertEqual(
            sorted(unitig.sequence for unitig in branch.unitigs),
            sorted(unitig.sequence for unitig in reversed_branch.unitigs),
        )

        cycle_reads = [
            Read("cycle_a", "AAAACCCC"),
            Read("cycle_b", "CCCCAAAA"),
            Read("isolated", "ACGTACGT"),
        ]
        first = assemble_exact_unitigs(cycle_reads, minimum_overlap=4)
        reversed_input = assemble_exact_unitigs([
            Read(read.name, reverse_complement(read.sequence))
            for read in reversed(cycle_reads)
        ], minimum_overlap=4)
        self.assertEqual(
            sorted(unitig.sequence for unitig in first.unitigs),
            sorted(unitig.sequence for unitig in reversed_input.unitigs),
        )
        self.assertIn("ACGTACGT", {unitig.sequence for unitig in first.unitigs})
        self.assertTrue(any(len(unitig.layout) == 2 for unitig in first.unitigs))

    def test_prefix_and_suffix_containments_never_become_edges(self):
        graph = assemble_exact_unitigs([
            Read("parent", "AAAACCCCGGGG"),
            Read("prefix", "AAAACC"),
            Read("suffix", "CCGGGG"),
            Read("internal", "AACCCC"),
        ], minimum_overlap=4)
        self.assertEqual(graph.contained_groups, 3)
        self.assertEqual(graph.topology_sequences, ("AAAACCCCGGGG",))
        self.assertFalse(graph.candidate_edges)
        self.assertEqual(
            {placement.raw_index for placement in graph.evidence_placements},
            set(range(4)),
        )

    def test_clean_merged_regression_recovers_reference(self):
        rng = random.Random(199001)
        reference = "".join(rng.choice("ACGT") for _ in range(2400))
        placement = random.Random(199002)
        left, right, fragments = [], [], []
        for index in range(round(20 * 2400 / 135)):
            length = placement.randint(90, 180)
            start = placement.randrange(2400 - length + 1)
            fragment = reference[start:start + length]
            fragments.append(fragment)
            left.append(Read(f"pair_{index}/1", fragment[:75], (35,) * 75))
            right.append(Read(
                f"pair_{index}/2", reverse_complement(fragment[-75:]), (35,) * 75
            ))
        merged, _, diagnostics = merge_overlapping_pairs(left, right)
        self.assertEqual(diagnostics.merged_pairs, 138)
        for read in merged:
            if "|merged" in read.name:
                index = int(read.name.split("_")[1].split("/")[0])
                self.assertEqual(read.sequence, fragments[index])

        graph = assemble_exact_unitigs(merged, minimum_overlap=30)
        # No fragment covers reference coordinate 0; the attainable exact layout
        # is the 0-based half-open interval [1, 2400).
        self.assertIn(reference[1:], {unitig.sequence for unitig in graph.unitigs})
        self.assertEqual(max(map(len, (unitig.sequence for unitig in graph.unitigs))), 2399)

    def test_damage_layer_reconnects_chain_and_projects_each_raw_molecule_once(self):
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
        profile = DamageProfile((.4,) * 5, (.4,) * 5)

        exact = assemble_exact_unitigs(reads)
        aware, diagnostics = assemble_damage_aware_unitigs(reads, profile)
        pool, placement_diagnostics = evidence_pool(reads, aware)
        polished, _ = project_raw_consensus(pool, profile)

        self.assertEqual(len(exact.unitigs), 2)
        expected = min(truth, reverse_complement(truth))
        self.assertEqual([unitig.sequence for unitig in aware.unitigs], [expected])
        self.assertEqual([read.sequence for read in polished], [expected])
        self.assertGreater(
            diagnostics["candidate_classifications"].get("damage_compatible", 0),
            0,
        )
        self.assertEqual(placement_diagnostics["placed_raw_reads"], 3)
        placements = pool.active_derived[0].raw_placements
        self.assertEqual({placement.read_index for placement in placements}, {0, 1, 2})
        self.assertEqual(len(placements), 3)

        transformed_inputs = [
            list(reversed(reads)),
            [
                Read(
                    read.name,
                    reverse_complement(read.sequence),
                    tuple(reversed(read.qualities)),
                )
                for read in reversed(reads)
            ],
        ]
        for transformed in transformed_inputs:
            with self.subTest(transformed=[read.name for read in transformed]):
                transformed_graph, _ = assemble_damage_aware_unitigs(
                    transformed, profile,
                )
                transformed_pool, transformed_diagnostics = evidence_pool(
                    transformed, transformed_graph,
                )
                transformed_polished, _ = project_raw_consensus(
                    transformed_pool, profile,
                )
                self.assertEqual(
                    [unitig.sequence for unitig in transformed_graph.unitigs],
                    [expected],
                )
                self.assertEqual(
                    [read.sequence for read in transformed_polished], [expected]
                )
                self.assertEqual(transformed_diagnostics["placed_raw_reads"], 3)

    def test_global_branch_ambiguity_is_not_allocated_or_duplicated(self):
        reads = [
            Read("stem", "AAAACGTC"),
            Read("first", "CGTCGGAA"),
            Read("second", "CGTCATAT"),
        ]
        graph = assemble_exact_unitigs(reads, minimum_overlap=4)
        pool, diagnostics = evidence_pool(reads, graph)
        allocated = [
            placement.read_index
            for record in pool.active_derived
            for placement in record.raw_placements
        ]
        self.assertEqual(diagnostics["ambiguous_raw_indices"], [0])
        self.assertNotIn(0, allocated)
        self.assertEqual(sorted(allocated), [1, 2])

    def test_unique_contained_evidence_is_composed_onto_unitig(self):
        reads = [
            Read("parent", "AAAACCCCGGGG"),
            Read("child", "AACCCC"),
        ]
        graph = assemble_exact_unitigs(reads, minimum_overlap=4)
        pool, diagnostics = evidence_pool(reads, graph)
        self.assertEqual(diagnostics["placed_raw_reads"], 2)
        placements = {p.read_index: p for p in pool.active_derived[0].raw_placements}
        self.assertEqual((placements[0].offset, placements[0].reverse), (0, False))
        self.assertEqual((placements[1].offset, placements[1].reverse), (2, False))

    def test_transitive_reduction_requires_equal_corrected_spelling(self):
        sequences = ("A" * 80, "A" * 80, "A" * 80)
        direct = ExactEdge((0, False), (2, False), 40, 40, ((30, "A"),))
        first = ExactEdge((0, False), (1, False), 20, 60, ((30, "C"),))
        second = ExactEdge((1, False), (2, False), 20, 60)
        edges = {
            (edge.source, edge.target): edge
            for edge in (direct, first, second)
        }
        reduced = _transitive_reduction(edges, sequences)
        self.assertIn((direct.source, direct.target), reduced)


if __name__ == "__main__":
    unittest.main()
