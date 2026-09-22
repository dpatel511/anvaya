import importlib.util
import random
import sys
import unittest
from pathlib import Path

from anvaya.damage_string_graph import assemble
from anvaya.overlap_assembly import _MasterOverlapEdge
from anvaya.overlap_graph import _project_master_edges
from anvaya.overlap_progressive_links import _add_bidirected_edge
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


path = Path(__file__).resolve().parents[1] / "experiments/layout_loss_trace.py"
spec = importlib.util.spec_from_file_location("layout_loss_trace", path)
trace = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = trace
spec.loader.exec_module(trace)


class LayoutLossTraceTests(unittest.TestCase):
    @staticmethod
    def _redundant_path_fixture():
        rng = random.Random(20260913)
        reference = "".join(rng.choice("ACGT") for _ in range(180))
        reads = [
            Read("a", reference[0:100]), Read("b", reference[20:110]),
            Read("c", reference[40:140]), Read("d", reference[60:160]),
        ]
        edges = {}
        for source, target, shift in ((0, 1, 20), (1, 2, 20), (2, 3, 20), (0, 3, 60)):
            overlap = len(reads[source].sequence) - shift
            _add_bidirected_edge(
                edges, reads,
                _MasterOverlapEdge((source, False), (target, False), shift, overlap),
            )
        return reference, reads, edges

    def test_three_edge_equivalent_witness_restores_clean_path(self):
        reference, reads, edges = self._redundant_path_fixture()
        direct = edges[((0, False), (3, False))]
        witnesses, _, limited = trace._bounded_equivalent_witnesses(reads, edges, direct)
        self.assertFalse(limited)
        self.assertEqual([len(path) for path in witnesses], [3])
        reverse = edges[((3, True), (0, True))]
        reverse_witnesses, _, limited = trace._bounded_equivalent_witnesses(
            reads, edges, reverse
        )
        self.assertFalse(limited)
        self.assertEqual([len(path) for path in reverse_witnesses], [3])

        physical = trace._physical(direct.source, direct.target)
        reduced = {
            key: edge for key, edge in edges.items()
            if trace._physical(edge.source, edge.target) != physical
        }
        layouts = {}
        projected, diagnostics = _project_master_edges(
            reads, reduced, name_prefix="regression", path_layouts=layouts,
            verify_transitive_alleles=True,
        )
        self.assertEqual([read.sequence for read in projected], [reference[:160]])
        self.assertEqual(diagnostics.merged_contigs, 4)
        self.assertEqual(
            sorted(layouts["regression_1"], key=lambda item: item[2]),
            [(0, False, 0), (1, False, 20), (2, False, 40), (3, False, 60)],
        )

    def test_equivalent_witness_rejects_correction_disagreement(self):
        _, reads, edges = self._redundant_path_fixture()
        direct = edges[((0, False), (3, False))]
        different = next(base for base in "ACGT" if base != reads[0].sequence[10])
        corrected = _MasterOverlapEdge(
            direct.source, direct.target, direct.shift, direct.overlap, ((10, different),)
        )
        witnesses, _, _ = trace._bounded_equivalent_witnesses(reads, edges, corrected)
        self.assertEqual(witnesses, [])

    def test_equivalent_witness_abstains_on_repeated_paths_and_cycles(self):
        reference, reads, edges = self._redundant_path_fixture()
        reads.append(Read("b_duplicate", reads[1].sequence))
        for source, target, shift in ((0, 4, 20), (4, 2, 20), (2, 0, 1)):
            overlap = len(reads[source].sequence) - shift
            _add_bidirected_edge(
                edges, reads,
                _MasterOverlapEdge((source, False), (target, False), shift, overlap),
            )
        direct = edges[((0, False), (3, False))]
        witnesses, _, limited = trace._bounded_equivalent_witnesses(reads, edges, direct)
        self.assertFalse(limited)
        self.assertEqual(len(witnesses), 2)
        self.assertEqual({trace._spell_edges(reads, path) for path in witnesses}, {reference[:160]})

    def test_maximal_nonbranching_paths_preserve_both_real_branches(self):
        rng = random.Random(20260914)
        prefix = "".join(rng.choice("ACGT") for _ in range(100))
        left = Read("left", prefix)
        first = Read("first", prefix[40:] + "A" * 40)
        second = Read("second", prefix[40:] + "C" * 40)
        reads = [left, first, second]
        edges = {}
        for target in (1, 2):
            _add_bidirected_edge(
                edges, reads,
                _MasterOverlapEdge((0, False), (target, False), 40, 60),
            )
        result = trace._maximal_nonbranching_variant(
            reads, edges, {}, [read.sequence for read in reads], set()
        )
        self.assertEqual(result["path_sequences"], 2)
        self.assertEqual(result["contigs"], 2)
        self.assertEqual(result["longest"], 140)
        self.assertGreater(result["nodes_reused_across_paths"], 0)

    def test_expected_edges_receive_terminal_dispositions(self):
        case = trace.trace_case(95031)
        self.assertEqual(
            sum(case["expected_edge_dispositions"].values()),
            case["expected_physical_dovetails"],
        )
        self.assertNotIn("unclassified", case["expected_edge_dispositions"])
        self.assertEqual(case["selected_edges_outside_truth_dovetails"], 0)

if __name__ == "__main__":
    unittest.main()
