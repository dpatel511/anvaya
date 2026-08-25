import unittest
from collections import Counter

from anvaya.graph import (
    branching_nodes,
    build_dbg,
    in_degree,
    out_degree,
    sink_nodes,
    source_nodes,
)


class BuildDbgTests(unittest.TestCase):
    def test_builds_expected_path(self) -> None:
        graph = build_dbg(["ATGCA"], 3)

        self.assertEqual(
            graph,
            {
                "AT": Counter({"TG": 1}),
                "TG": Counter({"GC": 1}),
                "GC": Counter({"CA": 1}),
                "CA": Counter(),
            },
        )

    def test_counts_repeated_edges(self) -> None:
        graph = build_dbg(["ATGC", "ATGC"], 3)

        self.assertEqual(graph["AT"]["TG"], 2)
        self.assertEqual(graph["TG"]["GC"], 2)

    def test_represents_a_branch(self) -> None:
        graph = build_dbg(["ATGCA", "ATGGA"], 3)

        self.assertEqual(graph["TG"], Counter({"GC": 1, "GG": 1}))
        self.assertIn("CA", graph)
        self.assertIn("GA", graph)

    def test_rejects_empty_reads(self) -> None:
        with self.assertRaises(ValueError):
            build_dbg([], 3)

    def test_rejects_k_less_than_two(self) -> None:
        with self.assertRaises(ValueError):
            build_dbg(["ATGC"], 1)

    def test_rejects_reads_shorter_than_k(self) -> None:
        with self.assertRaises(ValueError):
            build_dbg(["AT"], 3)


class GraphDegreeTests(unittest.TestCase):
    def test_degrees_for_linear_path(self) -> None:
        graph = build_dbg(["ATGCA"], 3)

        self.assertEqual((in_degree(graph, "AT"), out_degree(graph, "AT")), (0, 1))
        self.assertEqual((in_degree(graph, "TG"), out_degree(graph, "TG")), (1, 1))
        self.assertEqual((in_degree(graph, "CA"), out_degree(graph, "CA")), (1, 0))

    def test_multiplicity_does_not_change_degree(self) -> None:
        graph = build_dbg(["ATGC", "ATGC"], 3)

        self.assertEqual(out_degree(graph, "AT"), 1)
        self.assertEqual(in_degree(graph, "TG"), 1)

    def test_identifies_sources_sinks_and_branch(self) -> None:
        graph = build_dbg(["ATGCA", "ATGGA"], 3)

        self.assertEqual(source_nodes(graph), {"AT"})
        self.assertEqual(sink_nodes(graph), {"CA", "GA"})
        self.assertEqual(branching_nodes(graph), {"TG"})

    def test_identifies_merging_node_as_branching(self) -> None:
        graph = build_dbg(["ATGCA", "CTGCA"], 3)

        self.assertEqual(source_nodes(graph), {"AT", "CT"})
        self.assertEqual(sink_nodes(graph), {"CA"})
        self.assertEqual(branching_nodes(graph), {"TG"})

    def test_rejects_unknown_node(self) -> None:
        graph = build_dbg(["ATGC"], 3)

        with self.assertRaises(KeyError):
            in_degree(graph, "XX")
        with self.assertRaises(KeyError):
            out_degree(graph, "XX")


if __name__ == "__main__":
    unittest.main()
