import unittest

from anvaya.graph import build_dbg
from anvaya.metrics import GraphSummary, summarize_graph


class SummarizeGraphTests(unittest.TestCase):
    def test_summarizes_linear_graph(self) -> None:
        graph = build_dbg(["ATGCA"], 3)

        self.assertEqual(
            summarize_graph(graph),
            GraphSummary(
                nodes=4,
                edges=3,
                observations=3,
                sources=1,
                sinks=1,
                branching_nodes=0,
                unitigs=1,
            ),
        )

    def test_summarizes_branching_graph(self) -> None:
        graph = build_dbg(["ATGCA", "ATGGA"], 3)

        self.assertEqual(
            summarize_graph(graph),
            GraphSummary(
                nodes=6,
                edges=5,
                observations=6,
                sources=1,
                sinks=2,
                branching_nodes=1,
                unitigs=3,
            ),
        )

    def test_summarizes_empty_graph(self) -> None:
        self.assertEqual(
            summarize_graph({}),
            GraphSummary(0, 0, 0, 0, 0, 0, 0),
        )


if __name__ == "__main__":
    unittest.main()
