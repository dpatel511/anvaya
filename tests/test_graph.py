import unittest
from collections import Counter

from anvaya.graph import build_dbg


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


if __name__ == "__main__":
    unittest.main()
