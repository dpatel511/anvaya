import unittest
from collections import Counter

from anvaya.graph import build_dbg


class AmbiguousGraphTests(unittest.TestCase):
    def test_skips_kmers_containing_n(self) -> None:
        graph = build_dbg(["ACGNTA"], 3)

        self.assertEqual(
            graph,
            {
                "AC": Counter({"CG": 1}),
                "CG": Counter(),
            },
        )

    def test_all_ambiguous_kmers_produce_empty_graph(self) -> None:
        self.assertEqual(build_dbg(["AN"], 2), {})


if __name__ == "__main__":
    unittest.main()
