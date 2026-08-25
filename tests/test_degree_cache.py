"""Tests for linear-time graph degree calculation."""

import unittest

from anvaya.graph import build_dbg, node_degrees


class NodeDegreesTests(unittest.TestCase):
    def test_calculates_all_node_degrees(self) -> None:
        graph = build_dbg(["AACG", "AATG"], 3)

        in_degrees, out_degrees = node_degrees(graph)

        self.assertEqual(
            in_degrees, {"AA": 0, "AC": 1, "CG": 1, "AT": 1, "TG": 1}
        )
        self.assertEqual(
            out_degrees, {"AA": 2, "AC": 1, "CG": 0, "AT": 1, "TG": 0}
        )


if __name__ == "__main__":
    unittest.main()
