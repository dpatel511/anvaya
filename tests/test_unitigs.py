import unittest

from anvaya.graph import build_dbg
from anvaya.unitigs import extract_unitigs


class ExtractUnitigsTests(unittest.TestCase):
    def test_extracts_linear_unitig(self) -> None:
        graph = build_dbg(["ATGCA"], 3)

        self.assertEqual(extract_unitigs(graph), ["ATGCA"])

    def test_stops_and_restarts_at_outgoing_branch(self) -> None:
        graph = build_dbg(["ATGCA", "ATGGA"], 3)

        self.assertEqual(extract_unitigs(graph), ["ATG", "TGCA", "TGGA"])

    def test_stops_and_restarts_at_merge(self) -> None:
        graph = build_dbg(["ATGCA", "CTGCA"], 3)

        self.assertEqual(extract_unitigs(graph), ["ATG", "CTG", "TGCA"])

    def test_edge_multiplicity_does_not_duplicate_unitig(self) -> None:
        graph = build_dbg(["ATGCA", "ATGCA"], 3)

        self.assertEqual(extract_unitigs(graph), ["ATGCA"])

    def test_extracts_isolated_cycle(self) -> None:
        graph = build_dbg(["ATAT"], 3)

        self.assertEqual(extract_unitigs(graph), ["ATAT"])

    def test_extracts_self_loop(self) -> None:
        graph = build_dbg(["AAA"], 3)

        self.assertEqual(extract_unitigs(graph), ["AAA"])

    def test_empty_graph_has_no_unitigs(self) -> None:
        self.assertEqual(extract_unitigs({}), [])


if __name__ == "__main__":
    unittest.main()
