import unittest

from anvaya.bidirected import build_bidirected_dbg, extract_bidirected_unitigs
from anvaya.sequences import canonical_sequence, reverse_complement
from anvaya.unitig_graph import build_compacted_unitig_graph


class CompactedUnitigGraphTests(unittest.TestCase):
    def test_preserves_linear_unitig_and_support(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA"], 3, end_window=2)

        compacted = build_compacted_unitig_graph(graph)

        self.assertEqual(compacted.unitig_count, 1)
        self.assertEqual(
            canonical_sequence(compacted.sequences[0]),
            canonical_sequence("AACTGGA"),
        )
        support = compacted.support(0)
        self.assertEqual(support.edge_count, 5)
        self.assertEqual(support.observations, 5)
        self.assertEqual(support.minimum_edge_support, 1)
        self.assertEqual(support.mean_edge_support, 1.0)
        self.assertEqual(support.terminal_observations, 4)
        self.assertEqual(support.internal_observations, 1)

    def test_matches_existing_unitigs_on_branching_graph(self) -> None:
        graph = build_bidirected_dbg(
            ["AACTGGA", "AACTAGA", "TCCAGTT"], 3, min_count=1
        )
        expected = sorted(
            canonical_sequence(sequence)
            for sequence in extract_bidirected_unitigs(graph)
        )

        compacted = build_compacted_unitig_graph(graph)

        observed = sorted(
            canonical_sequence(sequence) for sequence in compacted.sequences
        )
        self.assertEqual(observed, expected)

    def test_links_are_reverse_symmetric(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA", "AACTAGA"], 3)
        compacted = build_compacted_unitig_graph(graph)

        self.assertGreater(compacted.oriented_link_count, 0)
        for source in range(compacted.unitig_count * 2):
            for target in compacted.outgoing(source):
                self.assertIn(source ^ 1, compacted.outgoing(target ^ 1))

    def test_oriented_sequences_are_reverse_complements(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA"], 3)
        compacted = build_compacted_unitig_graph(graph)

        self.assertEqual(
            compacted.oriented_sequence(1),
            reverse_complement(compacted.oriented_sequence(0)),
        )

    def test_preserves_isolated_cycle(self) -> None:
        graph = build_bidirected_dbg(["ATAT"], 3)

        compacted = build_compacted_unitig_graph(graph)

        self.assertEqual(
            sorted(map(canonical_sequence, compacted.sequences)),
            sorted(map(canonical_sequence, extract_bidirected_unitigs(graph))),
        )


if __name__ == "__main__":
    unittest.main()
