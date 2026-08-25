import io
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.bidirected import (
    build_bidirected_dbg,
    extract_bidirected_unitigs,
    oriented_sequence,
    summarize_bidirected_graph,
)
from anvaya.cli import main
from anvaya.sequences import canonical_sequence, reverse_complement


class BidirectedGraphTests(unittest.TestCase):
    def test_read_and_reverse_complement_share_physical_edges(self) -> None:
        read = "AACTGGA"
        graph = build_bidirected_dbg([read, reverse_complement(read)], 3)

        self.assertEqual(len(graph.edge_support), 5)
        self.assertTrue(
            all(
                support.forward == 1 and support.reverse == 1
                for support in graph.edge_support.values()
            )
        )

    def test_combined_strand_support_passes_filter(self) -> None:
        read = "AACTGGA"
        graph = build_bidirected_dbg(
            [read, reverse_complement(read)], 3, min_count=2
        )

        self.assertEqual(len(graph.edge_support), 5)
        unitigs = extract_bidirected_unitigs(graph)
        self.assertEqual(len(unitigs), 1)
        self.assertEqual(canonical_sequence(unitigs[0]), canonical_sequence(read))

    def test_spells_one_orientation_of_linear_path(self) -> None:
        read = "AACTGGA"
        graph = build_bidirected_dbg([read], 3)

        unitigs = extract_bidirected_unitigs(graph)

        self.assertEqual(len(unitigs), 1)
        self.assertEqual(canonical_sequence(unitigs[0]), canonical_sequence(read))

    def test_records_palindromic_kmer_support(self) -> None:
        graph = build_bidirected_dbg(["ATAT"], 4)

        support = next(iter(graph.edge_support.values()))
        self.assertEqual(support.palindromic, 1)
        self.assertEqual(support.total, 1)

    def test_oriented_handles_spell_reverse_complements(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA"], 3)

        for node_id, sequence in enumerate(graph.sequences):
            self.assertEqual(oriented_sequence(graph, node_id << 1), sequence)
            self.assertEqual(
                oriented_sequence(graph, (node_id << 1) | 1),
                reverse_complement(sequence),
            )

    def test_summary_counts_physical_nodes_and_edges(self) -> None:
        read = "AACTGGA"
        graph = build_bidirected_dbg([read, reverse_complement(read)], 3)

        summary = summarize_bidirected_graph(graph)

        self.assertEqual(summary.edges, 5)
        self.assertEqual(summary.observations, 10)
        self.assertEqual(summary.unitigs, 1)


class BidirectedCliTests(unittest.TestCase):
    def test_orientation_aware_cli_writes_single_unitig(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            output = root / "unitigs.fasta"
            read = "AACTGGA"
            reads.write_text(
                f">forward\n{read}\n>reverse\n{reverse_complement(read)}\n",
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                exit_code = main(
                    [
                        "assemble",
                        "-i",
                        str(reads),
                        "--k",
                        "3",
                        "--min-count",
                        "2",
                        "--orientation-aware",
                        "-o",
                        str(output),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertIn("orientation_aware=true", stdout.getvalue())
            self.assertIn("unitigs=1", stdout.getvalue())
            self.assertTrue(output.exists())


if __name__ == "__main__":
    unittest.main()
