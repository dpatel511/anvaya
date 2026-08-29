import io
import tempfile
import unittest
from array import array
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cli import main
from anvaya.paired_extension import (
    collect_paired_unitig_links,
    resolve_paired_extensions,
    spell_extended_paths,
)
from anvaya.sequences import reverse_complement
from anvaya.unitig_graph import CompactedUnitigGraph, build_compacted_unitig_graph


def _branch_graph() -> CompactedUnitigGraph:
    graph = CompactedUnitigGraph(k=3, sequences=["AAA", "AAC", "AAG"])
    graph.out_degrees = bytearray(6)
    graph.out_links = array("q", [-1]) * 24

    def add(source: int, target: int) -> None:
        offset = source * 4 + graph.out_degrees[source]
        graph.out_links[offset] = target
        graph.out_degrees[source] += 1

    add(0, 2)
    add(0, 4)
    add(3, 1)
    add(5, 1)
    return graph


class PairedExtensionTests(unittest.TestCase):
    def test_strong_winner_joins_reciprocal_unitig_link(self) -> None:
        graph = _branch_graph()

        successors, resolved = resolve_paired_extensions(
            graph,
            {(0, 2): 6, (0, 4): 1},
        )
        sequences = spell_extended_paths(graph, successors)

        self.assertEqual(resolved, 1)
        self.assertEqual(successors, {0: 2, 3: 1})
        self.assertEqual(len(sequences), 2)
        self.assertEqual(max(map(len, sequences)), 4)
        self.assertIn("AAAC", sequences)

    def test_tied_extensions_stop_without_joining(self) -> None:
        graph = _branch_graph()

        successors, resolved = resolve_paired_extensions(
            graph,
            {(0, 2): 5, (0, 4): 5},
        )

        self.assertEqual(resolved, 0)
        self.assertEqual(successors, {})
        self.assertEqual(len(spell_extended_paths(graph, successors)), 3)

    def test_shared_downstream_target_is_not_discriminating(self) -> None:
        graph = CompactedUnitigGraph(
            k=3,
            sequences=["AAA", "AAC", "AAG", "ACT"],
        )
        graph.out_degrees = bytearray(8)
        graph.out_links = array("q", [-1]) * 32

        def add(source: int, target: int) -> None:
            offset = source * 4 + graph.out_degrees[source]
            graph.out_links[offset] = target
            graph.out_degrees[source] += 1

        add(0, 2)
        add(0, 4)
        add(3, 1)
        add(5, 1)
        add(2, 6)
        add(4, 6)
        add(7, 3)
        add(7, 5)

        successors, resolved = resolve_paired_extensions(
            graph,
            {(0, 6): 10},
        )

        self.assertEqual(resolved, 0)
        self.assertNotIn(0, successors)

    def test_work_budget_rejects_expensive_junction(self) -> None:
        graph = _branch_graph()

        successors, resolved = resolve_paired_extensions(
            graph,
            {(0, 2): 5},
            maximum_states=1,
        )

        self.assertEqual(resolved, 0)
        self.assertEqual(successors, {})

    def test_invalid_work_budget_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "maximum_states"):
            resolve_paired_extensions(
                _branch_graph(),
                {(0, 2): 5},
                maximum_states=0,
            )

    def test_maps_one_physical_pair_to_reverse_symmetric_links(self) -> None:
        left = "AACTGGA"
        right = reverse_complement("GGACCTT")
        source = build_bidirected_dbg(
            [left, right],
            3,
            end_window=2,
            track_molecule_links=True,
            read_molecule_ids=[0, 0],
        )
        compacted = build_compacted_unitig_graph(
            source,
            track_edge_handles=True,
        )

        evidence, paired, mapped = collect_paired_unitig_links(
            source,
            compacted,
        )

        self.assertEqual(paired, 1)
        self.assertEqual(mapped, 1)
        self.assertEqual(sum(evidence.values()), 2)
        for (first, second), count in evidence.items():
            self.assertEqual(evidence[(second ^ 1, first ^ 1)], count)

    def test_cli_requires_paired_input(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            reads.write_text(">read\nAACTGGA\n", encoding="utf-8")

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(
                    [
                        "assemble",
                        "-i",
                        str(reads),
                        "--k",
                        "3",
                        "--orientation-aware",
                        "--end-window",
                        "2",
                        "--paired-unitig-extension",
                        "-o",
                        str(root / "contigs.fasta"),
                    ]
                )

    def test_cli_reports_paired_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left_path = root / "left.fasta"
            right_path = root / "right.fasta"
            left_path.write_text(">left\nAACTGGA\n", encoding="utf-8")
            right_path.write_text(
                f">right\n{reverse_complement('GGACCTT')}\n",
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                exit_code = main(
                    [
                        "assemble",
                        "-1",
                        str(left_path),
                        "-2",
                        str(right_path),
                        "--k",
                        "3",
                        "--orientation-aware",
                        "--end-window",
                        "2",
                        "--paired-unitig-extension",
                        "-o",
                        str(root / "contigs.fasta"),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertIn("paired_unitig_extension=true", stdout.getvalue())
            self.assertIn("paired_mapped_pairs=1", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
