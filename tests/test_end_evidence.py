import io
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.bidirected import (
    _outgoing_edges,
    build_bidirected_dbg,
    extract_bidirected_unitigs,
)
from anvaya.cli import main
from anvaya.cleaning import remove_weak_tips
from anvaya.sequences import reverse_complement


class EndEvidenceGraphTests(unittest.TestCase):
    def test_counts_terminal_and_internal_observations(self) -> None:
        graph = build_bidirected_dbg(["AACCTGGAA"], 3, end_window=2)

        self.assertEqual(graph.terminal_observations, 4)
        self.assertEqual(
            sum(graph.internal_support),
            graph.observations - graph.terminal_observations,
        )
        self.assertEqual(
            len(graph.left_end_support),
            graph.edge_count * graph.end_window,
        )
        self.assertEqual(
            len(graph.right_end_support),
            graph.edge_count * graph.end_window,
        )

    def test_normalizes_reverse_complement_end_evidence(self) -> None:
        read = "AACCTGGAA"
        graph = build_bidirected_dbg(
            [read, reverse_complement(read)],
            3,
            end_window=2,
        )

        self.assertEqual(graph.terminal_observations, 8)
        self.assertEqual(sum(graph.left_end_support), 4)
        self.assertEqual(sum(graph.right_end_support), 4)
        self.assertTrue(
            all(count % 2 == 0 for count in graph.left_end_support)
        )
        self.assertTrue(
            all(count % 2 == 0 for count in graph.right_end_support)
        )

    def test_records_five_prime_damage_path_as_terminal(self) -> None:
        graph = build_bidirected_dbg(
            ["CAACTGGA"] * 10 + ["TAACTGGA"] * 3,
            3,
            end_window=3,
        )
        damaged_edges = [
            edge_id
            for edge_id in range(graph.edge_count)
            if graph.strand_support(edge_id).total == 3
        ]

        self.assertEqual(len(damaged_edges), 1)
        evidence = graph.end_support(damaged_edges[0])
        self.assertEqual(evidence.left[0], 3)
        self.assertEqual(sum(evidence.right), 0)
        self.assertEqual(evidence.internal, 0)

    def test_links_terminal_edge_observation_to_source_read(self) -> None:
        graph = build_bidirected_dbg(
            ["CAACTGGA", "TAACTGGA"],
            3,
            end_window=3,
            track_molecule_links=True,
        )
        damaged_edge = next(
            edge_id
            for edge_id in range(graph.edge_count)
            if graph.strand_support(edge_id).total == 1
            and graph.end_support(edge_id).left[0] == 1
        )

        links = graph.molecule_end_links(damaged_edge)

        self.assertIn(
            (1, "left", 0),
            {
                (link.read_index, link.end, link.distance)
                for link in links
            },
        )
        self.assertEqual(len(graph.molecule_link_offsets), graph.edge_count + 1)

    def test_molecule_links_are_opt_in(self) -> None:
        graph = build_bidirected_dbg(["AACCTGGAA"], 3, end_window=2)

        with self.assertRaisesRegex(ValueError, "were not collected"):
            graph.molecule_end_links(0)

    def test_records_three_prime_damage_path_as_terminal(self) -> None:
        graph = build_bidirected_dbg(
            ["AACTGGAG"] * 10 + ["AACTGGAA"] * 3,
            3,
            end_window=3,
        )
        damaged_edges = [
            edge_id
            for edge_id in range(graph.edge_count)
            if graph.strand_support(edge_id).total == 3
        ]

        self.assertEqual(len(damaged_edges), 1)
        evidence = graph.end_support(damaged_edges[0])
        self.assertEqual(evidence.right[0], 3)
        self.assertEqual(sum(evidence.left), 0)
        self.assertEqual(evidence.internal, 0)

    def test_evidence_collection_does_not_change_unitigs(self) -> None:
        reads = ["GCTTGTTCCGGA"] * 10 + ["GCTTATTCCGGA"] * 3
        baseline = build_bidirected_dbg(reads, 4)
        annotated = build_bidirected_dbg(reads, 4, end_window=3)

        self.assertEqual(
            sorted(extract_bidirected_unitigs(annotated)),
            sorted(extract_bidirected_unitigs(baseline)),
        )
        self.assertEqual(annotated.observations, baseline.observations)
        self.assertEqual(annotated.edge_count, baseline.edge_count)

    def test_records_palindromic_terminal_support_as_ambiguous(self) -> None:
        graph = build_bidirected_dbg(["ATAT"], 4, end_window=2)

        evidence = graph.end_support(0)
        self.assertEqual(evidence.ambiguous_terminal, 1)
        self.assertEqual(evidence.internal, 0)
        self.assertEqual(sum(evidence.left), 0)
        self.assertEqual(sum(evidence.right), 0)

    def test_overlapping_end_windows_count_each_observation_once(self) -> None:
        graph = build_bidirected_dbg(["AACTG"], 3, end_window=3)

        self.assertEqual(graph.observations, 3)
        self.assertEqual(graph.terminal_observations, 3)
        self.assertEqual(sum(graph.internal_support), 0)

    def test_tip_removal_updates_active_terminal_total(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAT"] * 10 + ["AAGCCTAAA"],
            5,
            end_window=2,
        )
        before = graph.terminal_observations

        remove_weak_tips(graph)

        self.assertLess(graph.terminal_observations, before)
        active_edges = {
            edge_id
            for handle in range(graph.node_count * 2)
            for edge_id in _outgoing_edges(graph, handle)
        }
        self.assertEqual(
            graph.terminal_observations,
            sum(
                graph.strand_support(edge_id).total
                - graph.internal_support[edge_id]
                for edge_id in active_edges
            ),
        )

    def test_rejects_invalid_end_window(self) -> None:
        with self.assertRaisesRegex(ValueError, "must not be negative"):
            build_bidirected_dbg(["AACTGGA"], 3, end_window=-1)
        with self.assertRaisesRegex(TypeError, "integer"):
            build_bidirected_dbg(
                ["AACTGGA"],
                3,
                end_window=1.5,  # type: ignore[arg-type]
            )
        with self.assertRaisesRegex(ValueError, "require read-end"):
            build_bidirected_dbg(
                ["AACTGGA"],
                3,
                track_molecule_links=True,
            )
        with self.assertRaisesRegex(TypeError, "must be a boolean"):
            build_bidirected_dbg(
                ["AACTGGA"],
                3,
                end_window=2,
                track_molecule_links=1,  # type: ignore[arg-type]
            )

    def test_requires_evidence_before_access(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA"], 3)

        with self.assertRaisesRegex(ValueError, "was not collected"):
            graph.end_support(0)


class EndEvidenceCliTests(unittest.TestCase):
    def test_cli_reports_terminal_evidence_without_changing_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            baseline_output = root / "baseline.fasta"
            evidence_output = root / "evidence.fasta"
            reads.write_text(">read\nAACCTGGAA\n", encoding="utf-8")

            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "3",
                        "--orientation-aware", "-o", str(baseline_output),
                    ]
                )
            stdout = io.StringIO()
            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "3",
                        "--orientation-aware", "--end-window", "2",
                        "-o", str(evidence_output),
                    ]
                )

            self.assertIn("end_window=2", stdout.getvalue())
            self.assertIn("terminal_observations=4", stdout.getvalue())
            self.assertEqual(
                evidence_output.read_bytes(),
                baseline_output.read_bytes(),
            )

    def test_end_window_requires_orientation_aware_graph(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            reads.write_text(">read\nAACTGGA\n", encoding="utf-8")

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "3",
                        "--end-window", "2",
                        "-o", str(root / "unitigs.fasta"),
                    ]
                )


if __name__ == "__main__":
    unittest.main()
