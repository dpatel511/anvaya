import io
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.bubbles import find_simple_bubbles
from anvaya.cli import main
from anvaya.sequences import reverse_complement


class SimpleBubbleTests(unittest.TestCase):
    def test_finds_substitution_bubble_and_support(self) -> None:
        graph = build_bidirected_dbg(
            ["GCTTGTTCCGGA"] * 10 + ["GCTTATTCCGGA"] * 3,
            4,
        )
        before_degrees = bytes(graph.out_degrees)
        before_observations = graph.observations

        bubbles = find_simple_bubbles(graph)

        self.assertEqual(len(bubbles), 1)
        self.assertEqual(len(bubbles[0].paths), 2)
        self.assertEqual(
            sorted(path.minimum_edge_support for path in bubbles[0].paths),
            [3, 10],
        )
        self.assertEqual(bytes(graph.out_degrees), before_degrees)
        self.assertEqual(graph.observations, before_observations)

    def test_linear_graph_has_no_bubble(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA"] * 5, 3)
        self.assertEqual(find_simple_bubbles(graph), [])

    def test_limit_excludes_long_bubble(self) -> None:
        graph = build_bidirected_dbg(
            ["GCTTGTTCCGGA"] * 10 + ["GCTTATTCCGGA"] * 3,
            4,
        )
        self.assertEqual(find_simple_bubbles(graph, max_edges=1), [])
        self.assertEqual(len(find_simple_bubbles(graph, max_edges=8)), 1)

    def test_result_is_deterministic(self) -> None:
        graph = build_bidirected_dbg(
            ["GCTTGTTCCGGA"] * 10 + ["GCTTATTCCGGA"] * 3,
            4,
        )
        self.assertEqual(find_simple_bubbles(graph), find_simple_bubbles(graph))

    def test_reverse_complements_do_not_duplicate_bubble(self) -> None:
        major = "GCTTGTTCCGGA"
        minor = "GCTTATTCCGGA"
        graph = build_bidirected_dbg(
            [major, reverse_complement(major)] * 5
            + [minor, reverse_complement(minor)] * 2,
            4,
        )

        bubbles = find_simple_bubbles(graph)

        self.assertEqual(len(bubbles), 1)
        self.assertEqual(
            sorted(path.minimum_edge_support for path in bubbles[0].paths),
            [4, 10],
        )

    def test_rejects_invalid_limit(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA"], 3)
        with self.assertRaisesRegex(ValueError, "at least 1"):
            find_simple_bubbles(graph, max_edges=0)
        with self.assertRaisesRegex(TypeError, "integer"):
            find_simple_bubbles(graph, max_edges=1.5)  # type: ignore[arg-type]


class BubbleCliTests(unittest.TestCase):
    def test_cli_reports_bubbles_without_changing_unitigs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            plain_output = root / "plain.fasta"
            detected_output = root / "detected.fasta"
            records = ["GCTTGTTCCGGA"] * 10 + ["GCTTATTCCGGA"] * 3
            reads.write_text(
                "".join(
                    f">read_{index}\n{read}\n"
                    for index, read in enumerate(records)
                ),
                encoding="utf-8",
            )

            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "4",
                        "--orientation-aware", "-o", str(plain_output),
                    ]
                )
            stdout = io.StringIO()
            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "4",
                        "--orientation-aware", "--detect-bubbles",
                        "-o", str(detected_output),
                    ]
                )

            self.assertIn("bubble_detection=true", stdout.getvalue())
            self.assertIn("bubbles_detected=1", stdout.getvalue())
            self.assertEqual(
                detected_output.read_bytes(),
                plain_output.read_bytes(),
            )

    def test_detection_requires_orientation_aware_graph(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            reads.write_text(">read\nAACTGGA\n", encoding="utf-8")

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "3",
                        "--detect-bubbles", "-o", str(root / "unitigs.fasta"),
                    ]
                )


if __name__ == "__main__":
    unittest.main()
