import io
import tempfile
import unittest
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.assembly import assemble
from anvaya.cli import main
from anvaya.graph import build_dbg


class MinimumCountGraphTests(unittest.TestCase):
    def test_default_preserves_unfiltered_graph(self) -> None:
        self.assertEqual(
            build_dbg(["ATGC", "ATGA"], 3),
            build_dbg(["ATGC", "ATGA"], 3, min_count=1),
        )

    def test_removes_kmers_below_minimum_count(self) -> None:
        graph = build_dbg(["ATGC", "ATGA"], 3, min_count=2)

        self.assertEqual(
            graph,
            {
                "AT": Counter({"TG": 2}),
                "TG": Counter(),
            },
        )

    def test_preserves_retained_edge_multiplicity(self) -> None:
        graph = build_dbg(["ATGC", "ATGC", "ATGA"], 3, min_count=2)

        self.assertEqual(graph["AT"]["TG"], 3)
        self.assertEqual(graph["TG"]["GC"], 2)
        self.assertNotIn("GA", graph)

    def test_allows_an_empty_filtered_graph(self) -> None:
        self.assertEqual(build_dbg(["ATGC"], 3, min_count=2), {})

    def test_rejects_invalid_minimum_count(self) -> None:
        with self.assertRaises(ValueError):
            build_dbg(["ATGC"], 3, min_count=0)
        with self.assertRaises(TypeError):
            build_dbg(["ATGC"], 3, min_count=True)

    def test_assembly_uses_minimum_count(self) -> None:
        self.assertEqual(assemble(["ATGC", "ATGA"], 3, min_count=2), ["ATG"])


class MinimumCountCliTests(unittest.TestCase):
    def test_cli_reports_requested_minimum_count(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            output = root / "unitigs.fasta"
            reads.write_text(
                ">read-1\nATGC\n>read-2\nATGA\n",
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
                        "-o",
                        str(output),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertIn("min_count=2", stdout.getvalue())
            self.assertEqual(
                output.read_text(encoding="utf-8"),
                ">unitig_1 length=3\nATG\n",
            )


if __name__ == "__main__":
    unittest.main()
