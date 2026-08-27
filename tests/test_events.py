import csv
import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.bubbles import find_simple_bubbles
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.cli import main
from anvaya.events import write_event_report
from anvaya.incomplete_branches import find_incomplete_branch_candidates


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class EventReportTests(unittest.TestCase):
    def test_reports_bubble_sequences_support_and_substitution(self) -> None:
        graph = build_bidirected_dbg(
            ["GCTTGTTCCGGA"] * 10 + ["GCTTATTCCGGA"] * 3,
            4,
            end_window=10,
        )
        bubbles = find_simple_bubbles(graph)

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "events.tsv"
            summary = write_event_report(graph, [], bubbles, output)
            rows = _read_rows(output)

        self.assertEqual(summary.bubbles, 1)
        self.assertEqual(summary.paths, 2)
        self.assertEqual({row["event_type"] for row in rows}, {"bubble"})
        self.assertTrue(all(None not in row for row in rows))
        reference = next(row for row in rows if row["reference_path"] == "true")
        alternative = next(
            row for row in rows if row["reference_path"] == "false"
        )
        self.assertEqual(reference["minimum_edge_support"], "10")
        self.assertEqual(alternative["minimum_edge_support"], "3")
        self.assertNotEqual(reference["sequence"], alternative["sequence"])
        self.assertRegex(
            alternative["substitutions"],
            r"\d+:[ACGT]>[ACGT]",
        )
        self.assertEqual(alternative["damage_compatible"], "true")

    def test_reports_weak_tip_before_graph_cleaning(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAA"] * 10 + ["AAGCCTAAA"],
            5,
            end_window=3,
            track_molecule_links=True,
        )
        before = bytes(graph.out_degrees)
        tips = find_weak_tip_candidates(graph)

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "events.tsv"
            summary = write_event_report(graph, tips, [], output)
            rows = _read_rows(output)

        self.assertEqual(summary.tips, 1)
        self.assertEqual(summary.matched_tips, 1)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["event_type"], "tip")
        self.assertTrue(rows[0]["sequence"])
        self.assertGreater(int(rows[0]["left_terminal"]), 0)
        self.assertEqual(rows[0]["backbone_matched"], "true")
        self.assertEqual(rows[0]["tip_match_sequence"], "AGCCTAAA")
        self.assertEqual(rows[0]["backbone_sequence"], "AGCCCAAA")
        self.assertEqual(rows[0]["substitutions"], "5:C>T")
        self.assertEqual(rows[0]["damage_compatible"], "true")
        self.assertEqual(rows[0]["ry_identity"], "1.000000")
        self.assertEqual(rows[0]["tip_classification"], "damage-like")
        self.assertGreater(float(rows[0]["tip_damage_score"]), 0.8)
        self.assertEqual(
            rows[0]["substitution_terminal_observations"],
            "1",
        )
        self.assertEqual(rows[0]["mean_damage_distance"], "1.000000")
        self.assertEqual(rows[0]["molecule_links_collected"], "true")
        self.assertEqual(rows[0]["joint_molecule_observations"], "1")
        self.assertEqual(rows[0]["joint_molecule_fraction"], "1.000000")
        self.assertEqual(bytes(graph.out_degrees), before)

    def test_reports_classified_incomplete_branch(self) -> None:
        graph = build_bidirected_dbg(
            ["AAGCCCAAA"] * 10 + ["AAGCCTAAA", "AAGCCTAAC"],
            5,
            end_window=3,
        )
        before = bytes(graph.out_degrees)
        candidates = find_incomplete_branch_candidates(graph)

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "events.tsv"
            summary = write_event_report(
                graph,
                [],
                [],
                output,
                incomplete_branches=candidates,
            )
            rows = _read_rows(output)

        self.assertEqual(summary.incomplete_branches, 1)
        self.assertEqual(summary.matched_incomplete_branches, 1)
        self.assertEqual(summary.paths, 1)
        self.assertEqual(rows[0]["event_type"], "incomplete_branch")
        self.assertEqual(rows[0]["backbone_matched"], "true")
        self.assertEqual(rows[0]["substitutions"], "5:C>T")
        self.assertEqual(rows[0]["tip_classification"], "damage-like")
        self.assertEqual(bytes(graph.out_degrees), before)

    def test_requires_read_end_evidence(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA"], 3)

        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "requires read-end"):
                write_event_report(
                    graph,
                    [],
                    [],
                    Path(directory) / "events.tsv",
                )

    def test_profile_failure_does_not_write_event_report(self) -> None:
        graph = build_bidirected_dbg(["AACTGGA"], 3, end_window=2)

        with tempfile.TemporaryDirectory() as directory:
            events = Path(directory) / "events.tsv"
            profile = Path(directory) / "profile.json"
            with self.assertRaisesRegex(ValueError, "requires molecule links"):
                write_event_report(
                    graph,
                    [],
                    [],
                    events,
                    damage_profile_path=profile,
                )

            self.assertFalse(events.exists())
            self.assertFalse(profile.exists())


class EventReportCliTests(unittest.TestCase):
    def test_cli_reports_incomplete_branch_without_changing_assembly(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            baseline = root / "baseline.fasta"
            reported = root / "reported.fasta"
            report = root / "events.tsv"
            records = ["AAGCCCAAA"] * 10 + ["AAGCCTAAA", "AAGCCTAAC"]
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
                        "assemble", "-i", str(reads), "--k", "5",
                        "--orientation-aware", "--end-window", "3",
                        "-o", str(baseline),
                    ]
                )
            stdout = io.StringIO()
            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "5",
                        "--orientation-aware", "--end-window", "3",
                        "--event-report", str(report),
                        "-o", str(reported),
                    ]
                )

            self.assertIn("reported_incomplete_branches=1", stdout.getvalue())
            self.assertIn(
                "reported_incomplete_branch_matches=1",
                stdout.getvalue(),
            )
            self.assertEqual(reported.read_bytes(), baseline.read_bytes())

    def test_cli_writes_report_without_changing_assembly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            baseline = root / "baseline.fasta"
            reported = root / "reported.fasta"
            report = root / "events.tsv"
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
                        "--orientation-aware", "--end-window", "10",
                        "-o", str(baseline),
                    ]
                )
            stdout = io.StringIO()
            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "4",
                        "--orientation-aware", "--end-window", "10",
                        "--event-report", str(report),
                        "-o", str(reported),
                    ]
                )

            self.assertTrue(report.exists())
            self.assertIn("reported_bubbles=1", stdout.getvalue())
            self.assertIn("reported_tip_matches=0", stdout.getvalue())
            self.assertIn(f"event_report={report}", stdout.getvalue())
            self.assertEqual(reported.read_bytes(), baseline.read_bytes())

    def test_event_report_requires_positive_end_window(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            reads.write_text(">read\nAACTGGA\n", encoding="utf-8")

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "3",
                        "--orientation-aware",
                        "--event-report", str(root / "events.tsv"),
                        "-o", str(root / "unitigs.fasta"),
                    ]
                )

    def test_cli_writes_damage_profile_without_changing_assembly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            baseline = root / "baseline.fasta"
            reported = root / "reported.fasta"
            events = root / "events.tsv"
            profile = root / "damage-profile.json"
            records = ["AAGCCCAAA"] * 10 + ["AAGCCTAAA"]
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
                        "assemble", "-i", str(reads), "--k", "5",
                        "--orientation-aware", "--end-window", "3",
                        "-o", str(baseline),
                    ]
                )
            stdout = io.StringIO()
            with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "5",
                        "--orientation-aware", "--end-window", "3",
                        "--event-report", str(events),
                        "--damage-profile-report", str(profile),
                        "-o", str(reported),
                    ]
                )

            payload = json.loads(profile.read_text(encoding="utf-8"))
            self.assertEqual(payload["eligible_loci"], 1)
            self.assertIn(f"damage_profile_report={profile}", stdout.getvalue())
            self.assertEqual(reported.read_bytes(), baseline.read_bytes())

    def test_damage_profile_requires_event_report(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = root / "reads.fasta"
            reads.write_text(">read\nAACTGGA\n", encoding="utf-8")

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main(
                    [
                        "assemble", "-i", str(reads), "--k", "3",
                        "--orientation-aware", "--end-window", "3",
                        "--damage-profile-report", str(root / "profile.json"),
                        "-o", str(root / "unitigs.fasta"),
                    ]
                )


if __name__ == "__main__":
    unittest.main()
