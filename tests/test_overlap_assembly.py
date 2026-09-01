import tempfile
import unittest
from pathlib import Path

from anvaya.cli import main
from anvaya.overlap_assembly import (
    _Alignment,
    _consensus,
    _merge_contig_pass,
    assemble_overlap_contigs,
    write_overlap_correction_report,
)
from anvaya.reads import Read


class OverlapAssemblyTests(unittest.TestCase):
    def test_layout_retains_strong_internal_consensus(self) -> None:
        target = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        alternative = target[:10] + "A" + target[11:]
        members = [
            _Alignment(index, alternative, 0, 2)
            for index in range(5)
        ]

        result = _consensus(
            target,
            members,
            minimum_extension_support=5,
            dominance_ratio=4.0,
        )

        self.assertEqual(result.sequence, alternative)

    def test_polishes_supported_terminal_deamination(self) -> None:
        genome = "AACGTCCGTTGCACTGATCGATGCTAGCTACGATGGCCTATTGCA"
        truth = genome[5:40]
        damaged = "T" + truth[1:]
        reads = [Read("damaged-center", damaged)]
        reads.extend(Read(f"internal-{index}", genome[:35]) for index in range(3))
        reads.extend(Read(f"terminal-{index}", damaged) for index in range(2))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
        )

        self.assertEqual([contig.sequence for contig in contigs], [truth])
        self.assertEqual(summary.corrected_bases, 1)

    def test_records_accepted_correction_for_truth_audit(self) -> None:
        genome = "AACGTCCGTTGCACTGATCGATGCTAGCTACGATGGCCTATTGCA"
        truth = genome[5:40]
        reads = [Read("damaged-center", "T" + truth[1:])]
        reads.extend(Read(f"internal-{index}", genome[:35]) for index in range(3))
        reads.extend(Read(f"terminal-{index}", "T" + truth[1:]) for index in range(2))
        events = []

        assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
            correction_events=events,
        )

        self.assertEqual(len(events), 1)
        self.assertEqual((events[0].original_base, events[0].corrected_base), ("T", "C"))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "corrections.tsv"
            write_overlap_correction_report(events, path)
            report = path.read_text(encoding="utf-8")
        self.assertIn("sequence_before\tsequence_after", report.splitlines()[0])
        self.assertIn("overlap_contig_1", report)

    def test_polishing_does_not_change_frozen_layout(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        reads = [Read("center", center)]
        reads.extend(
            Read(f"right-{index}", center[5:] + "TTGCA")
            for index in range(5)
        )

        raw, raw_summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            damage_end_window=0,
        )
        polished, polished_summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            damage_end_window=5,
        )

        self.assertEqual(
            [len(contig.sequence) for contig in polished],
            [len(contig.sequence) for contig in raw],
        )
        self.assertEqual(polished_summary.output_contigs, raw_summary.output_contigs)
        self.assertEqual(polished_summary.contig_merges, raw_summary.contig_merges)
        self.assertEqual(polished_summary.added_bases, raw_summary.added_bases)

    def test_rejects_terminal_only_correction_support(self) -> None:
        truth = "CCGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        damaged = "T" + truth[1:]
        reads = [Read("damaged-center", damaged)]
        reads.extend(Read(f"truth-terminal-{index}", truth) for index in range(3))
        reads.extend(Read(f"damaged-terminal-{index}", damaged) for index in range(2))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
        )

        self.assertEqual([contig.sequence for contig in contigs], [damaged])
        self.assertEqual(summary.corrected_bases, 0)
        self.assertEqual(summary.correction_terminal_only, 1)

    def test_preserves_non_damage_terminal_difference(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        alternative = "G" + center[1:]
        reads = [Read("center", center)]
        reads.extend(Read(f"alternative-{index}", alternative) for index in range(3))
        reads.extend(Read(f"center-{index}", center) for index in range(2))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
        )

        self.assertEqual([contig.sequence for contig in contigs], [center])
        self.assertEqual(summary.corrected_bases, 0)

    def test_preserves_damage_like_difference_outside_terminal_window(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        alternative = center[:9] + "C" + center[10:]
        self.assertEqual(center[9], "T")
        reads = [Read("center", center)]
        reads.extend(Read(f"alternative-{index}", alternative) for index in range(3))
        reads.extend(Read(f"center-{index}", center) for index in range(2))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
        )

        self.assertEqual([contig.sequence for contig in contigs], [center])
        self.assertEqual(summary.corrected_bases, 0)

    def test_reuses_bridge_reads_then_merges_neighboring_contigs(self) -> None:
        genome = (
            "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTATTCGAGTCCGATACGTGACCTAG"
            "GATCCGTACGATCGTACGA"
        )
        reads = [Read("left", genome[:35]), Read("right", genome[30:65])]
        reads.extend(
            Read(f"bridge-{index}", genome[15:50])
            for index in range(5)
        )

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_rounds=2,
            maximum_contig_iterations=2,
        )

        self.assertEqual([contig.sequence for contig in contigs], [genome[:65]])
        self.assertEqual(summary.clusters, 2)
        self.assertEqual(summary.contig_merges, 1)

    def test_contig_phase_merges_a_unique_best_overlap(self) -> None:
        left = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        right = left[15:] + "TTGCAACGT"

        merged, merges, rounds, added, ambiguous = _merge_contig_pass(
            [Read("left", left), Read("right", right)],
            anchor_k=7,
            anchors_per_read=8,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=2,
            minimum_overlap=20,
            minimum_identity=0.90,
            minimum_ry_identity=0.99,
            maximum_rounds=3,
            minimum_overlap_margin=3,
        )

        self.assertEqual([contig.sequence for contig in merged], [left + "TTGCAACGT"])
        self.assertEqual((merges, rounds, added, ambiguous), (1, 1, 9, 0))

    def test_contig_phase_stops_at_equal_competing_extensions(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        first = center[15:] + "TTGCAACGT"
        second = center[15:] + "CCATGGTAC"

        merged, merges, _, _, ambiguous = _merge_contig_pass(
            [Read("center", center), Read("first", first), Read("second", second)],
            anchor_k=7,
            anchors_per_read=8,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=2,
            minimum_overlap=20,
            minimum_identity=0.90,
            minimum_ry_identity=0.99,
            maximum_rounds=3,
            minimum_overlap_margin=3,
        )

        self.assertNotIn(center + "TTGCAACGT", [contig.sequence for contig in merged])
        self.assertNotIn(center + "CCATGGTAC", [contig.sequence for contig in merged])
        self.assertEqual(merges, 0)
        self.assertGreaterEqual(ambiguous, 1)

    def test_emits_a_directly_extended_contig(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        extension = "TTGCA"
        reads = [Read("center", center)]
        reads.extend(
            Read(f"right-{index}", center[5:] + extension)
            for index in range(5)
        )

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
        )

        self.assertEqual([contig.sequence for contig in contigs], [center + extension])
        self.assertEqual(summary.output_contigs, 1)
        self.assertEqual(summary.clustered_reads, 6)
        self.assertEqual(summary.added_bases, 5)

    def test_does_not_emit_unsupported_singletons(self) -> None:
        reads = [Read("first", "AACCGGTTAACCGG"), Read("second", "TTTTGGGGCCCCAA")]

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=5,
            minimum_overlap=10,
        )

        self.assertEqual(contigs, [])
        self.assertEqual(summary.clusters, 0)

    def test_reindexes_an_extended_end_for_another_round(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        first = center + "TTGCA"
        reads = [Read("center", center)]
        reads.extend(
            Read(f"round-one-{index}", center[5:] + "TTGCA")
            for index in range(5)
        )
        reads.extend(
            Read(f"round-two-{index}", first[20:] + "AACGT")
            for index in range(5)
        )

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_rounds=2,
        )

        self.assertEqual(
            [contig.sequence for contig in contigs],
            [center + "TTGCAAACGT"],
        )
        self.assertEqual(summary.extension_rounds, 2)
        self.assertEqual(summary.added_bases, 10)

    def test_separates_ry_inconsistent_clusters(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        overlap = list(center[5:])
        overlap[10] = "C" if overlap[10] in "AG" else "A"
        conflicting = "".join(overlap) + "TTGCA"
        reads = [Read("center", center)]
        reads.extend(Read(f"truth-{index}", center) for index in range(5))
        reads.extend(Read(f"member-{index}", conflicting) for index in range(5))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
        )

        self.assertEqual(
            {contig.sequence for contig in contigs},
            {center, conflicting},
        )
        self.assertEqual(summary.clusters, 2)

    def test_cli_is_separate_from_dbg_assembly(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "reads.fasta"
            output_path = root / "contigs.fasta"
            records = [("center", center)] + [
                (f"right-{index}", center[5:] + "TTGCA")
                for index in range(5)
            ]
            input_path.write_text(
                "".join(f">{name}\n{sequence}\n" for name, sequence in records),
                encoding="utf-8",
            )

            result = main([
                "overlap-assemble",
                "-i",
                str(input_path),
                "-o",
                str(output_path),
                "--min-output-length",
                "1",
            ])

            self.assertEqual(result, 0)
            self.assertIn(center + "TTGCA", output_path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
