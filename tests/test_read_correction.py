import csv
import tempfile
import unittest
from pathlib import Path

from anvaya.read_correction import write_correction_report
from anvaya.reads import Read


class ReadCorrectionTests(unittest.TestCase):
    def test_reports_unique_complete_solid_kmer_rescue(self) -> None:
        reads = [
            Read(f"truth-{index}", "GGGAACCG", (35,) * 8)
            for index in range(3)
        ]
        reads.append(Read("error", "GGGAATCG", (35, 35, 35, 35, 35, 8, 35, 35)))

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "corrections.tsv"
            summary = write_correction_report(reads, path, k=3, min_count=2)
            with path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(summary.would_correct, 1)
        candidate = next(row for row in rows if row["decision"] == "would_correct")
        self.assertEqual(candidate["observed_base"], "T")
        self.assertEqual(candidate["proposed_base"], "C")

    def test_protects_terminal_damage_compatible_reversal(self) -> None:
        reads = [Read(f"truth-{index}", "AACCG", (35,) * 5) for index in range(3)]
        reads.append(Read("damage", "AATCG", (35, 35, 8, 35, 35)))

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "corrections.tsv"
            summary = write_correction_report(
                reads, path, k=3, min_count=2, end_window=3
            )

        self.assertEqual(summary.would_correct, 0)
        self.assertEqual(summary.protected_damage, 1)

    def test_fasta_reads_are_not_assigned_corrections(self) -> None:
        reads = [Read("read", "AACCG")]
        with tempfile.TemporaryDirectory() as directory:
            summary = write_correction_report(
                reads,
                Path(directory) / "corrections.tsv",
                k=3,
                min_count=2,
            )

        self.assertEqual(summary.reads_without_quality, 1)
        self.assertEqual(summary.would_correct, 0)

    def test_evenly_spaced_sample_keeps_full_spectrum(self) -> None:
        reads = [Read(f"read-{index}", "AACCG", (35,) * 5) for index in range(10)]
        with tempfile.TemporaryDirectory() as directory:
            summary = write_correction_report(
                reads,
                Path(directory) / "corrections.tsv",
                k=3,
                min_count=2,
                maximum_reads=3,
            )

        self.assertEqual(summary.reads, 10)
        self.assertEqual(summary.audited_reads, 3)

    def test_ambiguous_kmer_windows_are_zero_support(self) -> None:
        reads = [Read("ambiguous", "AANC", (35, 8, 35, 35))]
        with tempfile.TemporaryDirectory() as directory:
            summary = write_correction_report(
                reads,
                Path(directory) / "corrections.tsv",
                k=3,
                min_count=2,
            )

        self.assertEqual(summary.low_quality_bases, 1)


if __name__ == "__main__":
    unittest.main()
