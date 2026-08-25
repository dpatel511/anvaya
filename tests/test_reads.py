import gzip
import tempfile
import unittest
from pathlib import Path

from anvaya.reads import Read, load_reads


class LoadReadsTests(unittest.TestCase):
    def test_loads_multiline_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fasta"
            path.write_text(">read-1\nATG\nCA\n>read-2\nGGTA\n", encoding="utf-8")

            self.assertEqual(
                load_reads(path),
                [Read("read-1", "ATGCA"), Read("read-2", "GGTA")],
            )

    def test_loads_fastq_with_phred_scores(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq"
            path.write_text("@read-1\nATGC\n+\nI!5?\n", encoding="utf-8")

            self.assertEqual(load_reads(path), [Read("read-1", "ATGC", (40, 0, 20, 30))])

    def test_loads_gzipped_fastq(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq.gz"
            with gzip.open(path, mode="wt", encoding="utf-8") as handle:
                handle.write("@read-1\nATGC\n+\nIIII\n")

            self.assertEqual(load_reads(path), [Read("read-1", "ATGC", (40, 40, 40, 40))])

    def test_rejects_mismatched_fastq_lengths(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq"
            path.write_text("@read-1\nATGC\n+\nIII\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                load_reads(path)


if __name__ == "__main__":
    unittest.main()
