import tempfile
import unittest
from pathlib import Path

from anvaya.output import write_fasta


class WriteFastaTests(unittest.TestCase):
    def test_writes_named_wrapped_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "unitigs.fasta"

            write_fasta(["ACGTAC", "TT"], path, line_width=4)

            self.assertEqual(
                path.read_text(encoding="utf-8"),
                ">unitig_1 length=6\nACGT\nAC\n>unitig_2 length=2\nTT\n",
            )

    def test_rejects_non_positive_line_width(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "unitigs.fasta"

            with self.assertRaises(ValueError):
                write_fasta(["ACGT"], path, line_width=0)


if __name__ == "__main__":
    unittest.main()
