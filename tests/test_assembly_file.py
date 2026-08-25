import tempfile
import unittest
from pathlib import Path

from anvaya.assembly import assemble_file


class AssembleFileTests(unittest.TestCase):
    def test_assembles_reads_from_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fasta"
            path.write_text(">read-1\nATGCA\n", encoding="utf-8")

            self.assertEqual(assemble_file(path, 3), ["ATGCA"])


if __name__ == "__main__":
    unittest.main()
