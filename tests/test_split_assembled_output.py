import json
import tempfile
import unittest
from pathlib import Path

from experiments.split_assembled_output import split_assembled_output


class SplitAssembledOutputTests(unittest.TestCase):
    def test_partitions_records_by_distinct_molecules_without_changing_sequences(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "all.fasta"
            placements = root / "placements.tsv"
            assembled = root / "assembled.fasta"
            unresolved = root / "unresolved.fasta"
            report = root / "report.json"
            fasta.write_text(
                ">unitig_1 description\nAAAA\nCC\n>unitig_2\nTTT\n",
                encoding="utf-8",
            )
            placements.write_text(
                "contig_id\tmolecule_id\tread_start_0based\tread_stop_exclusive\n"
                "unitig_1\t7\t0\t4\nunitig_1\t7\t1\t5\n"
                "unitig_1\t8\t0\t6\nunitig_2\t9\t0\t3\n",
                encoding="utf-8",
            )

            result = split_assembled_output(
                fasta, placements, assembled, unresolved, report,
            )

            self.assertEqual(assembled.read_text(encoding="utf-8"),
                             ">unitig_1 description\nAAAA\nCC\n")
            self.assertEqual(unresolved.read_text(encoding="utf-8"),
                             ">unitig_2\nTTT\n")
            self.assertEqual(result["assembled"]["records"], 1)
            self.assertEqual(result["assembled"]["bases"], 6)
            self.assertEqual(result["unresolved"]["records"], 1)
            self.assertEqual(result["unresolved"]["bases"], 3)
            self.assertEqual(json.loads(report.read_text(encoding="utf-8")), result)

    def test_rejects_noncontiguous_placement_groups(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "all.fasta"
            placements = root / "placements.tsv"
            fasta.write_text(">a\nAAAA\n>b\nCCCC\n", encoding="utf-8")
            placements.write_text(
                "contig_id\tmolecule_id\na\t1\nb\t2\na\t3\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "not grouped"):
                split_assembled_output(
                    fasta, placements, root / "a.fa", root / "u.fa",
                    root / "report.json",
                )


if __name__ == "__main__":
    unittest.main()
