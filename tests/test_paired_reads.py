import unittest

from anvaya.paired_reads import merge_overlapping_pairs
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


class PairedReadMergeTests(unittest.TestCase):
    def test_merges_a_unique_pair_overlap(self) -> None:
        fragment = "AACCGGTTAACCGGTTACGATGCA"
        left = Read("pair/1", fragment[:16], tuple(range(16)))
        right_sequence = reverse_complement(fragment[8:])
        right = Read("pair/2", right_sequence, tuple(range(16, 32)))

        reads, molecules, diagnostics = merge_overlapping_pairs(
            [left], [right], minimum_overlap=8,
            minimum_identity=1.0, minimum_ry_identity=1.0,
        )

        self.assertEqual([read.sequence for read in reads], [fragment])
        self.assertEqual(molecules, [0])
        self.assertEqual(diagnostics.merged_pairs, 1)
        self.assertEqual(len(reads[0].qualities), len(fragment))

    def test_retains_both_mates_without_supported_overlap(self) -> None:
        left = Read("pair/1", "AACCGGTTAACC")
        right = Read("pair/2", "TTTTGGGGTTTT")

        reads, molecules, diagnostics = merge_overlapping_pairs(
            [left], [right], minimum_overlap=8,
            minimum_identity=1.0, minimum_ry_identity=1.0,
        )

        self.assertEqual(reads, [left, right])
        self.assertEqual(molecules, [0, 0])
        self.assertEqual(diagnostics.unmerged_pairs, 1)


if __name__ == "__main__":
    unittest.main()
