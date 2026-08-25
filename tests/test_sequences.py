import unittest

from anvaya.sequences import canonical_sequence, normalize_dna, reverse_complement


class SequenceUtilitiesTests(unittest.TestCase):
    def test_normalizes_dna(self) -> None:
        self.assertEqual(normalize_dna("acgtn"), "ACGTN")

    def test_reverse_complement(self) -> None:
        self.assertEqual(reverse_complement("ACGTN"), "NACGT")

    def test_canonical_sequence(self) -> None:
        self.assertEqual(canonical_sequence("TGC"), "GCA")

    def test_rejects_invalid_symbol(self) -> None:
        with self.assertRaises(ValueError):
            normalize_dna("ACGX")


if __name__ == "__main__":
    unittest.main()
