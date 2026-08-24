import unittest

from anvaya.kmers import kmers


class KmersTests(unittest.TestCase):
    def test_returns_overlapping_kmers(self) -> None:
        self.assertEqual(kmers("ATGCA", 3), ["ATG", "TGC", "GCA"])

    def test_normalizes_lowercase_sequence(self) -> None:
        self.assertEqual(kmers("atgn", 2), ["AT", "TG", "GN"])

    def test_accepts_k_equal_to_sequence_length(self) -> None:
        self.assertEqual(kmers("ACGT", 4), ["ACGT"])

    def test_rejects_empty_sequence(self) -> None:
        with self.assertRaises(ValueError):
            kmers("", 1)

    def test_rejects_non_positive_k(self) -> None:
        with self.assertRaises(ValueError):
            kmers("ACGT", 0)

    def test_rejects_k_larger_than_sequence(self) -> None:
        with self.assertRaises(ValueError):
            kmers("ACGT", 5)

    def test_rejects_non_integer_k(self) -> None:
        with self.assertRaises(TypeError):
            kmers("ACGT", 2.5)  # type: ignore[arg-type]

    def test_rejects_unsupported_symbols(self) -> None:
        with self.assertRaises(ValueError):
            kmers("ACGX", 2)


if __name__ == "__main__":
    unittest.main()
