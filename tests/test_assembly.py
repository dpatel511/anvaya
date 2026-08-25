import unittest

from anvaya.assembly import assemble


class AssembleTests(unittest.TestCase):
    def test_assembles_linear_reads(self) -> None:
        self.assertEqual(assemble(["ATGCA"], 3), ["ATGCA"])

    def test_assembles_branching_reads_into_unitigs(self) -> None:
        self.assertEqual(
            assemble(["ATGCA", "ATGGA"], 3),
            ["ATG", "TGCA", "TGGA"],
        )

    def test_preserves_graph_validation(self) -> None:
        with self.assertRaises(ValueError):
            assemble([], 3)


if __name__ == "__main__":
    unittest.main()
