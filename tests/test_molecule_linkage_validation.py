import importlib.util
import sys
import unittest
from pathlib import Path


EXPERIMENT = (
    Path(__file__).parents[1]
    / "experiments"
    / "08_molecule_linkage_validation.py"
)
SPEC = importlib.util.spec_from_file_location(
    "molecule_linkage_validation",
    EXPERIMENT,
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class MoleculeLinkageValidationTests(unittest.TestCase):
    def test_exact_mcnemar_pvalue(self) -> None:
        self.assertEqual(MODULE.exact_mcnemar_pvalue(0, 0), 1.0)
        self.assertEqual(MODULE.exact_mcnemar_pvalue(1, 0), 1.0)
        self.assertAlmostEqual(
            MODULE.exact_mcnemar_pvalue(0, 6),
            0.03125,
        )

    def test_paired_terminal_errors_change_both_ends(self) -> None:
        read = MODULE.previous.SimulatedRead(
            name="read",
            sequence="CAAAAAAAG",
            source_sequence="CAAAAAAAG",
            reference_start=0,
            reference_end=9,
        )

        changed = MODULE.add_paired_terminal_errors([read], 1.0, 1)[0]

        self.assertEqual(changed.sequence, "TAAAAAAAA")
        self.assertEqual(changed.error_positions, (0, 8))


if __name__ == "__main__":
    unittest.main()
