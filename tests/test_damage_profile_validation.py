import importlib.util
import sys
import unittest
from pathlib import Path


EXPERIMENT = (
    Path(__file__).parents[1]
    / "experiments"
    / "09_damage_profile_validation.py"
)
SPEC = importlib.util.spec_from_file_location(
    "damage_profile_validation",
    EXPERIMENT,
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class DamageProfileValidationTests(unittest.TestCase):
    def test_pearson_correlation(self) -> None:
        self.assertAlmostEqual(MODULE.pearson([3, 2, 1], [6, 4, 2]), 1.0)
        self.assertAlmostEqual(MODULE.pearson([3, 2, 1], [2, 4, 6]), -1.0)
        self.assertIsNone(MODULE.pearson([1, 1], [1, 2]))


if __name__ == "__main__":
    unittest.main()
