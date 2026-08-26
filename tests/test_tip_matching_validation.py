import importlib.util
import sys
import unittest
from pathlib import Path

from helpers.simulation import SimulatedRead


def _load_experiment():
    path = (
        Path(__file__).resolve().parents[1]
        / "experiments"
        / "06_tip_matching_validation.py"
    )
    specification = importlib.util.spec_from_file_location(
        "tip_matching_validation",
        path,
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("could not load validation experiment")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


experiment = _load_experiment()


class TipMatchingValidationTests(unittest.TestCase):
    def test_terminal_errors_do_not_change_internal_bases(self) -> None:
        read = SimulatedRead(
            name="read-1",
            sequence="ACGT" * 15,
            source_sequence="ACGT" * 15,
            reference_start=0,
            reference_end=60,
        )

        changed = experiment.add_terminal_errors([read], 1.0, 7)[0]

        self.assertEqual(
            changed.sequence[
                experiment.END_WINDOW : -experiment.END_WINDOW
            ],
            read.sequence[
                experiment.END_WINDOW : -experiment.END_WINDOW
            ],
        )
        self.assertEqual(len(changed.error_positions), 2 * experiment.END_WINDOW)
        self.assertTrue(
            all(
                position < experiment.END_WINDOW
                or position >= len(read.sequence) - experiment.END_WINDOW
                for position in changed.error_positions
            )
        )

    def test_related_reference_has_reproducible_spaced_snps(self) -> None:
        reference = "A" * experiment.REFERENCE_LENGTH

        first, first_positions = experiment.rare_reference(reference, 11)
        second, second_positions = experiment.rare_reference(reference, 11)

        self.assertEqual(first, second)
        self.assertEqual(first_positions, second_positions)
        self.assertEqual(len(first_positions), experiment.SNP_COUNT)
        self.assertEqual(
            sum(
                left != right
                for left, right in zip(reference, first, strict=True)
            ),
            experiment.SNP_COUNT,
        )

    def test_strain_label_distinguishes_rare_and_shared_sequence(self) -> None:
        major = "AAAACCCCGGGGTTTT"
        rare = "AAAATCCCGGGGTTTT"

        self.assertEqual(
            experiment.strain_label("AAAATCCC", major, rare),
            "rare",
        )
        self.assertEqual(
            experiment.strain_label("GGGGTTTT", major, rare),
            "shared",
        )


if __name__ == "__main__":
    unittest.main()
