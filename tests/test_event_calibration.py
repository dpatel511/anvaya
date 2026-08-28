import csv
import tempfile
import unittest
from pathlib import Path

from anvaya.event_calibration import (
    CALIBRATION_FIELDS,
    CalibrationExample,
    EventCalibrationScores,
    calibrate_event_report,
    calibrate_score_batch,
    fit_conformal_calibration,
    load_calibration_model,
    write_calibration_model,
)


def _scores(label: str, observations: int = 25) -> EventCalibrationScores:
    values = {"damage": -20.0, "sequencing_error": -20.0, "variation": -20.0}
    values[label] = -5.0
    return EventCalibrationScores(
        observations,
        values["damage"] - 0.1,
        values["damage"] + 0.1,
        values["sequencing_error"],
        values["variation"],
        5,
        observations - 5,
    )


def _model():
    examples = [
        CalibrationExample(label, _scores(label), f"{label}-{index}")
        for label in ("damage", "sequencing_error", "variation")
        for index in range(30)
    ]
    return fit_conformal_calibration(examples, minimum_per_class=20)


class EventCalibrationTests(unittest.TestCase):
    def test_clear_error_is_eligible_after_protective_null_rejection(self) -> None:
        result = calibrate_score_batch(
            [_scores("sequencing_error")],
            _model(),
            alpha=0.05,
        )[0]

        self.assertEqual(result.status, "calibrated")
        self.assertEqual(result.decision, "eligible_error")
        self.assertGreater(result.error_pvalue or 0.0, 0.05)
        self.assertLessEqual(result.damage_qvalue or 1.0, 0.05)
        self.assertLessEqual(result.variation_qvalue or 1.0, 0.05)

    def test_ambiguous_event_is_not_eligible_for_cleaning(self) -> None:
        ambiguous = EventCalibrationScores(25, -5.1, -4.9, -5.0, -20.0, 5, 20)
        result = calibrate_score_batch(
            [ambiguous],
            _model(),
            alpha=0.05,
        )[0]

        self.assertNotEqual(result.decision, "eligible_error")
        self.assertEqual(result.decision, "insufficient")
        self.assertIn("ambiguous_conformal_prediction_set", result.reasons)

    def test_profile_instability_blocks_a_decision(self) -> None:
        unstable = EventCalibrationScores(25, -30.0, 0.0, -5.0, -20.0, 5, 20)
        result = calibrate_score_batch(
            [unstable],
            _model(),
            alpha=0.05,
            maximum_damage_instability=1.0,
        )[0]

        self.assertEqual(result.decision, "insufficient")
        self.assertIn("damage_profile_unstable", result.reasons)

    def test_singleton_alternative_cannot_be_eligible_for_cleaning(self) -> None:
        singleton = EventCalibrationScores(100, -30.1, -29.9, -5.0, -20.0, 1, 99)
        result = calibrate_score_batch(
            [singleton],
            _model(),
            alpha=0.05,
        )[0]

        self.assertEqual(result.decision, "insufficient")
        self.assertIn("too_few_alternative_observations", result.reasons)

    def test_model_round_trip_and_report_calibration(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model_path = root / "model.json"
            input_path = root / "events.tsv"
            output_path = root / "calibrated.tsv"
            write_calibration_model(_model(), model_path)
            loaded = load_calibration_model(model_path)
            self.assertEqual(loaded.schema_version, 1)

            fields = (
                "event_id",
                "likelihood_status",
                "likelihood_observations",
                "likelihood_alternative_observations",
                "likelihood_reference_observations",
                "damage_model_log_likelihood",
                "damage_model_log_likelihood_min",
                "damage_model_log_likelihood_max",
                "error_model_log_likelihood",
                "variation_model_penalized_log_likelihood",
            )
            with input_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "event_id": "event-1",
                        "likelihood_status": "scored",
                        "likelihood_observations": 25,
                        "likelihood_alternative_observations": 5,
                        "likelihood_reference_observations": 20,
                        "damage_model_log_likelihood": -20.0,
                        "damage_model_log_likelihood_min": -20.1,
                        "damage_model_log_likelihood_max": -19.9,
                        "error_model_log_likelihood": -5.0,
                        "variation_model_penalized_log_likelihood": -20.0,
                    }
                )
            summary = calibrate_event_report(
                input_path,
                model_path,
                output_path,
                alpha=0.05,
            )
            self.assertEqual(summary.eligible_error, 1)
            with output_path.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(set(CALIBRATION_FIELDS).issubset(row))
            self.assertEqual(row["calibration_decision"], "eligible_error")

    def test_rejects_alpha_below_empirical_pvalue_resolution(self) -> None:
        with self.assertRaisesRegex(ValueError, "lacks p-value resolution"):
            calibrate_score_batch(
                [_scores("sequencing_error")],
                _model(),
                alpha=0.01,
            )


if __name__ == "__main__":
    unittest.main()
