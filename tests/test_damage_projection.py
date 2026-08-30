import unittest

from anvaya.damage_projection import decide_damage_projection
from anvaya.event_likelihood import EventLikelihood
from anvaya.tip_classification import TipClassification
from anvaya.tip_matching import TipBackboneMatch


def _match() -> TipBackboneMatch:
    return TipBackboneMatch(
        branch_handle=0,
        tip_edge_ids=(1,),
        backbone_edge_ids=(2,),
        tip_sequence="AAT",
        backbone_sequence="AAC",
        tip_observations=3,
        backbone_observations=30,
        backbone_minimum_edge_support=30,
        relative_coverage=0.1,
        sequence_identity=2 / 3,
        ry_identity=1.0,
        substitutions=(),
        damage_compatible=True,
    )


def _classification(label: str) -> TipClassification:
    return TipClassification(label, (), 0.9, 0.1, 0.1, None, ())  # type: ignore[arg-type]


def _likelihood(**overrides: object) -> EventLikelihood:
    values: dict[str, object] = {
        "status": "scored",
        "reasons": (),
        "profile_scope": "held_out_candidate_fold_1_of_5",
        "conditioning_scope": "two_path_fragment_ascertainment",
        "observations": 10,
        "alternative_observations": 2,
        "reference_observations": 8,
        "missing_quality_observations": 0,
        "damage_log_likelihood": -4.0,
        "damage_log_likelihood_min": -4.1,
        "damage_log_likelihood_max": -3.9,
        "damage_log_likelihood_spread": 0.2,
        "damage_conditioning_log_probability": -1.0,
        "error_log_likelihood": -8.0,
        "error_conditioning_log_probability": -1.0,
        "variation_log_likelihood": -5.0,
        "variation_conditioning_log_probability": -1.0,
        "variation_frequency": 0.2,
        "variation_complexity_penalty": 0.1,
        "variation_penalized_log_likelihood": -5.0,
        "best_explanation": "damage",
        "log_likelihood_margin": 1.0,
        "damage_error_log_contrast": 4.0,
        "damage_variation_log_contrast": 1.0,
        "error_variation_log_contrast": -3.0,
    }
    values.update(overrides)
    return EventLikelihood(**values)  # type: ignore[arg-type]


class DamageProjectionTests(unittest.TestCase):
    def test_would_project_only_when_all_gates_pass(self) -> None:
        decision = decide_damage_projection(
            _match(), _classification("damage-like"), _likelihood()
        )

        self.assertEqual(decision.decision, "would_project")
        self.assertEqual(decision.target_sequence, "AAC")

    def test_keeps_variation_competitor_unresolved(self) -> None:
        decision = decide_damage_projection(
            _match(),
            _classification("damage-like"),
            _likelihood(damage_variation_log_contrast=-0.1),
        )

        self.assertEqual(decision.decision, "leave_unresolved")
        self.assertEqual(decision.target_sequence, "")
        self.assertIn("damage_does_not_beat_variation", decision.reasons)

    def test_keeps_unscored_candidate_unresolved(self) -> None:
        decision = decide_damage_projection(
            _match(),
            _classification("damage-like"),
            _likelihood(status="insufficient_evidence", best_explanation=None),
        )

        self.assertEqual(decision.decision, "leave_unresolved")
        self.assertIn("likelihood_not_scored", decision.reasons)


if __name__ == "__main__":
    unittest.main()
