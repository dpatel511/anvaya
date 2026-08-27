import unittest

from anvaya.damage_likelihood import fit_candidate_damage_model


def _replicated_counts(
    alternatives: tuple[int, ...],
    references: tuple[int, ...],
    loci: int = 20,
):
    return tuple((alternatives, references) for _ in range(loci))


class CandidateDamageLikelihoodTests(unittest.TestCase):
    def test_fits_terminal_decay_against_constant_background(self) -> None:
        fit = fit_candidate_damage_model(
            _replicated_counts(
                (8, 5, 3, 2, 1),
                (22, 25, 27, 28, 29),
            ),
            end_window=5,
        )

        self.assertEqual(fit.status, "fitted")
        self.assertEqual(fit.observation_scope, "matched_candidate_loci")
        self.assertIsNotNone(fit.amplitude)
        self.assertIsNotNone(fit.decay)
        self.assertIsNotNone(fit.background)
        assert fit.amplitude is not None
        assert fit.decay is not None
        self.assertGreater(fit.amplitude.value, 0.1)
        self.assertLess(fit.decay.value, 0.9)
        self.assertLessEqual(
            fit.amplitude.lower_95_support,
            fit.amplitude.value,
        )
        self.assertGreaterEqual(
            fit.amplitude.upper_95_support,
            fit.amplitude.value,
        )
        self.assertGreater(fit.likelihood_ratio_statistic or 0.0, 10.0)
        self.assertTrue(
            all(
                left > right
                for left, right in zip(
                    fit.fitted_probabilities,
                    fit.fitted_probabilities[1:],
                )
            )
        )

    def test_constant_profile_does_not_gain_decay_evidence(self) -> None:
        fit = fit_candidate_damage_model(
            _replicated_counts(
                (3, 3, 3, 3, 3),
                (27, 27, 27, 27, 27),
            ),
            end_window=5,
        )

        self.assertEqual(fit.status, "fitted")
        self.assertLess(fit.likelihood_ratio_statistic or 0.0, 0.1)

    def test_refuses_fit_when_evidence_is_sparse(self) -> None:
        fit = fit_candidate_damage_model(
            (((1, 0, 0), (4, 0, 0)),),
            end_window=3,
        )

        self.assertEqual(fit.status, "insufficient_evidence")
        self.assertIn("fewer_than_10_observed_loci", fit.reasons)
        self.assertIn("fewer_than_3_observed_distances", fit.reasons)
        self.assertIsNone(fit.amplitude)
        self.assertEqual(fit.fitted_probabilities, ())

    def test_rejects_malformed_locus_counts(self) -> None:
        with self.assertRaisesRegex(ValueError, "match end_window"):
            fit_candidate_damage_model((((1,), (2,)),), end_window=2)
        with self.assertRaisesRegex(ValueError, "must not be negative"):
            fit_candidate_damage_model((((-1,), (2,)),), end_window=1)
        with self.assertRaisesRegex(TypeError, "must be integers"):
            fit_candidate_damage_model((((1.5,), (2,)),), end_window=1)


if __name__ == "__main__":
    unittest.main()
