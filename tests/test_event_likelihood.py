import unittest

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.damage_likelihood import fit_candidate_damage_model
from anvaya.damage_profile import DamageProfile, DamageProfileLocus
from anvaya.event_likelihood import (
    CrossFittedDamageModels,
    NucleotideObservation,
    compare_event_likelihoods,
    fit_cross_fitted_damage_models,
    score_matched_event,
)
from anvaya.tip_matching import match_tip_to_backbone


def _observations(
    alternative: int,
    reference: int,
    distance: int,
    alternative_quality: int = 35,
) -> list[NucleotideObservation]:
    return [
        NucleotideObservation(True, distance, alternative_quality)
        for _ in range(alternative)
    ] + [NucleotideObservation(False, distance, 35) for _ in range(reference)]


class EventLikelihoodTests(unittest.TestCase):
    def test_reports_distinct_molecules_across_repeated_base_evidence(self) -> None:
        observations = [
            NucleotideObservation(True, 0, 30, molecule_id=1),
            NucleotideObservation(True, 1, 30, molecule_id=1),
            NucleotideObservation(False, 0, 30, molecule_id=2),
            NucleotideObservation(False, 1, 30, molecule_id=2),
        ]

        result = compare_event_likelihoods(observations, (0.2, 0.1))

        self.assertEqual(result.observations, 2)
        self.assertEqual(result.alternative_observations, 1)
        self.assertEqual(result.reference_observations, 1)

    def test_excludes_molecule_with_conflicting_allele_evidence(self) -> None:
        observations = [
            NucleotideObservation(True, 0, 30, molecule_id=1),
            NucleotideObservation(False, 1, 30, molecule_id=1),
            NucleotideObservation(True, 0, 30, molecule_id=2),
            NucleotideObservation(False, 0, 30, molecule_id=3),
        ]

        result = compare_event_likelihoods(observations, (0.2, 0.1))

        self.assertEqual(result.observations, 2)
        self.assertEqual(result.alternative_observations, 1)
        self.assertEqual(result.reference_observations, 1)

    def test_ordinary_substitution_cannot_rank_as_damage_on_tie(self) -> None:
        result = compare_event_likelihoods(
            _observations(2, 8, 0),
            (0.0, 0.0),
            damage_explanation_applicable=False,
        )

        self.assertNotEqual(result.best_explanation, "damage")

    def test_damage_curve_beats_penalized_fitted_variation(self) -> None:
        result = compare_event_likelihoods(
            _observations(4, 16, 0),
            (0.2,),
        )

        self.assertEqual(result.status, "scored")
        self.assertEqual(result.best_explanation, "damage")
        self.assertAlmostEqual(result.variation_frequency or 0.0, 0.2, places=3)
        self.assertGreater(result.log_likelihood_margin or 0.0, 0.0)

    def test_low_quality_single_alternative_favors_error(self) -> None:
        result = compare_event_likelihoods(
            _observations(1, 99, 0, alternative_quality=8),
            (0.2,),
        )

        self.assertEqual(result.best_explanation, "sequencing_error")

    def test_position_independent_support_favors_variation(self) -> None:
        observations = []
        for distance in range(5):
            observations.extend(_observations(4, 16, distance))
        result = compare_event_likelihoods(
            observations,
            (0.35, 0.20, 0.10, 0.04, 0.01),
        )

        self.assertEqual(result.best_explanation, "variation")
        self.assertAlmostEqual(result.variation_frequency or 0.0, 0.2, places=3)

    def test_cross_fit_ensemble_exposes_damage_score_instability(self) -> None:
        observations = [
            NucleotideObservation(
                observation.alternative,
                observation.distance,
                observation.quality,
                damage_probability=0.2,
                damage_probability_ensemble=(0.05, 0.35),
            )
            for observation in _observations(4, 16, 0)
        ]
        result = compare_event_likelihoods(observations, ())

        self.assertEqual(result.status, "scored")
        self.assertIsNotNone(result.damage_log_likelihood_min)
        self.assertIsNotNone(result.damage_log_likelihood_max)
        self.assertGreater(result.damage_log_likelihood_spread or 0.0, 0.0)
        self.assertAlmostEqual(
            result.damage_error_log_contrast or 0.0,
            (result.damage_log_likelihood or 0.0)
            - (result.error_log_likelihood or 0.0),
        )

    def test_missing_quality_keeps_event_unscored(self) -> None:
        result = compare_event_likelihoods(
            _observations(2, 10, 0),
            (0.2,),
            missing_quality_observations=1,
        )

        self.assertEqual(result.status, "insufficient_evidence")
        self.assertIn("missing_base_quality", result.reasons)
        self.assertIsNone(result.best_explanation)

    def test_cross_fitting_builds_models_without_held_out_fold(self) -> None:
        loci = tuple(
            DamageProfileLocus(
                alternative_edge_id=index,
                reference_edge_id=100 + index,
                channel="five_prime_ct",
                alternative=(3, 2, 1, 0, 0),
                reference=(7, 8, 9, 10, 10),
            )
            for index in range(20)
        )
        full_fit = fit_candidate_damage_model(
            tuple((locus.alternative, locus.reference) for locus in loci),
            5,
        )
        profile = DamageProfile(
            end_window=5,
            matched_paths=20,
            matched_loci=20,
            observed_loci=20,
            coverage_cap=20,
            bins=(),
            loci=loci,
            candidate_damage_fit=full_fit,
        )

        models = fit_cross_fitted_damage_models(profile, folds=5)

        self.assertEqual(models.status, "fitted")
        self.assertEqual(models.fitted_folds, 5)
        self.assertEqual(len(models.locus_folds), 20)
        self.assertEqual(set(models.probabilities_by_fold), set(range(5)))

    def test_scores_exact_tip_and_backbone_molecules(self) -> None:
        backbone = "AAGCCCAAA"
        alternative = "AAGCCTAAA"
        reads = [backbone] * 30 + [alternative] * 3
        graph = build_bidirected_dbg(
            reads,
            5,
            end_window=11,
            track_molecule_links=True,
            read_qualities=[(35,) * 9] * len(reads),
        )
        tip = find_weak_tip_candidates(graph)[0]
        match = match_tip_to_backbone(graph, tip)
        self.assertIsNotNone(match)
        assert match is not None
        change = match.substitutions[0]
        edge_index = change.position - graph.node_length - 1
        key = (
            match.tip_edge_ids[edge_index],
            match.backbone_edge_ids[edge_index],
            "five_prime_ct",
        )
        models = CrossFittedDamageModels(
            status="fitted",
            reasons=(),
            folds=5,
            fitted_folds=5,
            locus_folds={key: 0},
            probabilities_by_fold={0: (0.1,) * 11},
            independent_probabilities=(0.1,) * 11,
        )

        result = score_matched_event(graph, match, models)

        self.assertEqual(result.status, "scored")
        self.assertEqual(result.alternative_observations, 3)
        self.assertEqual(result.reference_observations, 30)
        self.assertEqual(result.best_explanation, "damage")
        self.assertEqual(
            result.profile_scope,
            "held_out_candidate_fold_1_of_5",
        )


if __name__ == "__main__":
    unittest.main()
