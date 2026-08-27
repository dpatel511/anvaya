import json
import tempfile
import unittest
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.damage_profile import (
    DamageProfileLocus,
    _summarize_loci,
    infer_damage_profile,
    write_damage_profile,
)
from anvaya.tip_matching import match_tip_to_backbone


class DamageProfileTests(unittest.TestCase):
    def _match(self, backbone: str, alternative: str, end_window: int = 3):
        graph = build_bidirected_dbg(
            [backbone] * 10 + [alternative],
            5,
            end_window=end_window,
            track_molecule_links=True,
        )
        tip = find_weak_tip_candidates(graph)[0]
        match = match_tip_to_backbone(graph, tip)
        self.assertIsNotNone(match)
        assert match is not None
        return graph, match

    def test_infers_five_prime_candidate_locus_fraction(self) -> None:
        graph, match = self._match(
            "AAGCCCAAA", "AAGCCTAAA", end_window=6
        )

        profile = infer_damage_profile(graph, [match])

        self.assertEqual(profile.matched_paths, 1)
        self.assertEqual(profile.matched_loci, 1)
        self.assertEqual(profile.observed_loci, 1)
        self.assertEqual(profile.eligible_loci, 1)
        self.assertEqual(
            sum(value.five_prime_ct_alternative for value in profile.bins),
            1,
        )
        self.assertEqual(
            sum(value.five_prime_c_reference for value in profile.bins),
            10,
        )
        self.assertEqual(profile.bins[5].five_prime_ct_fraction, 1 / 11)
        self.assertEqual(profile.bins[5].damage_fraction, 1 / 11)
        self.assertEqual(profile.bins[5].alternative_edge_contributors, 1)
        self.assertEqual(profile.bins[5].reference_edge_contributors, 1)
        self.assertEqual(
            profile.bins[5].largest_alternative_edge_contribution, 1
        )
        self.assertEqual(
            profile.bins[5].largest_reference_edge_contribution, 10
        )
        self.assertEqual(profile.bins[5].equal_locus_contributors, 1)
        self.assertEqual(profile.bins[5].equal_locus_fraction, 1 / 11)
        self.assertEqual(profile.bins[5].coverage_capped_weight, 11)
        self.assertEqual(profile.bins[5].coverage_capped_fraction, 1 / 11)
        self.assertEqual(profile.coverage_cap, 20)
        self.assertEqual(profile.loci[0].channel, "five_prime_ct")
        self.assertEqual(profile.loci[0].alternative[5], 1)
        self.assertEqual(profile.loci[0].reference[5], 10)
        self.assertTrue(
            all(value.three_prime_ga_fraction is None for value in profile.bins)
        )

    def test_normalizes_three_prime_profile(self) -> None:
        graph, match = self._match("ACATACACG", "ACATACACA")

        profile = infer_damage_profile(graph, [match])

        self.assertEqual(profile.bins[0].three_prime_ga_alternative, 1)
        self.assertEqual(profile.bins[0].three_prime_g_reference, 10)
        self.assertEqual(profile.bins[0].three_prime_ga_fraction, 1 / 11)
        self.assertEqual(profile.bins[0].damage_fraction, 1 / 11)

    def test_counts_duplicate_graph_locus_once(self) -> None:
        graph, match = self._match(
            "AAGCCCAAA", "AAGCCTAAA", end_window=6
        )

        profile = infer_damage_profile(graph, [match, match])

        self.assertEqual(profile.matched_paths, 2)
        self.assertEqual(profile.eligible_loci, 1)
        self.assertEqual(
            sum(value.five_prime_ct_alternative for value in profile.bins),
            1,
        )

    def test_excludes_kmer_boundary_without_terminal_changed_base(self) -> None:
        graph, match = self._match("AAGCCCAAA", "AAGCCTAAA")

        profile = infer_damage_profile(graph, [match])

        self.assertEqual(profile.matched_loci, 1)
        self.assertEqual(profile.observed_loci, 0)
        self.assertEqual(profile.eligible_loci, 0)
        self.assertTrue(
            all(value.damage_fraction is None for value in profile.bins)
        )

    def test_requires_molecule_links(self) -> None:
        graph = build_bidirected_dbg(["AAGCCCAAA"], 5, end_window=3)

        with self.assertRaisesRegex(ValueError, "requires molecule links"):
            infer_damage_profile(graph, [])

    def test_rejects_invalid_coverage_cap(self) -> None:
        graph, match = self._match(
            "AAGCCCAAA", "AAGCCTAAA", end_window=6
        )

        with self.assertRaisesRegex(ValueError, "at least 1"):
            infer_damage_profile(graph, [match], coverage_cap=0)
        with self.assertRaisesRegex(TypeError, "must be an integer"):
            infer_damage_profile(graph, [match], coverage_cap=1.5)  # type: ignore[arg-type]

    def test_caps_influence_of_deep_locus(self) -> None:
        loci = (
            DamageProfileLocus(1, 2, "five_prime_ct", (1,), (99,)),
            DamageProfileLocus(3, 4, "five_prime_ct", (1,), (1,)),
        )

        contributors, equal, weight, capped = _summarize_loci(loci, 0, 20)

        self.assertEqual(contributors, 2)
        self.assertEqual(equal, 0.255)
        self.assertEqual(weight, 22)
        self.assertAlmostEqual(capped, 12 / 220)

    def test_writes_explicit_profile_schema(self) -> None:
        graph, match = self._match(
            "AAGCCCAAA", "AAGCCTAAA", end_window=6
        )
        profile = infer_damage_profile(graph, [match])

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "profile.json"
            write_damage_profile(profile, path)
            payload = json.loads(path.read_text(encoding="utf-8"))

        self.assertEqual(payload["schema_version"], 4)
        self.assertEqual(
            payload["interpretation"],
            "matched_candidate_locus_fraction",
        )
        self.assertEqual(payload["coverage_cap"], 20)
        self.assertEqual(payload["matched_loci"], 1)
        self.assertEqual(payload["observed_loci"], 1)
        self.assertEqual(len(payload["loci"]), 1)
        self.assertEqual(len(payload["bins"]), graph.end_window)


if __name__ == "__main__":
    unittest.main()
