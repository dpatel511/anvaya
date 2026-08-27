import json
import tempfile
import unittest
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.damage_profile import infer_damage_profile, write_damage_profile
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
        graph, match = self._match("AAGCCCAAA", "AAGCCTAAA")

        profile = infer_damage_profile(graph, [match])

        self.assertEqual(profile.matched_paths, 1)
        self.assertEqual(profile.eligible_loci, 1)
        self.assertEqual(
            sum(value.five_prime_ct_alternative for value in profile.bins),
            1,
        )
        self.assertEqual(
            sum(value.five_prime_c_reference for value in profile.bins),
            10,
        )
        self.assertEqual(profile.bins[1].five_prime_ct_fraction, 1 / 11)
        self.assertEqual(profile.bins[1].damage_fraction, 1 / 11)
        self.assertEqual(profile.bins[1].alternative_edge_contributors, 1)
        self.assertEqual(profile.bins[1].reference_edge_contributors, 1)
        self.assertEqual(
            profile.bins[1].largest_alternative_edge_contribution, 1
        )
        self.assertEqual(
            profile.bins[1].largest_reference_edge_contribution, 10
        )
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
        graph, match = self._match("AAGCCCAAA", "AAGCCTAAA")

        profile = infer_damage_profile(graph, [match, match])

        self.assertEqual(profile.matched_paths, 2)
        self.assertEqual(profile.eligible_loci, 1)
        self.assertEqual(
            sum(value.five_prime_ct_alternative for value in profile.bins),
            1,
        )

    def test_requires_molecule_links(self) -> None:
        graph = build_bidirected_dbg(["AAGCCCAAA"], 5, end_window=3)

        with self.assertRaisesRegex(ValueError, "requires molecule links"):
            infer_damage_profile(graph, [])

    def test_writes_explicit_profile_schema(self) -> None:
        graph, match = self._match("AAGCCCAAA", "AAGCCTAAA")
        profile = infer_damage_profile(graph, [match])

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "profile.json"
            write_damage_profile(profile, path)
            payload = json.loads(path.read_text(encoding="utf-8"))

        self.assertEqual(payload["schema_version"], 2)
        self.assertEqual(
            payload["interpretation"],
            "matched_candidate_locus_fraction",
        )
        self.assertEqual(len(payload["bins"]), graph.end_window)


if __name__ == "__main__":
    unittest.main()
