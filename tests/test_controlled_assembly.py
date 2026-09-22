import importlib.util
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from anvaya.overlap_progressive import ProgressiveSequencePool, iterate_progressive_raw_extension
from anvaya.raw_consensus import RawPlacement
from anvaya.reads import Read
from anvaya.sequences import reverse_complement

path = Path(__file__).resolve().parents[1] / "experiments/controlled_assembly.py"
spec = importlib.util.spec_from_file_location("controlled_assembly", path)
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


class ControlledAssemblyTests(unittest.TestCase):
    def _ambiguous_pool(self, sequence, raw_count):
        reads = [Read(f"raw_{index}", sequence) for index in range(raw_count)]
        pool = ProgressiveSequencePool.from_reads(reads)
        records = [record.consumed() for record in pool.records]
        records[0] = replace(
            pool.records[0].corrected(Read("contig", sequence)),
            raw_placements=tuple(
                RawPlacement(index, 0, False, 0, len(sequence))
                for index in range(raw_count)
            ),
        )
        return ProgressiveSequencePool(tuple(records))

    def _placed_pool(self, sequences):
        reads = [Read(f"raw_{index}", sequence)
                 for index, sequence in enumerate(sequences)]
        pool = ProgressiveSequencePool.from_reads(reads)
        for index, read in enumerate(reads):
            pool = pool.replace_record(replace(
                pool.records[index].corrected(read),
                raw_placements=(RawPlacement(index, 0, False, 0, len(read.sequence)),),
            ))
        return pool

    def test_reuse_extends_from_consumed_reads_without_retiring_centers(self):
        reference = benchmark.fixture(93004, "unique", 3, 0)[0][0]
        reads = [Read("center", reference[:80])]
        reads += [Read(f"support{i}", reference[40:120]) for i in range(5)]
        pool = ProgressiveSequencePool.from_reads(reads)
        pool = pool.replace_record(replace(pool.records[0].corrected(reads[0]),
            raw_placements=(RawPlacement(0, 0, False, 0, 80),)))
        for index in range(1, 6):
            pool = pool.replace_record(pool.records[index].consumed())
        baseline, _ = iterate_progressive_raw_extension(pool)
        result, diagnostics = iterate_progressive_raw_extension(pool,
            minimum_anchor_matches=1, reuse_raw_evidence=True)
        self.assertEqual(baseline, pool)
        self.assertEqual(result.records[0].current.sequence, reference[:120])
        self.assertEqual(len(result.active_derived), 1)
        self.assertEqual(diagnostics.consumed_reads, 0)
        self.assertEqual(result.raw_evidence, pool.raw_evidence)
        origins = [(0, 0, 80, False)] + [(0, 40, 120, False)] * 5
        self.assertEqual(benchmark.score(result, [reference], origins)["mismatches"], 0)
        duplicate_pool = ProgressiveSequencePool(tuple(
            replace(record, molecule_id=1) if record.index else record for record in pool.records))
        duplicate_result, duplicate_diagnostics = iterate_progressive_raw_extension(
            duplicate_pool, minimum_anchor_matches=1, reuse_raw_evidence=True)
        self.assertEqual(duplicate_result.records[0].current.sequence, reads[0].sequence)
        self.assertEqual(duplicate_diagnostics.added_bases, 0)

    def test_clean_and_damaged_origins_match(self):
        clean = benchmark.fixture(93001, "strains", 3, 0)
        damaged = benchmark.fixture(93001, "strains", 3, .4)
        self.assertEqual(clean[0], damaged[0])
        self.assertEqual(clean[2], damaged[2])
        self.assertEqual(clean[3], [])
        self.assertGreater(len(damaged[3]), 0)
        self.assertEqual(clean[4], damaged[4])

    def test_sequencing_error_control_is_separate_from_damage(self):
        clean = benchmark.fixture(93005, "unique", 3, 0)
        erroneous = benchmark.fixture(93005, "unique", 3, 0, error_rate=1.0)
        self.assertEqual(clean[0], erroneous[0])
        self.assertEqual(clean[2], erroneous[2])
        self.assertEqual(erroneous[3], [])
        self.assertEqual(clean[4], [])
        self.assertGreater(len(erroneous[4]), 0)
        self.assertNotEqual(clean[1], erroneous[1])

    def test_reverse_origin_and_base_error_scoring(self):
        reference = "ACGTTGCA"
        read = Read("r", reverse_complement(reference[2:7]))
        pool = ProgressiveSequencePool.from_reads([read])
        record = replace(pool.records[0].corrected(read), raw_placements=(RawPlacement(0, 0, False, 0, 5),))
        pool = pool.replace_record(record)
        metrics = benchmark.score(pool, [reference], [(0, 2, 7, True)])
        self.assertEqual(metrics["resolved_unique_reference_bases"], 5)
        self.assertEqual(metrics["mismatches"], 0)
        altered = Read("r", "A" + read.sequence[1:])
        self.assertEqual(benchmark.score(pool, [reference], [(0, 2, 7, True)], [altered])["mismatches"], 1)

    def test_score_reports_repeat_masked_sequence_coverage(self):
        reference = "ACGTTGCA"
        pool = self._ambiguous_pool(reference, 1)
        metrics = benchmark.score(
            pool,
            [reference],
            [(0, 0, 8, False)],
            reference_ambiguous_intervals=[(0, 2, 6)],
        )
        self.assertEqual(metrics["sequence_resolved_unique_reference_bases"], 8)
        self.assertEqual(metrics["sequence_resolved_nonambiguous_reference_bases"], 4)
        self.assertEqual(metrics["sequence_resolved_ambiguous_reference_bases"], 4)

    def test_score_deduplicates_concordant_reference_coordinate_observations(self):
        correct = benchmark.score(
            self._placed_pool(["A", "A"]),
            ["A"],
            [(0, 0, 1, False), (0, 0, 1, False)],
        )
        self.assertEqual(correct["nonduplicated_accuracy_projected_bases"], 2)
        self.assertEqual(correct["nonduplicated_accuracy_duplicate_observations"], 1)
        self.assertEqual(correct["nonduplicated_accuracy_evaluated_reference_bases"], 1)
        self.assertEqual(correct["nonduplicated_accuracy_mismatches"], 0)

        erroneous = benchmark.score(
            self._placed_pool(["C", "C"]),
            ["A"],
            [(0, 0, 1, False), (0, 0, 1, False)],
        )
        self.assertEqual(erroneous["nonduplicated_accuracy_evaluated_reference_bases"], 1)
        self.assertEqual(erroneous["nonduplicated_accuracy_mismatches"], 1)

    def test_score_abstains_on_discordant_reference_coordinate_observations(self):
        metrics = benchmark.score(
            self._placed_pool(["A", "C"]),
            ["A"],
            [(0, 0, 1, False), (0, 0, 1, False)],
        )
        self.assertEqual(metrics["nonduplicated_accuracy_projected_bases"], 2)
        self.assertEqual(metrics["nonduplicated_accuracy_duplicate_observations"], 1)
        self.assertEqual(metrics["nonduplicated_accuracy_evaluated_reference_bases"], 0)
        self.assertEqual(metrics["nonduplicated_accuracy_mismatches"], 0)
        self.assertEqual(metrics["nonduplicated_accuracy_discordant_reference_bases"], 1)
        self.assertEqual(metrics["nonduplicated_accuracy_discordant_sites"], [{
            "reference": 0,
            "position_0based": 0,
            "observations": [
                {"base": "A", "contig_id": "unitig_1"},
                {"base": "C", "contig_id": "unitig_2"},
            ],
        }])

    def test_score_optionally_reports_coordinate_states_and_resolution(self):
        metrics = benchmark.score(
            self._placed_pool(["A", "C"]),
            ["AA"],
            [(0, 0, 1, False), (0, 1, 2, False)],
            include_coordinate_states=True,
        )
        self.assertEqual(metrics["nonduplicated_accuracy_coordinate_states"], [
            {
                "reference": 0,
                "position_0based": 0,
                "reference_base": "A",
                "observed_bases": ["A"],
                "observation_count": 1,
                "origin_resolution": ["unique_origin"],
                "status": "concordant_correct",
            },
            {
                "reference": 0,
                "position_0based": 1,
                "reference_base": "A",
                "observed_bases": ["C"],
                "observation_count": 1,
                "origin_resolution": ["unique_origin"],
                "status": "concordant_incorrect",
            },
        ])
        self.assertNotIn(
            "nonduplicated_accuracy_coordinate_states",
            benchmark.score(self._placed_pool(["A"]), ["A"], [(0, 0, 1, False)]),
        )

    def test_score_projects_reverse_base_to_zero_based_reference_coordinate(self):
        metrics = benchmark.score(
            self._placed_pool(["T"]),
            ["A"],
            [(0, 0, 1, True)],
        )
        self.assertEqual(metrics["nonduplicated_accuracy_evaluated_reference_bases"], 1)
        self.assertEqual(metrics["nonduplicated_accuracy_mismatches"], 0)

    def test_score_separates_repeat_and_cross_reference_provenance(self):
        repeat = benchmark.score(
            self._ambiguous_pool("ACGT", 2),
            ["ACGTGGGGACGT"],
            [(0, 0, 4, False), (0, 8, 12, False)],
        )
        evidence = repeat["contig_evidence"][0]
        self.assertEqual(evidence["provenance_category"], "same_reference_multiple_transforms")
        self.assertEqual(evidence["sequence_category"], "identical_candidate_sequences")
        self.assertEqual(repeat["accuracy_evaluated_bases"], 4)
        self.assertEqual(repeat["accuracy_mismatches"], 0)
        self.assertEqual(repeat["unresolved_bases"], 4)
        self.assertEqual(repeat["sequence_resolved_unique_reference_bases"], 0)

        shared = benchmark.score(
            self._ambiguous_pool("ACGT", 2),
            ["ACGT", "ACGT"],
            [(0, 0, 4, False), (1, 0, 4, False)],
        )
        evidence = shared["contig_evidence"][0]
        self.assertEqual(evidence["provenance_category"], "cross_reference_transforms")
        self.assertEqual(evidence["sequence_category"], "identical_candidate_sequences")

    def test_score_distinguishes_compatible_origin_and_diagnostic_mosaic(self):
        compatible = benchmark.score(
            self._ambiguous_pool("AAAA", 2),
            ["AAAA", "AAAT"],
            [(0, 0, 4, False), (1, 0, 4, False)],
        )
        self.assertEqual(
            compatible["contig_evidence"][0]["sequence_category"],
            "matches_candidate_sequence",
        )
        self.assertEqual(compatible["accuracy_evaluated_bases"], 4)
        self.assertEqual(compatible["sequence_resolved_unique_reference_bases"], 4)

        mosaic = benchmark.score(
            self._ambiguous_pool("ACAT", 2),
            ["ACAA", "AAGT"],
            [(0, 0, 4, False), (1, 0, 4, False)],
        )
        evidence = mosaic["contig_evidence"][0]
        self.assertEqual(evidence["sequence_category"], "diagnostic_mosaic")
        self.assertEqual(evidence["diagnostic_positions_0based"], [1, 2, 3])
        self.assertEqual(mosaic["diagnostic_mosaic_contigs"], 1)
        self.assertEqual(mosaic["diagnostic_mosaic_bases"], 4)
        self.assertEqual(mosaic["accuracy_evaluated_bases"], 0)
        self.assertEqual(mosaic["accuracy_unresolved_bases"], 4)

    def test_score_separates_injected_damage_from_unexplained_error(self):
        pool = self._ambiguous_pool("TAAA", 2)
        references = ["CAAA", "GAAA"]
        origins = [(0, 0, 4, False), (1, 0, 4, False)]

        damaged = benchmark.score(
            pool,
            references,
            origins,
            damage_events=[(0, 0, "C", "T")],
        )
        evidence = damaged["contig_evidence"][0]
        self.assertEqual(evidence["sequence_category"], "matches_damage_event_candidate")
        self.assertEqual(evidence["damage_compatible_transform_count"], 1)
        self.assertEqual(evidence["candidate_mismatch_lower_bound"], 1)
        self.assertEqual(evidence["candidate_mismatch_upper_bound"], 1)
        self.assertEqual(damaged["accuracy_evaluated_bases"], 4)
        self.assertEqual(damaged["accuracy_mismatches"], 1)
        self.assertEqual(damaged["accuracy_unresolved_bases"], 0)
        self.assertEqual(damaged["sequence_resolved_unique_reference_bases"], 4)

        error = benchmark.score(pool, references, origins, damage_events=[])
        evidence = error["contig_evidence"][0]
        self.assertEqual(evidence["sequence_category"], "incompatible_candidate_sequences")
        self.assertEqual(evidence["damage_compatible_transform_count"], 0)
        self.assertEqual(error["accuracy_evaluated_bases"], 0)
        self.assertEqual(error["accuracy_unresolved_bases"], 4)

    def test_reverse_placement_projects_damage_event_to_contig_orientation(self):
        raw = Read("raw", "TTTA")
        pool = ProgressiveSequencePool.from_reads([raw])
        record = replace(
            pool.records[0].corrected(Read("contig", "TAAA")),
            raw_placements=(RawPlacement(0, 0, True, 0, 4),),
        )
        pool = pool.replace_record(record)
        scored = benchmark.score(
            pool,
            ["CAAA"],
            [(0, 0, 4, True)],
            damage_events=[(0, 3, "G", "A")],
        )
        transform = scored["contig_evidence"][0]["transforms"][0]
        self.assertTrue(transform["damage_compatible"])
        self.assertEqual(transform["damage_supported_mismatches"], 1)

    def test_small_real_pipeline(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = benchmark.run_case(93003, "unique", 3, 0, Path(temporary) / "case")
        self.assertEqual(result["reads"], 45)
        self.assertGreater(result["before"]["contigs"], 0)
        self.assertEqual(result["before"]["mismatches"], 0)
        self.assertEqual(result["before"], result["after"])
