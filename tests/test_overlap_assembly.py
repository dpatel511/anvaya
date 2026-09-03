import tempfile
import unittest
from pathlib import Path

import anvaya.overlap_assembly as overlap_assembly_compat
from anvaya.cli import main
from anvaya.damage_consensus import _anchor_index
from anvaya.overlap_assembly import (
    _Alignment,
    _audit_read_supported_contig_links,
    _audit_cross_cluster_recruitment,
    _consensus,
    _filter_raw_supported_extensions,
    _merge_contig_pass,
    _overlap_confidence,
    _projection_regresses,
    _raw_supports_target_mismatches,
    _unique_best_extensions,
    assemble_overlap_contigs,
    audit_two_tier_redundancy,
    write_overlap_correction_report,
)
from anvaya.overlap_graph import (
    audit_master_overlap_graph,
    audit_raw_confirmed_master_overlap_graph,
    audit_strain_safe_containment,
)
from anvaya.overlap_reclustering import audit_iterative_reclustering
from anvaya.reads import Read


class OverlapAssemblyTests(unittest.TestCase):
    def test_legacy_overlap_module_exports_delegate_to_new_modules(self) -> None:
        self.assertEqual(
            overlap_assembly_compat.audit_master_overlap_graph([]),
            audit_master_overlap_graph([]),
        )
        self.assertEqual(
            overlap_assembly_compat.audit_raw_confirmed_master_overlap_graph(
                [], []
            ),
            audit_raw_confirmed_master_overlap_graph([], []),
        )
        self.assertEqual(
            overlap_assembly_compat.audit_strain_safe_containment([], []),
            audit_strain_safe_containment([], []),
        )
        self.assertEqual(
            overlap_assembly_compat.audit_iterative_reclustering([]),
            audit_iterative_reclustering([]),
        )

    def test_master_overlap_graph_collapses_exact_linear_path(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        contigs = [
            Read("first", genome[0:12]),
            Read("second", genome[4:16]),
            Read("third", genome[8:20]),
        ]

        projected, diagnostics = audit_master_overlap_graph(
            contigs,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=4,
        )

        self.assertEqual([contig.sequence for contig in projected], [genome])
        self.assertEqual(diagnostics.linear_paths, 1)
        self.assertEqual(diagnostics.merged_contigs, 3)
        self.assertGreaterEqual(diagnostics.transitive_edges_removed, 1)

    def test_master_overlap_graph_abstains_at_branch(self) -> None:
        contigs = [
            Read("source", "ACGTTGCACTGA"),
            Read("first", "CACTGATTTT"),
            Read("second", "CACTGACCCC"),
        ]

        projected, diagnostics = audit_master_overlap_graph(
            contigs,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
        )

        self.assertCountEqual(
            [contig.sequence for contig in projected],
            [contig.sequence for contig in contigs],
        )
        self.assertGreaterEqual(diagnostics.ambiguous_ends, 1)
        self.assertEqual(diagnostics.linear_paths, 0)

    def test_raw_confirmed_master_graph_uses_supported_overlap_allele(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        first = genome[:14]
        damaged_first = first[:10] + "A" + first[11:]
        contigs = [
            Read("first", damaged_first),
            Read("second", genome[6:]),
        ]
        raw_reads = [
            Read(f"raw-{index}", genome, [30] * len(genome))
            for index in range(3)
        ]

        projected, diagnostics = audit_raw_confirmed_master_overlap_graph(
            contigs,
            raw_reads,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=0.80,
            damage_end_window=0,
        )

        self.assertEqual([contig.sequence for contig in projected], [genome])
        self.assertEqual(diagnostics.near_exact_accepted, 1)
        self.assertEqual(diagnostics.corrected_overlap_bases, 1)

    def test_raw_confirmed_master_graph_abstains_on_supported_strains(self) -> None:
        genome = "ACGTTGCACTGATCGGACCT"
        first = genome[:14]
        alternate = first[:10] + "A" + first[11:]
        second = genome[6:]
        alternate_genome = genome[:10] + "A" + genome[11:]
        raw_reads = [
            Read(f"primary-{index}", genome, [30] * len(genome))
            for index in range(2)
        ] + [
            Read(f"alternate-{index}", alternate_genome, [30] * len(genome))
            for index in range(2)
        ]

        projected, diagnostics = audit_raw_confirmed_master_overlap_graph(
            [Read("first", alternate), Read("second", second)],
            raw_reads,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_overlap=6,
            minimum_identity=0.80,
            damage_end_window=0,
        )

        self.assertCountEqual(
            [contig.sequence for contig in projected],
            [alternate, second],
        )
        self.assertGreaterEqual(diagnostics.strain_conflicts, 1)
        self.assertEqual(diagnostics.near_exact_accepted, 0)

    def test_strain_safe_containment_removes_exact_duplicate(self) -> None:
        contigs = [
            Read("primary", "ACGTTGCACTGATCGATGCT"),
            Read("duplicate", "ACGTTGCACTGATCGATGCT"),
        ]

        projected, diagnostics = audit_strain_safe_containment(
            contigs,
            [],
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_identity=0.95,
        )

        self.assertEqual(len(projected), 1)
        self.assertEqual(diagnostics.exact_duplicates, 1)
        self.assertEqual(diagnostics.removed_contigs, 1)

    def test_strain_safe_containment_protects_supported_allele(self) -> None:
        primary = "ACGTTGCACTGATCGATGCT"
        alternate = primary[:10] + "A" + primary[11:]
        qualities = [30] * len(primary)
        raw_reads = [
            Read(f"primary-{index}", primary, qualities) for index in range(3)
        ] + [
            Read(f"alternate-{index}", alternate, qualities) for index in range(2)
        ]

        projected, diagnostics = audit_strain_safe_containment(
            [Read("a-primary", primary), Read("z-alternate", alternate)],
            raw_reads,
            anchor_k=3,
            anchors_per_read=8,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=1,
            minimum_identity=0.95,
            damage_end_window=0,
        )

        self.assertEqual(len(projected), 2)
        self.assertEqual(diagnostics.strain_protected, 1)
        self.assertEqual(diagnostics.removed_contigs, 0)

    def test_strain_safe_containment_removes_unsupported_allele(self) -> None:
        primary = "ACGTTGCACTGATCGATGCT"
        alternate = primary[:10] + "A" + primary[11:]
        qualities = [30] * len(primary)
        raw_reads = [
            Read(f"primary-{index}", primary, qualities) for index in range(3)
        ]

        projected, diagnostics = audit_strain_safe_containment(
            [Read("a-primary", primary), Read("z-alternate", alternate)],
            raw_reads,
            anchor_k=3,
            anchors_per_read=8,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=1,
            minimum_identity=0.95,
            damage_end_window=0,
        )

        self.assertEqual([contig.sequence for contig in projected], [primary])
        self.assertEqual(diagnostics.unsupported_variants, 1)
        self.assertEqual(diagnostics.removed_contigs, 1)

    def test_strain_safe_containment_preserves_partial_overlap(self) -> None:
        primary = "ACGTTGCACTGATCGATGCT"
        partial = "T" + primary[:-1]

        projected, diagnostics = audit_strain_safe_containment(
            [Read("primary", primary), Read("partial", partial)],
            [],
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_identity=0.95,
            minimum_coverage=0.95,
        )

        self.assertEqual(len(projected), 2)
        self.assertGreaterEqual(diagnostics.insufficient_evidence, 1)

    def test_strain_safe_containment_ignores_terminal_damage_support(self) -> None:
        primary = "ACGTTGCACTCATCGATGCT"
        alternate = primary[:10] + "T" + primary[11:]
        full_quality = [30] * len(primary)
        terminal_alternate = alternate[10:]
        raw_reads = [
            Read(f"primary-{index}", primary, full_quality) for index in range(3)
        ] + [
            Read(
                f"terminal-alternate-{index}",
                terminal_alternate,
                [30] * len(terminal_alternate),
            )
            for index in range(2)
        ]

        projected, diagnostics = audit_strain_safe_containment(
            [Read("a-primary", primary), Read("z-alternate", alternate)],
            raw_reads,
            anchor_k=3,
            anchors_per_read=8,
            minimum_anchor_matches=1,
            minimum_identity=0.95,
            damage_end_window=1,
        )

        self.assertEqual([contig.sequence for contig in projected], [primary])
        self.assertEqual(diagnostics.unsupported_variants, 1)

    def test_iterative_reclustering_preserves_pool_until_convergence(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        reads = [Read("center", center)]
        reads.extend(
            Read(f"extension-{index}", center[5:] + "TTGCA")
            for index in range(5)
        )

        projected, diagnostics = audit_iterative_reclustering(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            minimum_cluster_size=2,
            minimum_consensus_support=2,
            reciprocal_best_extension=False,
            maximum_iterations=3,
            minimum_output_length=1,
        )

        self.assertTrue(diagnostics.converged)
        self.assertEqual(diagnostics.iterations, 2)
        self.assertEqual(diagnostics.extended_sequences, 1)
        self.assertEqual(diagnostics.projected_contigs, 1)
        self.assertEqual(diagnostics.projected_longest_contig, len(center) + 5)
        self.assertEqual(len(projected), 1)

    def test_iterative_reclustering_gates_derived_extensions_by_raw_support(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        extension = center[5:] + "TTGCA"
        candidates = [
            _Alignment(index, extension, 5, 2)
            for index in range(3)
        ]

        retained, rejected = _filter_raw_supported_extensions(
            center,
            candidates,
            [candidates[0]],
            [0, 1, 2],
            [True, False, True],
            2,
        )

        self.assertEqual(retained, [])
        self.assertEqual(rejected, 1)

    def test_raw_support_keeps_original_evidence_after_layout_extension(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        extension = center[5:] + "TTGCA"
        candidates = [
            _Alignment(index, extension, 5, 2)
            for index in range(3)
        ]

        retained, rejected = _filter_raw_supported_extensions(
            center,
            candidates,
            [candidates[0]],
            [0, 1, 2],
            [False, False, False],
            2,
        )

        self.assertEqual(retained, [candidates[0]])
        self.assertEqual(rejected, 0)

    def test_raw_majority_confirms_one_derived_overlap_mismatch(self) -> None:
        target = "ACGTTGCACTGATCGATGCT"
        candidate = _Alignment(0, target[:10] + "A" + target[11:] + "TT", 0, 2)
        raw_candidates = [
            _Alignment(index, target, 0, 2)
            for index in range(1, 4)
        ]
        raw_candidates.append(candidate)

        self.assertTrue(_raw_supports_target_mismatches(
            target,
            candidate,
            raw_candidates,
            [Read("candidate", candidate.sequence)] + [
                Read(f"raw-{index}", target) for index in range(1, 4)
            ],
            [0, 1, 2, 3],
            3,
            2,
            20,
        ))

    def test_raw_majority_rejects_multiple_derived_overlap_mismatches(self) -> None:
        target = "ACGTTGCACTGATCGATGCT"
        sequence = target[:5] + "A" + target[6:10] + "A" + target[11:] + "TT"
        candidate = _Alignment(0, sequence, 0, 2)

        self.assertFalse(_raw_supports_target_mismatches(
            target,
            candidate,
            [_Alignment(index, target, 0, 2) for index in range(1, 4)],
            [Read("candidate", sequence)] + [
                Read(f"raw-{index}", target) for index in range(1, 4)
            ],
            [0, 1, 2, 3],
            3,
            2,
            20,
        ))

    def test_raw_mismatch_support_ignores_low_quality_bases(self) -> None:
        target = "ACGTTGCACTGATCGATGCT"
        candidate = _Alignment(0, target[:10] + "A" + target[11:] + "TT", 0, 2)
        raw_candidates = [
            _Alignment(index, target, 0, 2)
            for index in range(1, 4)
        ]
        raw_reads = [Read("candidate", candidate.sequence)] + [
            Read(f"raw-{index}", target, (10,) * len(target))
            for index in range(1, 4)
        ]

        self.assertFalse(_raw_supports_target_mismatches(
            target,
            candidate,
            raw_candidates,
            raw_reads,
            [0, 1, 2, 3],
            3,
            2,
            20,
        ))

    def test_projection_rollback_requires_loss_without_n50_gain(self) -> None:
        self.assertTrue(_projection_regresses(100, 1_000, [100] * 9))
        self.assertFalse(_projection_regresses(100, 1_000, [101] * 9))
        self.assertFalse(_projection_regresses(100, 1_000, [100] * 10))

    def test_two_tier_audit_separates_redundancy_from_novel_rescue(self) -> None:
        primary_sequence = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        novel_sequence = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTT"
        primary = [Read("primary", primary_sequence)]
        rescue = [
            Read("extends-primary", primary_sequence + "TTGCA"),
            Read("novel-a", novel_sequence),
            Read("novel-b", novel_sequence),
            Read("contained", primary_sequence[2:-2]),
        ]

        diagnostics = audit_two_tier_redundancy(
            primary,
            rescue,
            anchor_k=5,
            anchors_per_read=12,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=2,
        )

        self.assertEqual(diagnostics.rescue_contigs, 4)
        self.assertEqual(diagnostics.contained_by_primary, 1)
        self.assertEqual(diagnostics.redundant_with_rescue, 1)
        self.assertEqual(diagnostics.extends_primary, 1)
        self.assertEqual(diagnostics.ambiguous_primary_extensions, 0)
        self.assertEqual(diagnostics.replacement_contigs, 1)
        self.assertEqual(diagnostics.replacement_added_bases, 5)
        self.assertEqual(
            diagnostics.replacement_projected_longest_contig,
            len(primary_sequence) + 5,
        )
        self.assertEqual(diagnostics.novel_contigs, 1)
        self.assertEqual(
            diagnostics.projected_bases,
            len(primary_sequence) + len(novel_sequence),
        )
        self.assertEqual(
            diagnostics.projected_longest_contig,
            max(len(primary_sequence), len(novel_sequence)),
        )

    def test_two_tier_audit_rejects_multi_primary_extension(self) -> None:
        sequence = (
            "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTATTCGAGTCCGATACGTGACCTAG"
        )
        primary = [
            Read("primary-left", sequence[:35]),
            Read("primary-right", sequence[15:50]),
        ]

        diagnostics = audit_two_tier_redundancy(
            primary,
            [Read("ambiguous-rescue", sequence[:55])],
            anchor_k=5,
            anchors_per_read=16,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=2,
        )

        self.assertEqual(diagnostics.extends_primary, 1)
        self.assertEqual(diagnostics.ambiguous_primary_extensions, 1)
        self.assertEqual(diagnostics.replacement_contigs, 0)
        self.assertEqual(diagnostics.replacement_added_bases, 0)

    def test_damage_aware_ranking_can_prefer_a_cleaner_overlap(self) -> None:
        target = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        clean = _Alignment(1, target[8:] + "TTGCA", 8, 2)
        noisy_sequence = list(target[3:] + "AACGT")
        for position in (6, 13, 20):
            noisy_sequence[position] = (
                "A" if noisy_sequence[position] != "A" else "C"
            )
        noisy = _Alignment(2, "".join(noisy_sequence), 3, 2)

        length_ranked, _, _ = _unique_best_extensions(
            target,
            [clean, noisy],
            minimum_overlap_margin=3,
        )
        confidence_ranked, _, _ = _unique_best_extensions(
            target,
            [clean, noisy],
            minimum_overlap_margin=3,
            damage_aware_ranking=True,
        )

        self.assertEqual([item.read_index for item in length_ranked], [2])
        self.assertEqual([item.read_index for item in confidence_ranked], [1])

    def test_terminal_damage_transition_has_a_smaller_confidence_penalty(self) -> None:
        target = "CCGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        damage_like = _Alignment(1, "T" + target[1:] + "A", 0, 2)
        other = _Alignment(2, "G" + target[1:] + "A", 0, 2)

        damage_confidence = _overlap_confidence(
            target,
            damage_like,
            damage_mismatch_penalty=0.25,
            damage_end_window=5,
        )
        other_confidence = _overlap_confidence(
            target,
            other,
            damage_mismatch_penalty=0.25,
            damage_end_window=5,
        )

        self.assertGreater(damage_confidence, other_confidence)

    def test_damage_aware_ranking_requires_ranked_extension(self) -> None:
        with self.assertRaisesRegex(ValueError, "requires ranked_extension"):
            assemble_overlap_contigs(
                [Read("read", "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA")],
                damage_aware_ranking=True,
            )

    def test_cross_cluster_recruitment_audit_requires_frozen_contigs(self) -> None:
        with self.assertRaisesRegex(ValueError, "maximum_contig_iterations=0"):
            assemble_overlap_contigs(
                [Read("read", "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA")],
                cross_cluster_recruitment_audit=True,
            )

    def test_contig_link_audit_requires_frozen_contigs(self) -> None:
        with self.assertRaisesRegex(ValueError, "maximum_contig_iterations=0"):
            assemble_overlap_contigs(
                [Read("read", "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA")],
                read_supported_contig_link_audit=True,
            )

    def test_layout_retains_strong_internal_consensus(self) -> None:
        target = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        alternative = target[:10] + "A" + target[11:]
        members = [
            _Alignment(index, alternative, 0, 2)
            for index in range(5)
        ]

        result = _consensus(
            target,
            members,
            minimum_extension_support=5,
            dominance_ratio=4.0,
        )

        self.assertEqual(result.sequence, alternative)

    def test_polishes_supported_terminal_deamination(self) -> None:
        genome = "AACGTCCGTTGCACTGATCGATGCTAGCTACGATGGCCTATTGCA"
        truth = genome[5:40]
        damaged = "T" + truth[1:]
        reads = [Read("damaged-center", damaged)]
        reads.extend(Read(f"internal-{index}", genome[:35]) for index in range(3))
        reads.extend(Read(f"terminal-{index}", damaged) for index in range(2))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
        )

        self.assertEqual([contig.sequence for contig in contigs], [truth])
        self.assertEqual(summary.corrected_bases, 1)

    def test_records_accepted_correction_for_truth_audit(self) -> None:
        genome = "AACGTCCGTTGCACTGATCGATGCTAGCTACGATGGCCTATTGCA"
        truth = genome[5:40]
        reads = [Read("damaged-center", "T" + truth[1:])]
        reads.extend(Read(f"internal-{index}", genome[:35]) for index in range(3))
        reads.extend(Read(f"terminal-{index}", "T" + truth[1:]) for index in range(2))
        events = []

        assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
            correction_events=events,
        )

        self.assertEqual(len(events), 1)
        self.assertEqual((events[0].original_base, events[0].corrected_base), ("T", "C"))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "corrections.tsv"
            write_overlap_correction_report(events, path)
            report = path.read_text(encoding="utf-8")
        self.assertIn("sequence_before\tsequence_after", report.splitlines()[0])
        self.assertIn("overlap_contig_1", report)

    def test_polishing_does_not_change_frozen_layout(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        reads = [Read("center", center)]
        reads.extend(
            Read(f"right-{index}", center[5:] + "TTGCA")
            for index in range(5)
        )

        raw, raw_summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            damage_end_window=0,
        )
        polished, polished_summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            damage_end_window=5,
        )

        self.assertEqual(
            [len(contig.sequence) for contig in polished],
            [len(contig.sequence) for contig in raw],
        )
        self.assertEqual(polished_summary.output_contigs, raw_summary.output_contigs)
        self.assertEqual(polished_summary.contig_merges, raw_summary.contig_merges)
        self.assertEqual(polished_summary.added_bases, raw_summary.added_bases)

    def test_rejects_terminal_only_correction_support(self) -> None:
        truth = "CCGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        damaged = "T" + truth[1:]
        reads = [Read("damaged-center", damaged)]
        reads.extend(Read(f"truth-terminal-{index}", truth) for index in range(3))
        reads.extend(Read(f"damaged-terminal-{index}", damaged) for index in range(2))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
        )

        self.assertEqual([contig.sequence for contig in contigs], [damaged])
        self.assertEqual(summary.corrected_bases, 0)
        self.assertEqual(summary.correction_terminal_only, 1)

    def test_preserves_non_damage_terminal_difference(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        alternative = "G" + center[1:]
        reads = [Read("center", center)]
        reads.extend(Read(f"alternative-{index}", alternative) for index in range(3))
        reads.extend(Read(f"center-{index}", center) for index in range(2))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
        )

        self.assertEqual([contig.sequence for contig in contigs], [center])
        self.assertEqual(summary.corrected_bases, 0)

    def test_preserves_damage_like_difference_outside_terminal_window(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        alternative = center[:9] + "C" + center[10:]
        self.assertEqual(center[9], "T")
        reads = [Read("center", center)]
        reads.extend(Read(f"alternative-{index}", alternative) for index in range(3))
        reads.extend(Read(f"center-{index}", center) for index in range(2))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
        )

        self.assertEqual([contig.sequence for contig in contigs], [center])
        self.assertEqual(summary.corrected_bases, 0)

    def test_reuses_bridge_reads_then_merges_neighboring_contigs(self) -> None:
        genome = (
            "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTATTCGAGTCCGATACGTGACCTAG"
            "GATCCGTACGATCGTACGA"
        )
        reads = [Read("left", genome[:35]), Read("right", genome[30:65])]
        reads.extend(
            Read(f"bridge-{index}", genome[15:50])
            for index in range(5)
        )

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_rounds=2,
            maximum_contig_iterations=2,
        )

        self.assertEqual([contig.sequence for contig in contigs], [genome[:65]])
        self.assertEqual(summary.clusters, 2)
        self.assertEqual(summary.contig_merges, 1)

    def test_reports_clusters_below_minimum_support(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        reads = [Read("center", center)]
        reads.extend(Read(f"member-{index}", center) for index in range(2))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            minimum_cluster_size=5,
        )

        self.assertEqual(contigs, [])
        self.assertGreater(summary.clusters_below_minimum_size, 0)
        self.assertGreater(summary.candidate_alignments, 0)

    def test_ranked_extension_uses_unique_candidate_from_supported_cluster(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        reads = [Read("center", center)]
        reads.extend(Read(f"contained-{index}", center) for index in range(4))
        reads.append(Read("extension", center[5:] + "TTGCA"))

        baseline, _ = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
        )
        ranked, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
            ranked_extension=True,
        )

        self.assertEqual([contig.sequence for contig in baseline], [center])
        self.assertEqual([contig.sequence for contig in ranked], [center + "TTGCA"])
        self.assertEqual(summary.added_bases, 5)

    def test_extension_consensus_recalls_supported_appended_base(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        reads = [Read("center", center)]
        reads.extend(
            Read(f"truth-{index}", center[5:] + "TTGCA")
            for index in range(5)
        )
        reads.append(Read("long-error", center[5:] + "TCGCAA"))

        contigs, _ = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
            ranked_extension=True,
            extension_consensus=True,
        )

        self.assertEqual([contig.sequence for contig in contigs], [center + "TTGCAA"])

    def test_ranked_extension_can_require_boundary_support(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        reads = [Read("center", center)]
        reads.extend(Read(f"contained-{index}", center) for index in range(4))
        reads.append(Read("unsupported", center[5:] + "TTGCA"))

        contigs, _ = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
            ranked_extension=True,
            minimum_ranked_extension_support=2,
        )

        self.assertEqual([contig.sequence for contig in contigs], [center])

    def test_reciprocal_best_preserves_consistent_ranked_extension(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        reads = [Read("center", center)]
        reads.extend(Read(f"contained-{index}", center) for index in range(4))
        reads.append(Read("extension", center[5:] + "TTGCA"))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
            ranked_extension=True,
            reciprocal_best_extension=True,
        )

        self.assertEqual([contig.sequence for contig in contigs], [center + "TTGCA"])
        self.assertEqual(summary.reciprocal_extension_checks, 1)
        self.assertEqual(summary.reciprocal_extension_rejections, 0)

    def test_reciprocal_best_rejects_conflicting_predecessor(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        extension = center[5:] + "TTGCA"
        reads = [Read("center", center)]
        reads.extend(Read(f"contained-{index}", center) for index in range(4))
        reads.append(Read("extension", extension))
        reads.append(Read("conflicting-predecessor", "GGGGG" + extension[:33]))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_contig_iterations=0,
            ranked_extension=True,
            reciprocal_best_extension=True,
        )

        self.assertEqual([contig.sequence for contig in contigs], [center])
        self.assertEqual(summary.reciprocal_extension_checks, 1)
        self.assertEqual(summary.reciprocal_extension_rejections, 1)

    def test_contig_phase_merges_a_unique_best_overlap(self) -> None:
        left = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        right = left[15:] + "TTGCAACGT"

        merged, merges, rounds, added, ambiguous = _merge_contig_pass(
            [Read("left", left), Read("right", right)],
            anchor_k=7,
            anchors_per_read=8,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=2,
            minimum_overlap=20,
            minimum_identity=0.90,
            minimum_ry_identity=0.99,
            maximum_rounds=3,
            minimum_overlap_margin=3,
        )

        self.assertEqual([contig.sequence for contig in merged], [left + "TTGCAACGT"])
        self.assertEqual((merges, rounds, added, ambiguous), (1, 1, 9, 0))

    def test_contig_phase_stops_at_equal_competing_extensions(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        first = center[15:] + "TTGCAACGT"
        second = center[15:] + "CCATGGTAC"

        merged, merges, _, _, ambiguous = _merge_contig_pass(
            [Read("center", center), Read("first", first), Read("second", second)],
            anchor_k=7,
            anchors_per_read=8,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=2,
            minimum_overlap=20,
            minimum_identity=0.90,
            minimum_ry_identity=0.99,
            maximum_rounds=3,
            minimum_overlap_margin=3,
        )

        self.assertNotIn(center + "TTGCAACGT", [contig.sequence for contig in merged])
        self.assertNotIn(center + "CCATGGTAC", [contig.sequence for contig in merged])
        self.assertEqual(merges, 0)
        self.assertGreaterEqual(ambiguous, 1)

    def test_emits_a_directly_extended_contig(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        extension = "TTGCA"
        reads = [Read("center", center)]
        reads.extend(
            Read(f"right-{index}", center[5:] + extension)
            for index in range(5)
        )

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
        )

        self.assertEqual([contig.sequence for contig in contigs], [center + extension])
        self.assertEqual(summary.output_contigs, 1)
        self.assertEqual(summary.clustered_reads, 6)
        self.assertEqual(summary.added_bases, 5)

    def test_does_not_emit_unsupported_singletons(self) -> None:
        reads = [Read("first", "AACCGGTTAACCGG"), Read("second", "TTTTGGGGCCCCAA")]

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=5,
            minimum_overlap=10,
        )

        self.assertEqual(contigs, [])
        self.assertEqual(summary.clusters, 0)

    def test_reindexes_an_extended_end_for_another_round(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        first = center + "TTGCA"
        reads = [Read("center", center)]
        reads.extend(
            Read(f"round-one-{index}", center[5:] + "TTGCA")
            for index in range(5)
        )
        reads.extend(
            Read(f"round-two-{index}", first[20:] + "AACGT")
            for index in range(5)
        )

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
            maximum_rounds=2,
        )

        self.assertEqual(
            [contig.sequence for contig in contigs],
            [center + "TTGCAAACGT"],
        )
        self.assertEqual(summary.extension_rounds, 2)
        self.assertEqual(summary.added_bases, 10)

    def test_audits_consensus_supported_deferred_end_recruitment(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        first_extension = "TTGCAACGTCCGATA"
        first = center + first_extension
        reads = [Read("accepted", first[10:])]
        reads.extend(
            Read(f"deferred-{index}", first[25:] + "GACTA")
            for index in range(5)
        )
        diagnostics = _audit_cross_cluster_recruitment(
            [Read("frozen", first)],
            [{0}],
            reads,
            [0, None, None, None, None, None],
            list(range(len(reads))),
            _anchor_index(reads, 7, 0, 32, 100),
            anchor_k=7,
            anchors_per_read=32,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=1,
            minimum_overlap=25,
            minimum_identity=0.90,
            minimum_ry_identity=0.99,
            position_bits=max(len(read.sequence) for read in reads).bit_length(),
            target_window=max(len(read.sequence) for read in reads),
            minimum_overlap_margin=3,
            minimum_consensus_support=5,
            dominance_ratio=4.0,
            damage_aware_ranking=True,
            minimum_confidence_margin=0.0,
            damage_mismatch_penalty=0.25,
            ranking_damage_end_window=5,
        )

        self.assertEqual(diagnostics.deferred_reads, 5)
        self.assertEqual(diagnostics.deferred_candidate_reads, 5)
        self.assertEqual(diagnostics.deferred_extension_candidates, 5)
        self.assertEqual(diagnostics.deferred_ambiguous_assignments, 0)
        self.assertEqual(diagnostics.deferred_unique_assignments, 5)
        self.assertEqual(diagnostics.deferred_reciprocal_assignments, 5)
        self.assertEqual(diagnostics.supported_contigs, 1)
        self.assertEqual(diagnostics.supported_bases, 5)
        self.assertEqual(diagnostics.deferred_supported_contigs, 1)
        self.assertEqual(diagnostics.deferred_supported_bases, 5)

    def test_audits_consensus_supported_cross_cluster_recruitment(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        first = center + "TTGCAACGTCCGATA"
        bridge = first[25:] + "GACTA"
        reads = [Read("target-member", first[10:])]
        reads.extend(Read(f"bridge-{index}", bridge) for index in range(5))

        diagnostics = _audit_cross_cluster_recruitment(
            [Read("target", first), Read("source", bridge)],
            [{0}, {1, 2, 3, 4, 5}],
            reads,
            [0, 1, 1, 1, 1, 1],
            list(range(len(reads))),
            _anchor_index(reads, 7, 0, 32, 100),
            anchor_k=7,
            anchors_per_read=32,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=1,
            minimum_overlap=25,
            minimum_identity=0.90,
            minimum_ry_identity=0.99,
            position_bits=max(len(read.sequence) for read in reads).bit_length(),
            target_window=max(len(read.sequence) for read in reads),
            minimum_overlap_margin=3,
            minimum_consensus_support=5,
            dominance_ratio=4.0,
            damage_aware_ranking=True,
            minimum_confidence_margin=0.0,
            damage_mismatch_penalty=0.25,
            ranking_damage_end_window=5,
        )

        self.assertEqual(diagnostics.deferred_reads, 0)
        self.assertGreaterEqual(diagnostics.cross_cluster_extension_candidates, 5)
        self.assertGreaterEqual(diagnostics.cross_cluster_unique_assignments, 5)
        self.assertGreaterEqual(diagnostics.cross_cluster_reciprocal_assignments, 5)
        self.assertEqual(diagnostics.supported_bases, 5)
        self.assertEqual(diagnostics.cross_cluster_supported_contigs, 1)
        self.assertEqual(diagnostics.cross_cluster_supported_bases, 5)

    def test_cross_cluster_recruitment_audit_does_not_change_output(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        reads = [Read("center", center)]
        reads.extend(Read(f"contained-{index}", center) for index in range(4))
        reads.append(Read("extension", center[5:] + "TTGCA"))
        parameters = {
            "anchor_k": 7,
            "minimum_overlap": 20,
            "maximum_contig_iterations": 0,
            "ranked_extension": True,
        }

        baseline, _ = assemble_overlap_contigs(reads, **parameters)
        audited, _ = assemble_overlap_contigs(
            reads,
            cross_cluster_recruitment_audit=True,
            **parameters,
        )

        self.assertEqual(audited, baseline)

    def test_audits_read_supported_reciprocal_contig_link(self) -> None:
        left = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        extension = "TTGCAACGTCCGATA"
        right = left[15:] + extension
        merged = left + extension
        reads = [
            Read("left-bridge", merged[10:45]),
            Read("right-bridge", merged[10:45]),
        ]

        diagnostics = _audit_read_supported_contig_links(
            [Read("left", left), Read("right", right)],
            [{0}, {1}],
            reads,
            [0, 1],
            anchor_k=7,
            anchors_per_read=16,
            maximum_anchor_occurrences=100,
            minimum_anchor_matches=1,
            minimum_overlap=20,
            minimum_identity=0.90,
            minimum_ry_identity=0.99,
            minimum_overlap_margin=3,
            minimum_read_support=2,
        )

        self.assertGreaterEqual(diagnostics.candidate_overlaps, 2)
        self.assertGreaterEqual(diagnostics.reciprocal_overlaps, 1)
        self.assertEqual(diagnostics.read_supported_overlaps, 1)
        self.assertEqual(diagnostics.support_at_least_1, 1)
        self.assertEqual(diagnostics.support_at_least_2, 1)
        self.assertEqual(diagnostics.support_at_least_3, 0)
        self.assertEqual(diagnostics.support_at_least_5, 0)
        self.assertEqual(diagnostics.max_read_support, 2)
        self.assertEqual(diagnostics.ambiguous_ends, 0)
        self.assertEqual(diagnostics.cyclic_components, 0)
        self.assertEqual(diagnostics.linear_chains, 1)
        self.assertEqual(diagnostics.projected_joins, 1)
        self.assertEqual(diagnostics.projected_merged_bases, len(merged))
        self.assertEqual(diagnostics.projected_longest_contig, len(merged))

    def test_contig_link_audit_does_not_change_output(self) -> None:
        genome = (
            "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTATTCGAGTCCGATACGTGACCTAG"
            "GATCCGTACGATCGTACGA"
        )
        reads = [Read("left", genome[:35]), Read("right", genome[30:65])]
        reads.extend(Read(f"bridge-{index}", genome[15:50]) for index in range(5))
        parameters = {
            "anchor_k": 7,
            "minimum_overlap": 20,
            "maximum_rounds": 2,
            "maximum_contig_iterations": 0,
        }

        baseline, _ = assemble_overlap_contigs(reads, **parameters)
        audited, _ = assemble_overlap_contigs(
            reads,
            read_supported_contig_link_audit=True,
            **parameters,
        )

        self.assertEqual(audited, baseline)

    def test_separates_ry_inconsistent_clusters(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        overlap = list(center[5:])
        overlap[10] = "C" if overlap[10] in "AG" else "A"
        conflicting = "".join(overlap) + "TTGCA"
        reads = [Read("center", center)]
        reads.extend(Read(f"truth-{index}", center) for index in range(5))
        reads.extend(Read(f"member-{index}", conflicting) for index in range(5))

        contigs, summary = assemble_overlap_contigs(
            reads,
            anchor_k=7,
            minimum_overlap=20,
        )

        self.assertEqual(
            {contig.sequence for contig in contigs},
            {center, conflicting},
        )
        self.assertEqual(summary.clusters, 2)
        self.assertGreater(summary.candidate_ry_rejected, 0)

    def test_cli_is_separate_from_dbg_assembly(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "reads.fasta"
            output_path = root / "contigs.fasta"
            records = [("center", center)] + [
                (f"right-{index}", center[5:] + "TTGCA")
                for index in range(5)
            ]
            input_path.write_text(
                "".join(f">{name}\n{sequence}\n" for name, sequence in records),
                encoding="utf-8",
            )

            result = main([
                "overlap-assemble",
                "-i",
                str(input_path),
                "-o",
                str(output_path),
                "--min-output-length",
                "1",
            ])

            self.assertEqual(result, 0)
            self.assertIn(center + "TTGCA", output_path.read_text(encoding="utf-8"))

    def test_two_tier_audit_preserves_primary_cli_output(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "reads.fasta"
            baseline_path = root / "baseline.fasta"
            audited_path = root / "audited.fasta"
            input_path.write_text(
                "".join(
                    f">read-{index}\n{center}\n"
                    for index in range(5)
                ),
                encoding="utf-8",
            )
            common = [
                "overlap-assemble",
                "-i", str(input_path),
                "--anchor-k", "7",
                "--min-overlap", "20",
                "--max-contig-iterations", "0",
                "--min-output-length", "1",
            ]

            self.assertEqual(main(common + ["-o", str(baseline_path)]), 0)
            self.assertEqual(main(common + [
                "-o", str(audited_path),
                "--two-tier-redundancy-audit",
                "--rescue-min-cluster-size", "3",
            ]), 0)

            self.assertEqual(audited_path.read_bytes(), baseline_path.read_bytes())

    def test_iterative_reclustering_audit_preserves_primary_cli_output(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "reads.fasta"
            baseline_path = root / "baseline.fasta"
            audited_path = root / "audited.fasta"
            report_path = root / "reclustering.tsv"
            projection_path = root / "reclustering.fasta"
            input_path.write_text(
                "".join(f">read-{index}\n{center}\n" for index in range(5)),
                encoding="utf-8",
            )
            common = [
                "overlap-assemble",
                "-i", str(input_path),
                "--anchor-k", "7",
                "--min-overlap", "20",
                "--max-contig-iterations", "0",
                "--ranked-extension",
                "--reciprocal-best-extension",
                "--damage-aware-ranking",
                "--min-output-length", "1",
            ]

            self.assertEqual(main(common + ["-o", str(baseline_path)]), 0)
            self.assertEqual(main(common + [
                "-o", str(audited_path),
                "--iterative-reclustering-audit",
                "--iterative-reclustering-report", str(report_path),
                "--iterative-reclustering-projection", str(projection_path),
            ]), 0)

            self.assertEqual(audited_path.read_bytes(), baseline_path.read_bytes())
            self.assertIn("projected_n50", report_path.read_text(encoding="utf-8"))
            self.assertTrue(projection_path.exists())

    def test_master_overlap_graph_audit_preserves_other_cli_outputs(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "reads.fasta"
            baseline_path = root / "baseline.fasta"
            audited_path = root / "audited.fasta"
            reclustered_path = root / "reclustered.fasta"
            master_path = root / "master.fasta"
            raw_master_path = root / "raw-master.fasta"
            containment_path = root / "containment.fasta"
            input_path.write_text(
                "".join(f">read-{index}\n{center}\n" for index in range(5)),
                encoding="utf-8",
            )
            common = [
                "overlap-assemble",
                "-i", str(input_path),
                "--anchor-k", "7",
                "--min-overlap", "20",
                "--max-contig-iterations", "0",
                "--ranked-extension",
                "--reciprocal-best-extension",
                "--damage-aware-ranking",
                "--min-output-length", "1",
            ]

            self.assertEqual(main(common + ["-o", str(baseline_path)]), 0)
            self.assertEqual(main(common + [
                "-o", str(audited_path),
                "--iterative-reclustering-audit",
                "--iterative-reclustering-projection", str(reclustered_path),
                "--master-overlap-graph-audit",
                "--master-overlap-graph-projection", str(master_path),
                "--raw-confirmed-master-graph-audit",
                "--raw-confirmed-master-graph-projection", str(raw_master_path),
                "--strain-safe-containment-audit",
                "--strain-safe-containment-projection", str(containment_path),
            ]), 0)

            self.assertEqual(audited_path.read_bytes(), baseline_path.read_bytes())
            self.assertTrue(reclustered_path.exists())
            self.assertTrue(master_path.exists())
            self.assertTrue(raw_master_path.exists())
            self.assertTrue(containment_path.exists())

    def test_cli_accepts_candidate_discovery_parameters(self) -> None:
        center = "ACGTTGCACTGATCGATGCTAGCTACGATGGCCTA"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "reads.fasta"
            output_path = root / "contigs.fasta"
            input_path.write_text(
                "".join(
                    f">read-{index}\n{center}\n"
                    for index in range(5)
                ),
                encoding="utf-8",
            )

            result = main([
                "overlap-assemble",
                "-i", str(input_path),
                "-o", str(output_path),
                "--anchor-k", "7",
                "--anchors-per-read", "16",
                "--max-anchor-occurrences", "50",
                "--min-anchor-matches", "2",
                "--min-overlap", "20",
                "--ranked-extension",
                "--reciprocal-best-extension",
                "--damage-aware-ranking",
                "--min-overlap-confidence-margin", "0.02",
                "--damage-mismatch-penalty", "0.25",
                "--ranking-damage-end-window", "5",
                "--max-contig-iterations", "0",
                "--cross-cluster-recruitment-audit",
                "--read-supported-contig-link-audit",
                "--min-contig-link-read-support", "2",
            ])

            self.assertEqual(result, 0)
            self.assertTrue(output_path.exists())


if __name__ == "__main__":
    unittest.main()
