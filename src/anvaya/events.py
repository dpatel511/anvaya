"""Sequence and read-end evidence reports for graph alternatives."""

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from anvaya.bidirected import (
    BidirectedDeBruijnGraph,
    _edge_support_total,
    _successor,
    collect_target_edge_links,
    spell_edge_path,
)
from anvaya.bubbles import Bubble, BubblePath
from anvaya.cleaning import TipCandidate
from anvaya.classification import format_substitutions
from anvaya.damage_profile import infer_damage_profile, write_damage_profile
from anvaya.damage_projection import decide_damage_projection
from anvaya.event_likelihood import (
    CrossFittedDamageModels,
    fit_cross_fitted_damage_models,
    score_matched_event,
)
from anvaya.incomplete_branches import (
    IncompleteBranchCandidate,
    match_incomplete_branch_to_backbone,
)
from anvaya.tip_classification import classify_tip_match
from anvaya.tip_matching import (
    TipBackboneMatch,
    match_competing_paths,
    match_tip_to_backbone,
)


@dataclass(slots=True, frozen=True)
class EventReportSummary:
    """Counts written to one graph-event report."""

    tips: int = 0
    matched_tips: int = 0
    incomplete_branches: int = 0
    matched_incomplete_branches: int = 0
    bubbles: int = 0
    paths: int = 0
    projection_candidates: int = 0


_TIP_MATCH_FIELDS = (
    "tip_match_sequence",
    "backbone_matched",
    "backbone_sequence",
    "backbone_edge_count",
    "backbone_observations",
    "backbone_minimum_edge_support",
    "relative_coverage",
    "sequence_identity",
    "ry_identity",
)

_TIP_CLASSIFICATION_FIELDS = (
    "tip_classification",
    "tip_classification_reasons",
    "tip_damage_score",
    "tip_error_score",
    "tip_variation_score",
    "tip_terminal_fraction",
    "tip_strand_balance",
    "tip_terminal_enrichment",
    "substitution_terminal_observations",
    "substitution_other_terminal_observations",
    "substitution_internal_observations",
    "substitution_terminal_fraction",
    "mean_damage_distance",
    "molecule_links_collected",
    "joint_molecule_observations",
    "joint_molecule_fraction",
    "quality_observations",
    "missing_quality_observations",
    "mean_base_quality",
    "sequencing_error_log_likelihood",
    "likelihood_status",
    "likelihood_reasons",
    "likelihood_profile_scope",
    "likelihood_conditioning_scope",
    "likelihood_observations",
    "likelihood_alternative_observations",
    "likelihood_reference_observations",
    "likelihood_missing_quality_observations",
    "damage_model_log_likelihood",
    "damage_model_log_likelihood_min",
    "damage_model_log_likelihood_max",
    "damage_model_log_likelihood_spread",
    "damage_model_conditioning_log_probability",
    "error_model_log_likelihood",
    "error_model_conditioning_log_probability",
    "variation_model_log_likelihood",
    "variation_model_conditioning_log_probability",
    "variation_model_frequency",
    "variation_model_complexity_penalty",
    "variation_model_penalized_log_likelihood",
    "likelihood_best_explanation",
    "likelihood_log_margin",
    "damage_error_log_contrast",
    "damage_variation_log_contrast",
    "error_variation_log_contrast",
    "projection_decision",
    "projection_reasons",
    "projection_target_sequence",
)


def _oriented_edge_evidence(
    graph: BidirectedDeBruijnGraph,
    start: int,
    edge_ids: Sequence[int],
) -> list[tuple[int, int, int, int]]:
    values: list[tuple[int, int, int, int]] = []
    current = start

    for edge_id in edge_ids:
        evidence = graph.end_support(edge_id)
        edge_left = sum(evidence.left)
        edge_right = sum(evidence.right)
        if current != graph.edge_sources[edge_id]:
            edge_left, edge_right = edge_right, edge_left
        values.append(
            (
                edge_left,
                edge_right,
                evidence.internal,
                evidence.ambiguous_terminal,
            )
        )
        current = _successor(graph, current, edge_id)

    return values


def _path_evidence(
    graph: BidirectedDeBruijnGraph,
    start: int,
    edge_ids: Sequence[int],
) -> tuple[int, int, int, int]:
    values = _oriented_edge_evidence(graph, start, edge_ids)
    return tuple(sum(value[index] for value in values) for index in range(4))


def _substitutions(reference: str, alternative: str) -> list[tuple[int, str, str]]:
    if len(reference) != len(alternative):
        return []
    return [
        (position, reference_base, alternative_base)
        for position, (reference_base, alternative_base) in enumerate(
            zip(reference, alternative, strict=True),
            start=1,
        )
        if reference_base != alternative_base
    ]


def _format_substitutions(
    substitutions: Sequence[tuple[int, str, str]],
) -> str:
    return ";".join(
        f"{position}:{reference}>{alternative}"
        for position, reference, alternative in substitutions
    )


def _damage_compatibility(
    graph: BidirectedDeBruijnGraph,
    start: int,
    edge_ids: Sequence[int],
    substitutions: Sequence[tuple[int, str, str]],
) -> str:
    if not substitutions:
        return "unknown"
    edge_evidence = _oriented_edge_evidence(graph, start, edge_ids)

    def is_compatible(
        position: int,
        reference: str,
        alternative: str,
    ) -> bool:
        edge_index = position - graph.node_length - 1
        if not 0 <= edge_index < len(edge_evidence):
            return False
        left_terminal, right_terminal, _, _ = edge_evidence[edge_index]
        return (
            reference == "C"
            and alternative == "T"
            and left_terminal > 0
        ) or (
            reference == "G"
            and alternative == "A"
            and right_terminal > 0
        )

    compatible = all(
        is_compatible(position, reference, alternative)
        for position, reference, alternative in substitutions
    )
    return str(compatible).lower()


def _bubble_rows(
    graph: BidirectedDeBruijnGraph,
    event_id: str,
    bubble: Bubble,
    matches: Sequence[TipBackboneMatch | None],
    likelihood_models: CrossFittedDamageModels | None,
    edge_links: dict[int, dict[int, int | None]] | None = None,
) -> list[list[object]]:
    sequences = [
        spell_edge_path(graph, bubble.start, path.edge_ids)
        for path in bubble.paths
    ]
    reference_index = max(
        range(len(bubble.paths)),
        key=lambda index: (
            bubble.paths[index].minimum_edge_support,
            bubble.paths[index].observations,
            -index,
        ),
    )
    reference_sequence = sequences[reference_index]
    rows: list[list[object]] = []

    for index, (path, sequence) in enumerate(
        zip(bubble.paths, sequences, strict=True)
    ):
        if index != reference_index:
            match = matches[index]
            if match is not None:
                rows.append(
                    _matched_path_row(
                        graph,
                        event_id,
                        "bubble",
                        bubble.start,
                        bubble.end,
                        path.edge_ids,
                        path.observations,
                        path.minimum_edge_support,
                        sequence,
                        match,
                        likelihood_models,
                        edge_links,
                        path_index=index,
                    )
                )
                continue
        left, right, internal, ambiguous = _path_evidence(
            graph,
            bubble.start,
            path.edge_ids,
        )
        substitutions = (
            []
            if index == reference_index
            else _substitutions(reference_sequence, sequence)
        )
        compatibility = (
            "reference"
            if index == reference_index
            else _damage_compatibility(
                graph,
                bubble.start,
                path.edge_ids,
                substitutions,
            )
        )
        rows.append(
            [
                event_id,
                "bubble",
                index,
                str(index == reference_index).lower(),
                bubble.start,
                bubble.end,
                len(path.edge_ids),
                path.observations,
                path.minimum_edge_support,
                left,
                right,
                internal,
                ambiguous,
                sequence,
                _format_substitutions(substitutions),
                compatibility,
            ]
            + [""]
            * (len(_TIP_MATCH_FIELDS) + len(_TIP_CLASSIFICATION_FIELDS))
        )
    return rows


def _match_columns(
    graph: BidirectedDeBruijnGraph,
    match: TipBackboneMatch | None,
    likelihood_models: CrossFittedDamageModels | None,
    edge_links: dict[int, dict[int, int | None]] | None = None,
) -> tuple[str, str, list[object]]:
    if match is None:
        return (
            "",
            "unknown",
            [""]
            + ["false"]
            + [""] * (len(_TIP_MATCH_FIELDS) - 2)
            + [""] * len(_TIP_CLASSIFICATION_FIELDS),
        )

    decision = classify_tip_match(graph, match)
    likelihood = score_matched_event(graph, match, likelihood_models, edge_links)
    projection = decide_damage_projection(match, decision, likelihood)
    columns: list[object] = [
        match.tip_sequence,
        "true",
        match.backbone_sequence,
        len(match.backbone_edge_ids),
        match.backbone_observations,
        match.backbone_minimum_edge_support,
        f"{match.relative_coverage:.6f}",
        f"{match.sequence_identity:.6f}",
        f"{match.ry_identity:.6f}",
        decision.label,
        ";".join(decision.reasons),
        f"{decision.damage_score:.6f}",
        f"{decision.error_score:.6f}",
        f"{decision.variation_score:.6f}",
        f"{decision.evidence.tip_terminal_fraction:.6f}",
        (
            ""
            if decision.evidence.tip_strand_balance is None
            else f"{decision.evidence.tip_strand_balance:.6f}"
        ),
        f"{decision.evidence.terminal_enrichment:.6f}",
        decision.evidence.substitution_terminal_observations,
        decision.evidence.substitution_other_terminal_observations,
        decision.evidence.substitution_internal_observations,
        f"{decision.evidence.substitution_terminal_fraction:.6f}",
        (
            ""
            if decision.evidence.mean_damage_distance is None
            else f"{decision.evidence.mean_damage_distance:.6f}"
        ),
        str(decision.evidence.molecule_links_collected).lower(),
        (
            ""
            if decision.evidence.joint_molecule_observations is None
            else decision.evidence.joint_molecule_observations
        ),
        (
            ""
            if decision.evidence.joint_molecule_fraction is None
            else f"{decision.evidence.joint_molecule_fraction:.6f}"
        ),
        decision.evidence.quality_observations,
        decision.evidence.missing_quality_observations,
        (
            ""
            if decision.evidence.mean_base_quality is None
            else f"{decision.evidence.mean_base_quality:.6f}"
        ),
        (
            ""
            if decision.evidence.sequencing_error_log_likelihood is None
            else f"{decision.evidence.sequencing_error_log_likelihood:.6f}"
        ),
        likelihood.status,
        ";".join(likelihood.reasons),
        likelihood.profile_scope,
        likelihood.conditioning_scope,
        likelihood.observations,
        likelihood.alternative_observations,
        likelihood.reference_observations,
        likelihood.missing_quality_observations,
        (
            ""
            if likelihood.damage_log_likelihood is None
            else f"{likelihood.damage_log_likelihood:.6f}"
        ),
        (
            ""
            if likelihood.damage_log_likelihood_min is None
            else f"{likelihood.damage_log_likelihood_min:.6f}"
        ),
        (
            ""
            if likelihood.damage_log_likelihood_max is None
            else f"{likelihood.damage_log_likelihood_max:.6f}"
        ),
        (
            ""
            if likelihood.damage_log_likelihood_spread is None
            else f"{likelihood.damage_log_likelihood_spread:.6f}"
        ),
        (
            ""
            if likelihood.damage_conditioning_log_probability is None
            else f"{likelihood.damage_conditioning_log_probability:.6f}"
        ),
        (
            ""
            if likelihood.error_log_likelihood is None
            else f"{likelihood.error_log_likelihood:.6f}"
        ),
        (
            ""
            if likelihood.error_conditioning_log_probability is None
            else f"{likelihood.error_conditioning_log_probability:.6f}"
        ),
        (
            ""
            if likelihood.variation_log_likelihood is None
            else f"{likelihood.variation_log_likelihood:.6f}"
        ),
        (
            ""
            if likelihood.variation_conditioning_log_probability is None
            else f"{likelihood.variation_conditioning_log_probability:.6f}"
        ),
        (
            ""
            if likelihood.variation_frequency is None
            else f"{likelihood.variation_frequency:.6f}"
        ),
        (
            ""
            if likelihood.variation_complexity_penalty is None
            else f"{likelihood.variation_complexity_penalty:.6f}"
        ),
        (
            ""
            if likelihood.variation_penalized_log_likelihood is None
            else f"{likelihood.variation_penalized_log_likelihood:.6f}"
        ),
        likelihood.best_explanation or "",
        (
            ""
            if likelihood.log_likelihood_margin is None
            else f"{likelihood.log_likelihood_margin:.6f}"
        ),
        (
            ""
            if likelihood.damage_error_log_contrast is None
            else f"{likelihood.damage_error_log_contrast:.6f}"
        ),
        (
            ""
            if likelihood.damage_variation_log_contrast is None
            else f"{likelihood.damage_variation_log_contrast:.6f}"
        ),
        (
            ""
            if likelihood.error_variation_log_contrast is None
            else f"{likelihood.error_variation_log_contrast:.6f}"
        ),
        projection.decision,
        ";".join(projection.reasons),
        projection.target_sequence,
    ]
    return (
        format_substitutions(match.substitutions),
        str(match.damage_compatible).lower(),
        columns,
    )


def _matched_path_row(
    graph: BidirectedDeBruijnGraph,
    event_id: str,
    event_type: str,
    start: int,
    end: int,
    edge_ids: tuple[int, ...],
    observations: int,
    minimum_support: int,
    sequence: str,
    match: TipBackboneMatch | None,
    likelihood_models: CrossFittedDamageModels | None,
    edge_links: dict[int, dict[int, int | None]] | None = None,
    *,
    path_index: int = 0,
) -> list[object]:
    left, right, internal, ambiguous = _path_evidence(
        graph,
        start,
        edge_ids,
    )
    substitutions, compatibility, match_columns = _match_columns(
        graph,
        match,
        likelihood_models,
        edge_links,
    )
    return [
        event_id,
        event_type,
        path_index,
        "false",
        start,
        end,
        len(edge_ids),
        observations,
        minimum_support,
        left,
        right,
        internal,
        ambiguous,
        sequence,
        substitutions,
        compatibility,
        *match_columns,
    ]


def _tip_row(
    graph: BidirectedDeBruijnGraph,
    event_id: str,
    tip: TipCandidate,
    match: TipBackboneMatch | None,
    likelihood_models: CrossFittedDamageModels | None,
    edge_links: dict[int, dict[int, int | None]] | None = None,
) -> tuple[list[object], TipBackboneMatch | None]:
    observations = sum(
        _edge_support_total(graph, edge_id) for edge_id in tip.edge_ids
    )
    minimum_support = min(
        _edge_support_total(graph, edge_id) for edge_id in tip.edge_ids
    )
    return (
        _matched_path_row(
            graph,
            event_id,
            "tip",
            tip.start,
            tip.end,
            tip.edge_ids,
            observations,
            minimum_support,
            spell_edge_path(graph, tip.start, tip.edge_ids),
            match,
            likelihood_models,
            edge_links,
        ),
        match,
    )


def _incomplete_branch_row(
    graph: BidirectedDeBruijnGraph,
    event_id: str,
    candidate: IncompleteBranchCandidate,
    match: TipBackboneMatch | None,
    likelihood_models: CrossFittedDamageModels | None,
    edge_links: dict[int, dict[int, int | None]] | None = None,
) -> tuple[list[object], TipBackboneMatch | None]:
    return (
        _matched_path_row(
            graph,
            event_id,
            "incomplete_branch",
            candidate.start,
            candidate.end,
            candidate.edge_ids,
            candidate.observations,
            candidate.minimum_edge_support,
            spell_edge_path(graph, candidate.start, candidate.edge_ids),
            match,
            likelihood_models,
            edge_links,
        ),
        match,
    )


def write_event_report(
    graph: BidirectedDeBruijnGraph,
    tips: Sequence[TipCandidate],
    bubbles: Sequence[Bubble],
    path: str | Path,
    incomplete_branches: Sequence[IncompleteBranchCandidate] = (),
    damage_profile_path: str | Path | None = None,
    edge_links: dict[int, dict[int, int | None]] | None = None,
    reads: Sequence[str] | None = None,
    read_qualities: Sequence[Sequence[int] | None] | None = None,
    read_molecule_ids: Sequence[int] | None = None,
) -> EventReportSummary:
    """Write non-destructive tip, branch, and bubble-path evidence as TSV."""
    if graph.end_window == 0:
        raise ValueError("event reporting requires read-end evidence")
    if damage_profile_path is not None and not graph.molecule_links_collected:
        raise ValueError("damage profiling requires molecule links")

    tip_matches = [match_tip_to_backbone(graph, tip) for tip in tips]
    incomplete_matches = [
        match_incomplete_branch_to_backbone(graph, candidate)
        for candidate in incomplete_branches
    ]
    bubble_matches: list[tuple[TipBackboneMatch | None, ...]] = []
    for bubble in bubbles:
        reference_index = max(
            range(len(bubble.paths)),
            key=lambda index: (
                bubble.paths[index].minimum_edge_support,
                bubble.paths[index].observations,
                -index,
            ),
        )
        reference = bubble.paths[reference_index].edge_ids
        bubble_matches.append(
            tuple(
                None
                if index == reference_index
                else match_competing_paths(
                    graph,
                    bubble.start,
                    bubble_path.edge_ids,
                    reference,
                )
                for index, bubble_path in enumerate(bubble.paths)
            )
        )
    profile_matches = [match for match in tip_matches if match is not None]
    if reads is not None and edge_links is None:
        evidence_matches = list(profile_matches)
        evidence_matches.extend(
            match
            for match in incomplete_matches
            if match is not None
        )
        for matches in bubble_matches:
            evidence_matches.extend(
                match
                for match in matches
                if match is not None
            )
        target_edges: set[int] = set()
        for match in evidence_matches:
            for change in match.substitutions:
                if (change.reference, change.alternative) in {
                    ("C", "T"),
                    ("G", "A"),
                }:
                    continue
                edge_index = change.position - graph.node_length - 1
                if (
                    0 <= edge_index < len(match.tip_edge_ids)
                    and edge_index < len(match.backbone_edge_ids)
                ):
                    target_edges.add(match.tip_edge_ids[edge_index])
                    target_edges.add(match.backbone_edge_ids[edge_index])
        edge_links = collect_target_edge_links(
            graph,
            reads,
            tuple(target_edges),
            read_qualities,
            read_molecule_ids,
        )
    profile = None
    likelihood_models = None
    if damage_profile_path is not None:
        profile = infer_damage_profile(graph, profile_matches)
        likelihood_models = fit_cross_fitted_damage_models(profile)

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "event_id",
                "event_type",
                "path_index",
                "reference_path",
                "start_handle",
                "end_handle",
                "edge_count",
                "observations",
                "minimum_edge_support",
                "left_terminal",
                "right_terminal",
                "internal",
                "ambiguous_terminal",
                "sequence",
                "substitutions",
                "damage_compatible",
                *_TIP_MATCH_FIELDS,
                *_TIP_CLASSIFICATION_FIELDS,
            ]
        )
        path_count = 0
        projection_candidates = 0
        matched_tip_count = 0
        for index, (tip, match) in enumerate(
            zip(tips, tip_matches, strict=True),
            start=1,
        ):
            row, match = _tip_row(
                graph,
                f"tip-{index}",
                tip,
                match,
                likelihood_models,
                edge_links,
            )
            writer.writerow(row)
            projection_candidates += row[-3] == "would_project"
            if match is not None:
                matched_tip_count += 1
            path_count += 1
        matched_incomplete_count = 0
        for index, (candidate, match) in enumerate(
            zip(incomplete_branches, incomplete_matches, strict=True),
            start=1,
        ):
            row, match = _incomplete_branch_row(
                graph,
                f"incomplete-branch-{index}",
                candidate,
                match,
                likelihood_models,
                edge_links,
            )
            writer.writerow(row)
            projection_candidates += row[-3] == "would_project"
            if match is not None:
                matched_incomplete_count += 1
            path_count += 1
        for index, (bubble, matches) in enumerate(
            zip(bubbles, bubble_matches, strict=True),
            start=1,
        ):
            rows = _bubble_rows(
                graph,
                f"bubble-{index}",
                bubble,
                matches,
                likelihood_models,
                edge_links,
            )
            writer.writerows(rows)
            projection_candidates += sum(row[-3] == "would_project" for row in rows)
            path_count += len(rows)

    if damage_profile_path is not None:
        assert profile is not None
        assert likelihood_models is not None
        write_damage_profile(
            profile,
            damage_profile_path,
            cross_fit_summary={
                "status": likelihood_models.status,
                "reasons": likelihood_models.reasons,
                "folds": likelihood_models.folds,
                "fitted_folds": likelihood_models.fitted_folds,
                "fold_models": [
                    {
                        "fold": fold + 1,
                        "status": (
                            "fitted"
                            if fold in likelihood_models.probabilities_by_fold
                            else "insufficient_evidence"
                        ),
                        "fitted_probabilities": (
                            likelihood_models.probabilities_by_fold.get(
                                fold,
                                (),
                            )
                        ),
                    }
                    for fold in range(likelihood_models.folds)
                ],
                "independent_tip_profile_probabilities": (
                    likelihood_models.independent_probabilities
                ),
                "assignment": (
                    "sorted unique candidate loci assigned round-robin"
                ),
                "scoring_rule": (
                    "each tip locus uses a model fitted without its fold; "
                    "incomplete branches use the independent tip profile"
                ),
            },
        )

    return EventReportSummary(
        tips=len(tips),
        matched_tips=matched_tip_count,
        incomplete_branches=len(incomplete_branches),
        matched_incomplete_branches=matched_incomplete_count,
        bubbles=len(bubbles),
        paths=path_count,
        projection_candidates=projection_candidates,
    )
