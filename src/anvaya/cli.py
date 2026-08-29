"""Command-line interface for Anvaya."""

import argparse
import sys
import time
from collections.abc import Sequence
from pathlib import Path

from anvaya.bidirected import (
    build_bidirected_dbg,
    summarize_bidirected_graph,
)
from anvaya.bubbles import find_simple_bubbles
from anvaya.cleaning import (
    TipCleaningSummary,
    remove_damage_aware_tips,
    remove_weak_tips,
)
from anvaya.cleaning import find_weak_tip_candidates
from anvaya.events import EventReportSummary, write_event_report
from anvaya.event_calibration import calibrate_event_report
from anvaya.graph import build_dbg
from anvaya.incomplete_branches import find_incomplete_branch_candidates
from anvaya.metrics import summarize_graph
from anvaya.output import write_fasta
from anvaya.paired_extension import (
    PairedExtensionSummary,
    extend_paired_unitig_paths,
)
from anvaya.reads import load_reads
from anvaya.unitigs import extract_unitigs
from anvaya.unitig_bubbles import (
    UnitigBubbleReportSummary,
    find_unitig_bubbles,
    write_unitig_bubble_report,
)
from anvaya.unitig_graph import build_compacted_unitig_graph


def _kmer_size(value: str) -> int:
    try:
        k = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("k must be an integer") from error
    if k < 2:
        raise argparse.ArgumentTypeError("k must be at least 2")
    return k


def _minimum_count(value: str) -> int:
    try:
        count = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("minimum count must be an integer") from error
    if count < 1:
        raise argparse.ArgumentTypeError("minimum count must be at least 1")
    return count


def _end_window(value: str) -> int:
    try:
        window = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "end window must be an integer"
        ) from error
    if window < 0:
        raise argparse.ArgumentTypeError("end window must not be negative")
    return window


def _probability(value: str) -> float:
    try:
        probability = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("probability must be numeric") from error
    if not 0.0 < probability < 1.0:
        raise argparse.ArgumentTypeError("probability must be between zero and one")
    return probability


def _dominance_ratio(value: str) -> float:
    try:
        ratio = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "dominance ratio must be numeric"
        ) from error
    if ratio < 1.0:
        raise argparse.ArgumentTypeError(
            "dominance ratio must be at least 1"
        )
    return ratio


def _progress(message: str) -> None:
    print(f"[anvaya] {message}", file=sys.stderr, flush=True)


def build_parser() -> argparse.ArgumentParser:
    """Create the Anvaya argument parser."""
    parser = argparse.ArgumentParser(
        prog="anvaya",
        description="Damage-aware de Bruijn graph assembly research prototype",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    assemble_parser = subparsers.add_parser(
        "assemble",
        help="assemble FASTA/FASTQ reads into unitigs",
    )
    input_group = assemble_parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input",
        "-i",
        type=Path,
        help="single-end FASTA/FASTQ input, optionally gzipped",
    )
    input_group.add_argument(
        "--left",
        "-1",
        type=Path,
        help="left paired-end FASTA/FASTQ input, optionally gzipped",
    )
    assemble_parser.add_argument(
        "--right",
        "-2",
        type=Path,
        help="right paired-end FASTA/FASTQ input, optionally gzipped",
    )
    assemble_parser.add_argument(
        "--k",
        required=True,
        type=_kmer_size,
        help="k-mer size (at least 2)",
    )
    assemble_parser.add_argument(
        "--min-count",
        type=_minimum_count,
        default=1,
        help="minimum k-mer support to retain (default: 1)",
    )
    assemble_parser.add_argument(
        "--orientation-aware",
        action="store_true",
        help="use the experimental orientation-aware graph",
    )
    cleaning_group = assemble_parser.add_mutually_exclusive_group()
    cleaning_group.add_argument(
        "--clean-tips",
        action="store_true",
        help="experimentally remove short, weak dead-end paths",
    )
    cleaning_group.add_argument(
        "--damage-aware-clean-tips",
        action="store_true",
        help=(
            "iteratively remove only one-sided error-like weak tips and "
            "protect damage-like, variation-like, ambiguous, unmatched, "
            "or bidirectionally supported tips"
        ),
    )
    assemble_parser.add_argument(
        "--detect-bubbles",
        action="store_true",
        help="report simple graph bubbles without changing the assembly",
    )
    assemble_parser.add_argument(
        "--paired-unitig-extension",
        action="store_true",
        help=(
            "extend compacted unitigs only across reciprocal paired-read "
            "strong-winner links"
        ),
    )
    assemble_parser.add_argument(
        "--paired-min-support",
        type=_minimum_count,
        default=5,
        help="minimum independent read pairs for extension (default: 5)",
    )
    assemble_parser.add_argument(
        "--paired-dominance-ratio",
        type=_dominance_ratio,
        default=3.0,
        help="winner-to-runner-up support ratio (default: 3.0)",
    )
    assemble_parser.add_argument(
        "--paired-max-distance",
        type=_minimum_count,
        default=1000,
        help="maximum graph distance for mate support (default: 1000)",
    )
    assemble_parser.add_argument(
        "--paired-max-search-states",
        type=_minimum_count,
        default=1_000,
        help="maximum graph states searched per junction (default: 1000)",
    )
    assemble_parser.add_argument(
        "--end-window",
        type=_end_window,
        default=0,
        help=(
            "collect canonicalized read-end evidence within this many "
            "k-mer positions (default: 0, disabled)"
        ),
    )
    assemble_parser.add_argument(
        "--event-report",
        type=Path,
        help="write pre-cleaning tip and bubble evidence as TSV",
    )
    assemble_parser.add_argument(
        "--damage-profile-report",
        type=Path,
        help=(
            "write a candidate-locus terminal damage profile as JSON "
            "(requires --event-report)"
        ),
    )
    assemble_parser.add_argument(
        "--unitig-bubble-report",
        type=Path,
        help="score compacted-graph bubble paths and write them as TSV",
    )
    assemble_parser.add_argument(
        "--output",
        "-o",
        required=True,
        type=Path,
        help="output unitig FASTA file",
    )

    calibration_parser = subparsers.add_parser(
        "calibrate-events",
        help="append report-only conformal confidence to an event TSV",
    )
    calibration_parser.add_argument("--input", "-i", required=True, type=Path)
    calibration_parser.add_argument("--model", required=True, type=Path)
    calibration_parser.add_argument("--output", "-o", required=True, type=Path)
    calibration_parser.add_argument(
        "--alpha",
        type=_probability,
        default=0.01,
        help="class rejection level before BH correction (default: 0.01)",
    )
    calibration_parser.add_argument(
        "--min-observations",
        type=_minimum_count,
        default=5,
        help="minimum molecules required for a decision (default: 5)",
    )
    calibration_parser.add_argument(
        "--min-alternative-observations",
        type=_minimum_count,
        default=2,
        help="minimum alternative molecules required for a decision (default: 2)",
    )
    calibration_parser.add_argument(
        "--min-reference-observations",
        type=_minimum_count,
        default=5,
        help="minimum reference molecules required for a decision (default: 5)",
    )
    return parser


def _resolve_input_paths(
    parser: argparse.ArgumentParser, arguments: argparse.Namespace
) -> list[Path]:
    if arguments.input is not None:
        if arguments.right is not None:
            parser.error("--right/-2 can only be used with --left/-1")
        return [arguments.input]

    if arguments.right is None:
        parser.error("paired-end input requires both --left/-1 and --right/-2")
    return [arguments.left, arguments.right]


def _run_assemble(
    input_paths: list[Path],
    output_path: Path,
    k: int,
    min_count: int,
    orientation_aware: bool,
    clean_tips: bool,
    damage_aware_clean_tips: bool,
    paired_unitig_extension: bool,
    paired_min_support: int,
    paired_dominance_ratio: float,
    paired_max_distance: int,
    paired_max_search_states: int,
    detect_bubbles: bool,
    end_window: int,
    event_report: Path | None,
    damage_profile_report: Path | None,
    unitig_bubble_report: Path | None,
) -> int:
    started = time.perf_counter()
    if clean_tips and not orientation_aware:
        raise ValueError("--clean-tips requires --orientation-aware")
    if damage_aware_clean_tips and not orientation_aware:
        raise ValueError(
            "--damage-aware-clean-tips requires --orientation-aware"
        )
    if damage_aware_clean_tips and end_window == 0:
        raise ValueError(
            "--damage-aware-clean-tips requires a positive --end-window"
        )
    if paired_unitig_extension and len(input_paths) != 2:
        raise ValueError("--paired-unitig-extension requires paired-end input")
    if paired_unitig_extension and not orientation_aware:
        raise ValueError(
            "--paired-unitig-extension requires --orientation-aware"
        )
    if paired_unitig_extension and end_window == 0:
        raise ValueError(
            "--paired-unitig-extension requires a positive --end-window"
        )
    if detect_bubbles and not orientation_aware:
        raise ValueError("--detect-bubbles requires --orientation-aware")
    if end_window and not orientation_aware:
        raise ValueError("--end-window requires --orientation-aware")
    if event_report is not None and end_window == 0:
        raise ValueError("--event-report requires a positive --end-window")
    if damage_profile_report is not None and event_report is None:
        raise ValueError("--damage-profile-report requires --event-report")
    if unitig_bubble_report is not None and not orientation_aware:
        raise ValueError(
            "--unitig-bubble-report requires --orientation-aware"
        )

    stage_started = time.perf_counter()
    _progress(f"Loading reads from {', '.join(map(str, input_paths))}")
    read_groups = [load_reads(path) for path in input_paths]
    if len(read_groups) == 2 and len(read_groups[0]) != len(read_groups[1]):
        raise ValueError("paired-end input files must contain the same number of reads")
    reads = [read for group in read_groups for read in group]
    molecule_ids = (
        list(range(len(read_groups[0]))) * 2
        if len(read_groups) == 2
        else list(range(len(reads)))
    )
    sequences = [read.sequence for read in reads]
    _progress(f"Loaded {len(reads)} reads in {time.perf_counter() - stage_started:.2f}s")

    stage_started = time.perf_counter()
    graph_kind = "orientation-aware" if orientation_aware else "directed"
    _progress(
        f"Building {graph_kind} de Bruijn graph with "
        f"k={k}, min_count={min_count}"
    )
    if orientation_aware:
        graph = build_bidirected_dbg(
            sequences,
            k,
            min_count,
            end_window=end_window,
            track_molecule_links=(
                event_report is not None or paired_unitig_extension
            ),
            read_qualities=[read.qualities for read in reads],
            read_molecule_ids=molecule_ids,
        )
        node_count = graph.node_count
        edge_count = graph.edge_count
    else:
        graph = build_dbg(sequences, k, min_count)
        node_count = len(graph)
        edge_count = sum(len(successors) for successors in graph.values())
    _progress(
        f"Built graph with {node_count} nodes and {edge_count} edges "
        f"in {time.perf_counter() - stage_started:.2f}s"
    )

    event_summary = EventReportSummary()
    if event_report is not None:
        stage_started = time.perf_counter()
        _progress("Reporting pre-cleaning graph events")
        tip_candidates = find_weak_tip_candidates(graph)
        incomplete_branch_candidates = find_incomplete_branch_candidates(graph)
        bubble_candidates = find_simple_bubbles(graph)
        event_summary = write_event_report(
            graph,
            tip_candidates,
            bubble_candidates,
            event_report,
            incomplete_branches=incomplete_branch_candidates,
            damage_profile_path=damage_profile_report,
            reads=sequences,
            read_qualities=[read.qualities for read in reads],
            read_molecule_ids=molecule_ids,
        )
        _progress(
            f"Reported {event_summary.tips} tips "
            f"({event_summary.matched_tips} matched to backbones) and "
            f"{event_summary.incomplete_branches} incomplete branches "
            f"({event_summary.matched_incomplete_branches} matched) and "
            f"{event_summary.bubbles} bubbles to {event_report} in "
            f"{time.perf_counter() - stage_started:.2f}s"
        )

    cleaning_summary = TipCleaningSummary()
    if clean_tips or damage_aware_clean_tips:
        stage_started = time.perf_counter()
        if damage_aware_clean_tips:
            _progress(
                "Removing one-sided error-like tips with damage-aware "
                "protection"
            )
            cleaning_summary = remove_damage_aware_tips(graph)
        else:
            _progress("Removing short, weak tips")
            cleaning_summary = remove_weak_tips(graph)
        _progress(
            f"Removed {cleaning_summary.tips_removed} tips containing "
            f"{cleaning_summary.edges_removed} edges in "
            f"{time.perf_counter() - stage_started:.2f}s"
        )

    bubbles = []
    if detect_bubbles:
        stage_started = time.perf_counter()
        _progress("Detecting simple bubbles")
        bubbles = find_simple_bubbles(graph)
        _progress(
            f"Detected {len(bubbles)} simple bubbles in "
            f"{time.perf_counter() - stage_started:.2f}s"
        )

    stage_started = time.perf_counter()
    _progress("Extracting unitigs")
    if orientation_aware:
        unitig_graph = build_compacted_unitig_graph(
            graph,
            track_edge_handles=paired_unitig_extension,
        )
        unitigs = unitig_graph.sequences
    else:
        unitig_graph = None
        unitigs = extract_unitigs(graph)
    _progress(f"Extracted {len(unitigs)} unitigs in {time.perf_counter() - stage_started:.2f}s")

    paired_extension_summary = PairedExtensionSummary()
    if paired_unitig_extension:
        stage_started = time.perf_counter()
        _progress("Extending unitigs with reciprocal paired-read evidence")
        paired_result = extend_paired_unitig_paths(
            graph,
            unitig_graph,
            minimum_support=paired_min_support,
            dominance_ratio=paired_dominance_ratio,
            maximum_distance=paired_max_distance,
            maximum_states=paired_max_search_states,
        )
        unitigs = list(paired_result.sequences)
        paired_extension_summary = paired_result.summary
        _progress(
            f"Joined {paired_extension_summary.joined_unitig_links} "
            f"unitig links into {len(unitigs)} contigs in "
            f"{time.perf_counter() - stage_started:.2f}s"
        )

    unitig_bubble_summary = UnitigBubbleReportSummary()
    if unitig_bubble_report is not None:
        stage_started = time.perf_counter()
        _progress("Scoring and classifying compacted-graph bubbles")
        unitig_bubbles = find_unitig_bubbles(unitig_graph)
        unitig_bubble_summary = write_unitig_bubble_report(
            unitig_graph,
            unitig_bubbles,
            unitig_bubble_report,
        )
        _progress(
            f"Reported {unitig_bubble_summary.bubbles} compacted-graph "
            f"bubbles to {unitig_bubble_report} in "
            f"{time.perf_counter() - stage_started:.2f}s"
        )

    stage_started = time.perf_counter()
    _progress("Calculating graph statistics")
    if orientation_aware:
        summary = summarize_bidirected_graph(graph, unitigs)
    else:
        summary = summarize_graph(graph, unitigs)
    _progress(f"Calculated graph statistics in {time.perf_counter() - stage_started:.2f}s")

    stage_started = time.perf_counter()
    _progress(f"Writing unitigs to {output_path}")
    write_fasta(unitigs, output_path)
    _progress(f"Wrote output in {time.perf_counter() - stage_started:.2f}s")
    _progress(f"Completed in {time.perf_counter() - started:.2f}s")

    print(f"input_files={len(input_paths)}")
    print(f"reads={len(reads)}")
    print(f"k={k}")
    print(f"min_count={min_count}")
    print(f"orientation_aware={str(orientation_aware).lower()}")
    print(
        "tip_cleaning="
        f"{str(clean_tips or damage_aware_clean_tips).lower()}"
    )
    print(
        "damage_aware_tip_cleaning="
        f"{str(damage_aware_clean_tips).lower()}"
    )
    print(f"tips_removed={cleaning_summary.tips_removed}")
    print(f"tip_edges_removed={cleaning_summary.edges_removed}")
    print(f"tip_observations_removed={cleaning_summary.observations_removed}")
    print(f"tip_cleaning_rounds={cleaning_summary.rounds}")
    print(f"tips_protected={cleaning_summary.tips_protected}")
    print(
        "damage_like_tips_protected="
        f"{cleaning_summary.damage_like_protected}"
    )
    print(
        "variation_like_tips_protected="
        f"{cleaning_summary.variation_like_protected}"
    )
    print(
        "ambiguous_tips_protected="
        f"{cleaning_summary.ambiguous_protected}"
    )
    print(
        "unmatched_tips_protected="
        f"{cleaning_summary.unmatched_protected}"
    )
    print(
        "bidirectional_tips_protected="
        f"{cleaning_summary.bidirectional_protected}"
    )
    print(f"bubble_detection={str(detect_bubbles).lower()}")
    print(f"bubbles_detected={len(bubbles)}")
    print(f"end_window={end_window}")
    print(
        "terminal_observations="
        f"{graph.terminal_observations if orientation_aware else 0}"
    )
    print(f"reported_tips={event_summary.tips}")
    print(f"reported_tip_matches={event_summary.matched_tips}")
    print(f"reported_incomplete_branches={event_summary.incomplete_branches}")
    print(
        "reported_incomplete_branch_matches="
        f"{event_summary.matched_incomplete_branches}"
    )
    print(f"reported_bubbles={event_summary.bubbles}")
    print(f"reported_paths={event_summary.paths}")
    print(f"event_report={event_report or ''}")
    print(f"damage_profile_report={damage_profile_report or ''}")
    print(f"unitig_bubbles_detected={unitig_bubble_summary.bubbles}")
    print(f"unitig_bubble_paths={unitig_bubble_summary.paths}")
    print(f"unitig_error_like={unitig_bubble_summary.error_like}")
    print(f"unitig_damage_like={unitig_bubble_summary.damage_like}")
    print(f"unitig_variation_like={unitig_bubble_summary.variation_like}")
    print(f"unitig_ambiguous={unitig_bubble_summary.ambiguous}")
    print(f"unitig_bubble_report={unitig_bubble_report or ''}")
    print(
        "paired_unitig_extension="
        f"{str(paired_unitig_extension).lower()}"
    )
    print(
        "paired_molecules="
        f"{paired_extension_summary.paired_molecules}"
    )
    print(f"paired_mapped_pairs={paired_extension_summary.mapped_pairs}")
    print(
        "paired_evidence_links="
        f"{paired_extension_summary.evidence_links}"
    )
    print(
        "paired_resolved_oriented_junctions="
        f"{paired_extension_summary.resolved_oriented_junctions}"
    )
    print(
        "paired_evaluated_oriented_junctions="
        f"{paired_extension_summary.evaluated_oriented_junctions}"
    )
    print(
        "paired_work_limited_oriented_junctions="
        f"{paired_extension_summary.work_limited_oriented_junctions}"
    )
    print(
        "paired_searched_states="
        f"{paired_extension_summary.searched_states}"
    )
    print(
        "paired_joined_unitig_links="
        f"{paired_extension_summary.joined_unitig_links}"
    )
    print(f"nodes={summary.nodes}")
    print(f"edges={summary.edges}")
    print(f"observations={summary.observations}")
    print(f"branching_nodes={summary.branching_nodes}")
    print(f"unitigs={summary.unitigs}")
    print(
        "unitig_links="
        f"{unitig_graph.oriented_link_count if unitig_graph is not None else 0}"
    )
    print(f"output={output_path}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Run the Anvaya command-line interface."""
    parser = build_parser()
    arguments = parser.parse_args(argv)

    try:
        if arguments.command == "assemble":
            input_paths = _resolve_input_paths(parser, arguments)
            return _run_assemble(
                input_paths,
                arguments.output,
                arguments.k,
                arguments.min_count,
                arguments.orientation_aware,
                arguments.clean_tips,
                arguments.damage_aware_clean_tips,
                arguments.paired_unitig_extension,
                arguments.paired_min_support,
                arguments.paired_dominance_ratio,
                arguments.paired_max_distance,
                arguments.paired_max_search_states,
                arguments.detect_bubbles,
                arguments.end_window,
                arguments.event_report,
                arguments.damage_profile_report,
                arguments.unitig_bubble_report,
            )
        if arguments.command == "calibrate-events":
            summary = calibrate_event_report(
                arguments.input,
                arguments.model,
                arguments.output,
                alpha=arguments.alpha,
                minimum_observations=arguments.min_observations,
                minimum_alternative_observations=(
                    arguments.min_alternative_observations
                ),
                minimum_reference_observations=(
                    arguments.min_reference_observations
                ),
            )
            print(f"rows={summary.rows}")
            print(f"scored={summary.scored}")
            print(f"eligible_error={summary.eligible_error}")
            print(f"protect_damage={summary.protect_damage}")
            print(f"protect_variation={summary.protect_variation}")
            print(f"insufficient={summary.insufficient}")
            print(f"output={arguments.output}")
            return 0
    except (OSError, ValueError) as error:
        parser.error(str(error))

    parser.error(f"unknown command: {arguments.command}")
