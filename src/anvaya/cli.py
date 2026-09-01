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
from anvaya.dead_end_attribution import (
    DeadEndAttributionSummary,
    write_dead_end_attribution_report,
)
from anvaya.damage_consensus import (
    DamageConsensusSummary,
    apply_damage_aware_consensus,
)
from anvaya.events import EventReportSummary, write_event_report
from anvaya.event_calibration import calibrate_event_report
from anvaya.fragmentation import FragmentationSummary, write_fragmentation_report
from anvaya.graph import build_dbg
from anvaya.incomplete_branches import find_incomplete_branch_candidates
from anvaya.metrics import summarize_graph
from anvaya.output import write_fasta
from anvaya.overlap_assembly import (
    OverlapCorrectionEvent,
    assemble_overlap_contigs,
    write_overlap_correction_report,
)
from anvaya.paired_extension import (
    PairedExtensionSummary,
    extend_paired_unitig_paths,
    spell_extended_paths,
)
from anvaya.reads import load_reads
from anvaya.read_correction import CorrectionSummary, write_correction_report
from anvaya.read_threading import (
    ReadThreadingSummary,
    audit_read_threads,
    resolve_read_thread_extensions,
    write_read_thread_report,
)
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
        "--read-thread-report",
        type=Path,
        help=(
            "audit direct read crossings of compacted-graph junctions as TSV "
            "without changing the assembly"
        ),
    )
    assemble_parser.add_argument(
        "--read-thread-extension",
        action="store_true",
        help=(
            "join only reciprocal strong-winner unitig transitions directly "
            "crossed by source reads"
        ),
    )
    assemble_parser.add_argument(
        "--thread-min-support",
        type=_minimum_count,
        default=5,
        help="minimum molecules for an auditable junction winner (default: 5)",
    )
    assemble_parser.add_argument(
        "--thread-dominance-ratio",
        type=_dominance_ratio,
        default=3.0,
        help="threaded winner-to-runner-up ratio (default: 3.0)",
    )
    assemble_parser.add_argument(
        "--thread-repeat-coverage-ratio",
        type=_dominance_ratio,
        default=2.0,
        help=(
            "break joined chains through unitigs with at least this local "
            "coverage multiple (default: 2.0)"
        ),
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
        "--fragmentation-report",
        type=Path,
        help=(
            "write report-only compacted-unitig end topology as TSV "
            "without changing the assembly"
        ),
    )
    assemble_parser.add_argument(
        "--dead-end-report",
        type=Path,
        help=(
            "attribute compacted-graph dead ends to raw read boundaries, "
            "filtered continuations, quality, or coverage as TSV"
        ),
    )
    assemble_parser.add_argument(
        "--correction-report",
        type=Path,
        help=(
            "audit conservative solid-k-mer substitution candidates as TSV "
            "without changing reads or assembly"
        ),
    )
    assemble_parser.add_argument(
        "--damage-consensus",
        action="store_true",
        help=(
            "correct strongly supported terminal deaminations from "
            "damage-tolerant read overlaps before graph construction"
        ),
    )
    assemble_parser.add_argument(
        "--damage-consensus-report",
        type=Path,
        help="write damage-aware consensus decisions as TSV",
    )
    assemble_parser.add_argument(
        "--correction-max-quality",
        type=_end_window,
        default=20,
        help="maximum Phred quality considered for correction (default: 20)",
    )
    assemble_parser.add_argument(
        "--correction-max-reads",
        type=_end_window,
        default=0,
        help=(
            "audit an evenly spaced read sample after building the full "
            "k-mer spectrum (default: 0, all reads)"
        ),
    )
    assemble_parser.add_argument(
        "--output",
        "-o",
        required=True,
        type=Path,
        help="output unitig FASTA file",
    )

    overlap_parser = subparsers.add_parser(
        "overlap-assemble",
        help="assemble merged fragments directly by conservative overlaps",
    )
    overlap_parser.add_argument(
        "--input",
        "-i",
        required=True,
        type=Path,
        help="merged-fragment FASTA/FASTQ input, optionally gzipped",
    )
    overlap_parser.add_argument(
        "--output", "-o", required=True, type=Path,
        help="output overlap-contig FASTA file",
    )
    overlap_parser.add_argument(
        "--max-rounds",
        type=_minimum_count,
        default=3,
        help="maximum extension rounds per seed (default: 3)",
    )
    overlap_parser.add_argument(
        "--max-contig-iterations",
        type=_end_window,
        default=3,
        help="maximum contig reindex-and-merge iterations (default: 3)",
    )
    overlap_parser.add_argument(
        "--min-cluster-size",
        type=_minimum_count,
        default=5,
        help="minimum independent fragments per contig (default: 5)",
    )
    overlap_parser.add_argument(
        "--anchor-k",
        type=_minimum_count,
        default=15,
        help="candidate-discovery anchor length (default: 15)",
    )
    overlap_parser.add_argument(
        "--anchors-per-read",
        type=_minimum_count,
        default=8,
        help="maximum sketched anchors per read (default: 8)",
    )
    overlap_parser.add_argument(
        "--max-anchor-occurrences",
        type=_minimum_count,
        default=100,
        help="ignore anchors occurring more often than this (default: 100)",
    )
    overlap_parser.add_argument(
        "--min-anchor-matches",
        type=_minimum_count,
        default=2,
        help="minimum agreeing anchors for an overlap candidate (default: 2)",
    )
    overlap_parser.add_argument(
        "--min-overlap",
        type=_minimum_count,
        default=30,
        help="minimum validated overlap length (default: 30)",
    )
    overlap_parser.add_argument(
        "--ranked-extension",
        action="store_true",
        help="extend supported clusters with a unique best candidate",
    )
    overlap_parser.add_argument(
        "--extension-consensus",
        action="store_true",
        help="recall supported bases after ranked extension",
    )
    overlap_parser.add_argument(
        "--min-ranked-extension-support",
        type=_minimum_count,
        default=1,
        help="minimum candidates agreeing on the first appended base (default: 1)",
    )
    overlap_parser.add_argument(
        "--reciprocal-best-extension",
        action="store_true",
        help="require selected read and contig extensions to point back uniquely",
    )
    overlap_parser.add_argument(
        "--min-output-length",
        type=_minimum_count,
        default=0,
        help="minimum emitted contig length (default: 0)",
    )
    overlap_parser.add_argument(
        "--damage-end-window",
        type=_end_window,
        default=0,
        help="terminal bases eligible for experimental damage polishing (default: 0)",
    )
    overlap_parser.add_argument(
        "--correction-report",
        type=Path,
        help="optional accepted-correction audit TSV",
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
    read_thread_report: Path | None,
    read_thread_extension: bool,
    thread_min_support: int,
    thread_dominance_ratio: float,
    thread_repeat_coverage_ratio: float,
    detect_bubbles: bool,
    end_window: int,
    event_report: Path | None,
    damage_profile_report: Path | None,
    unitig_bubble_report: Path | None,
    fragmentation_report: Path | None,
    dead_end_report: Path | None,
    correction_report: Path | None,
    correction_max_quality: int,
    correction_max_reads: int,
    damage_consensus: bool,
    damage_consensus_report: Path | None,
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
    if read_thread_report is not None and not orientation_aware:
        raise ValueError("--read-thread-report requires --orientation-aware")
    if read_thread_extension and not orientation_aware:
        raise ValueError("--read-thread-extension requires --orientation-aware")
    if read_thread_extension and end_window == 0:
        raise ValueError(
            "--read-thread-extension requires a positive --end-window"
        )
    if read_thread_extension and paired_unitig_extension:
        raise ValueError(
            "--read-thread-extension and --paired-unitig-extension are "
            "currently mutually exclusive"
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
    if fragmentation_report is not None and not orientation_aware:
        raise ValueError("--fragmentation-report requires --orientation-aware")
    if dead_end_report is not None and not orientation_aware:
        raise ValueError("--dead-end-report requires --orientation-aware")
    if correction_report is not None and min_count < 2:
        raise ValueError("--correction-report requires --min-count of at least 2")
    if damage_consensus and end_window == 0:
        raise ValueError("--damage-consensus requires a positive --end-window")
    if damage_consensus_report is not None and not damage_consensus:
        raise ValueError("--damage-consensus-report requires --damage-consensus")

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
    raw_sequences = [read.sequence for read in reads]
    sequences = raw_sequences
    _progress(f"Loaded {len(reads)} reads in {time.perf_counter() - stage_started:.2f}s")

    damage_consensus_summary = DamageConsensusSummary()
    if damage_consensus:
        stage_started = time.perf_counter()
        _progress("Correcting terminal deaminations from read consensus")
        consensus_reads, damage_consensus_summary = apply_damage_aware_consensus(
            reads,
            molecule_ids=molecule_ids,
            report_path=damage_consensus_report,
            end_window=end_window,
        )
        sequences = [read.sequence for read in consensus_reads]
        _progress(
            f"Corrected {damage_consensus_summary.corrected_bases} bases in "
            f"{damage_consensus_summary.corrected_reads} reads in "
            f"{time.perf_counter() - stage_started:.2f}s"
        )

    correction_summary = CorrectionSummary()
    if correction_report is not None:
        stage_started = time.perf_counter()
        _progress("Auditing conservative solid-k-mer correction candidates")
        correction_summary = write_correction_report(
            reads,
            correction_report,
            k=k,
            min_count=min_count,
            maximum_quality=correction_max_quality,
            end_window=max(end_window, 5),
            maximum_reads=correction_max_reads,
        )
        _progress(
            f"Reported {correction_summary.would_correct} correction "
            f"candidates to {correction_report} in "
            f"{time.perf_counter() - stage_started:.2f}s"
        )

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
            reads=raw_sequences,
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
            track_edge_handles=(
                paired_unitig_extension
                or read_thread_report is not None
                or read_thread_extension
            ),
        )
        unitigs = unitig_graph.sequences
    else:
        unitig_graph = None
        unitigs = extract_unitigs(graph)
    _progress(f"Extracted {len(unitigs)} unitigs in {time.perf_counter() - stage_started:.2f}s")

    fragmentation_summary = FragmentationSummary()
    if fragmentation_report is not None:
        stage_started = time.perf_counter()
        _progress("Attributing compacted-unitig boundaries")
        fragmentation_summary = write_fragmentation_report(
            unitig_graph, fragmentation_report
        )
        _progress(
            f"Reported {fragmentation_summary.ends} unitig ends to "
            f"{fragmentation_report} in "
            f"{time.perf_counter() - stage_started:.2f}s"
        )

    dead_end_summary = DeadEndAttributionSummary()
    if dead_end_report is not None:
        stage_started = time.perf_counter()
        _progress("Attributing dead ends from raw read evidence")
        dead_end_summary = write_dead_end_attribution_report(
            unitig_graph,
            reads,
            dead_end_report,
            min_count=min_count,
            end_window=max(end_window, 5),
        )
        _progress(
            f"Attributed {dead_end_summary.dead_ends} dead ends to "
            f"{dead_end_report} in {time.perf_counter() - stage_started:.2f}s"
        )

    read_thread_summary = ReadThreadingSummary()
    threaded_links = 0
    if read_thread_report is not None or read_thread_extension:
        stage_started = time.perf_counter()
        _progress("Auditing direct read paths across unitig junctions")
        molecules = (
            zip(*read_groups, strict=True)
            if len(read_groups) == 2
            else ((read,) for read in read_groups[0])
        )
        transitions, read_thread_summary = audit_read_threads(
            graph,
            unitig_graph,
            molecules,
            end_window=end_window,
            minimum_support=thread_min_support,
            dominance_ratio=thread_dominance_ratio,
        )
        if read_thread_report is not None:
            write_read_thread_report(transitions, read_thread_report)
        if read_thread_extension:
            successors = resolve_read_thread_extensions(
                transitions,
                repeat_coverage_ratio=thread_repeat_coverage_ratio,
            )
            threaded_links = len(successors) // 2
            unitigs = list(spell_extended_paths(unitig_graph, successors))
        _progress(
            f"Observed {read_thread_summary.supported_transitions} direct "
            f"transitions and {read_thread_summary.resolvable_junctions} "
            f"strong-winner junctions ({read_thread_summary.reciprocal_resolvable_links} "
            f"reciprocal physical links) in "
            f"{time.perf_counter() - stage_started:.2f}s"
        )

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
    print(f"projection_candidates={event_summary.projection_candidates}")
    print(f"event_report={event_report or ''}")
    print(f"damage_profile_report={damage_profile_report or ''}")
    print(f"unitig_bubbles_detected={unitig_bubble_summary.bubbles}")
    print(f"unitig_bubble_paths={unitig_bubble_summary.paths}")
    print(f"unitig_error_like={unitig_bubble_summary.error_like}")
    print(f"unitig_damage_like={unitig_bubble_summary.damage_like}")
    print(f"unitig_variation_like={unitig_bubble_summary.variation_like}")
    print(f"unitig_ambiguous={unitig_bubble_summary.ambiguous}")
    print(f"unitig_bubble_report={unitig_bubble_report or ''}")
    print(f"fragmentation_report={fragmentation_report or ''}")
    print(f"fragmentation_ends={fragmentation_summary.ends}")
    print(f"fragmentation_dead_ends={fragmentation_summary.dead_ends}")
    print(
        "fragmentation_unique_continuations="
        f"{fragmentation_summary.unique_continuations}"
    )
    print(
        "fragmentation_ambiguous_branches="
        f"{fragmentation_summary.ambiguous_branches}"
    )
    print(f"fragmentation_isolated_unitigs={fragmentation_summary.isolated_unitigs}")
    print(f"fragmentation_one_sided_unitigs={fragmentation_summary.one_sided_unitigs}")
    print(f"fragmentation_connected_unitigs={fragmentation_summary.connected_unitigs}")
    print(f"dead_end_report={dead_end_report or ''}")
    print(f"dead_end_total={dead_end_summary.dead_ends}")
    print(f"dead_end_read_boundaries={dead_end_summary.read_boundaries}")
    print(f"dead_end_coverage_gaps={dead_end_summary.coverage_gaps}")
    print(f"dead_end_filtered_unique={dead_end_summary.filtered_unique}")
    print(f"dead_end_filtered_conflicts={dead_end_summary.filtered_conflicts}")
    print(
        "dead_end_retained_context_elsewhere="
        f"{dead_end_summary.retained_context_elsewhere}"
    )
    print(f"dead_end_low_quality={dead_end_summary.low_quality}")
    print(f"dead_end_terminal_only={dead_end_summary.terminal_only}")
    print(f"dead_end_missing_quality={dead_end_summary.missing_quality}")
    print(f"correction_report={correction_report or ''}")
    print(f"correction_total_reads={correction_summary.reads}")
    print(f"correction_audited_reads={correction_summary.audited_reads}")
    print(f"correction_low_quality_bases={correction_summary.low_quality_bases}")
    print(f"correction_would_correct={correction_summary.would_correct}")
    print(f"correction_protected_damage={correction_summary.protected_damage}")
    print(f"correction_ambiguous={correction_summary.ambiguous}")
    print(f"correction_incomplete_rescue={correction_summary.incomplete_rescue}")
    print(f"damage_consensus={str(damage_consensus).lower()}")
    print(f"damage_consensus_report={damage_consensus_report or ''}")
    print(f"damage_consensus_candidates={damage_consensus_summary.candidate_bases}")
    print(f"damage_consensus_corrected_bases={damage_consensus_summary.corrected_bases}")
    print(f"damage_consensus_corrected_reads={damage_consensus_summary.corrected_reads}")
    print(
        "damage_consensus_insufficient_support="
        f"{damage_consensus_summary.insufficient_support}"
    )
    print(f"damage_consensus_ambiguous={damage_consensus_summary.ambiguous}")
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
    print(f"read_thread_report={read_thread_report or ''}")
    print(f"read_thread_extension={str(read_thread_extension).lower()}")
    print(f"thread_joined_unitig_links={threaded_links}")
    print(f"thread_junctions={read_thread_summary.junctions}")
    print(
        "thread_candidate_transitions="
        f"{read_thread_summary.candidate_transitions}"
    )
    print(
        "thread_supported_transitions="
        f"{read_thread_summary.supported_transitions}"
    )
    print(f"threaded_molecules={read_thread_summary.threaded_molecules}")
    print(f"threaded_reads={read_thread_summary.threaded_reads}")
    print(
        "thread_resolvable_junctions="
        f"{read_thread_summary.resolvable_junctions}"
    )
    print(
        "thread_reciprocal_resolvable_links="
        f"{read_thread_summary.reciprocal_resolvable_links}"
    )
    print(f"thread_repeat_coverage_ratio={thread_repeat_coverage_ratio}")
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
                arguments.read_thread_report,
                arguments.read_thread_extension,
                arguments.thread_min_support,
                arguments.thread_dominance_ratio,
                arguments.thread_repeat_coverage_ratio,
                arguments.detect_bubbles,
                arguments.end_window,
                arguments.event_report,
                arguments.damage_profile_report,
                arguments.unitig_bubble_report,
                arguments.fragmentation_report,
                arguments.dead_end_report,
                arguments.correction_report,
                arguments.correction_max_quality,
                arguments.correction_max_reads,
                arguments.damage_consensus,
                arguments.damage_consensus_report,
            )
        if arguments.command == "overlap-assemble":
            started = time.perf_counter()
            _progress(f"Loading merged fragments from {arguments.input}")
            reads = load_reads(arguments.input)
            _progress(f"Loaded {len(reads)} fragments")
            correction_events: list[OverlapCorrectionEvent] = []
            contigs, summary = assemble_overlap_contigs(
                reads,
                anchor_k=arguments.anchor_k,
                anchors_per_read=arguments.anchors_per_read,
                maximum_anchor_occurrences=arguments.max_anchor_occurrences,
                minimum_anchor_matches=arguments.min_anchor_matches,
                minimum_overlap=arguments.min_overlap,
                ranked_extension=arguments.ranked_extension,
                extension_consensus=arguments.extension_consensus,
                minimum_ranked_extension_support=arguments.min_ranked_extension_support,
                reciprocal_best_extension=arguments.reciprocal_best_extension,
                maximum_rounds=arguments.max_rounds,
                maximum_contig_iterations=arguments.max_contig_iterations,
                minimum_cluster_size=arguments.min_cluster_size,
                minimum_output_length=arguments.min_output_length,
                damage_end_window=arguments.damage_end_window,
                correction_events=correction_events,
            )
            if arguments.correction_report is not None:
                write_overlap_correction_report(
                    correction_events,
                    arguments.correction_report,
                )
            write_fasta([contig.sequence for contig in contigs], arguments.output)
            _progress(f"Completed in {time.perf_counter() - started:.2f}s")
            print(f"reads={summary.input_reads}")
            print(f"overlap_contigs={summary.output_contigs}")
            print(f"overlap_clusters={summary.clusters}")
            print(f"overlap_clustered_reads={summary.clustered_reads}")
            print(f"overlap_extension_rounds={summary.extension_rounds}")
            print(f"overlap_extended_contigs={summary.extended_contigs}")
            print(f"overlap_added_bases={summary.added_bases}")
            print(f"overlap_contig_iterations={summary.contig_iterations}")
            print(f"overlap_contig_merges={summary.contig_merges}")
            print(f"overlap_corrected_bases={summary.corrected_bases}")
            print(f"overlap_correction_report={arguments.correction_report or ''}")
            print(
                "overlap_correction_candidates="
                f"{summary.correction_candidates}"
            )
            print(
                "overlap_correction_terminal_only="
                f"{summary.correction_terminal_only}"
            )
            print(
                "overlap_correction_insufficient_support="
                f"{summary.correction_insufficient_support}"
            )
            print(
                "overlap_correction_ambiguous="
                f"{summary.correction_ambiguous}"
            )
            print(
                "overlap_ambiguous_extensions="
                f"{summary.ambiguous_extensions}"
            )
            print(f"overlap_candidate_offsets={summary.candidate_offsets}")
            print(f"overlap_candidate_below_anchor_support={summary.candidate_below_anchor_support}")
            print(f"overlap_candidate_short_overlap={summary.candidate_short_overlap}")
            print(f"overlap_candidate_dna_rejected={summary.candidate_dna_rejected}")
            print(f"overlap_candidate_ry_rejected={summary.candidate_ry_rejected}")
            print(f"overlap_candidate_molecule_ambiguous={summary.candidate_molecule_ambiguous}")
            print(f"overlap_candidate_alignments={summary.candidate_alignments}")
            print(f"overlap_unavailable_anchor_hits={summary.unavailable_anchor_hits}")
            print(f"overlap_targets_without_candidates={summary.targets_without_candidates}")
            print(f"overlap_clusters_below_minimum_size={summary.clusters_below_minimum_size}")
            print(f"overlap_consensus_without_extension={summary.consensus_without_extension}")
            print(f"overlap_reciprocal_extension_checks={summary.reciprocal_extension_checks}")
            print(
                "overlap_reciprocal_extension_rejections="
                f"{summary.reciprocal_extension_rejections}"
            )
            print(f"output={arguments.output}")
            return 0
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
