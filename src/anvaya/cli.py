"""Command-line interface for Anvaya."""

import argparse
import json
import sys
import time
from contextlib import ExitStack
from dataclasses import asdict
from collections.abc import Sequence
from pathlib import Path

from anvaya.output import write_fasta, write_named_fasta
from anvaya.raw_consensus import DamageProfile, project_raw_consensus
from anvaya.overlap_assembly import (
    IterativeReclusteringDiagnostics,
    MasterOverlapGraphDiagnostics,
    OverlapCorrectionEvent,
    RawConfirmedMasterGraphDiagnostics,
    StrainSafeContainmentDiagnostics,
    TwoTierRedundancyDiagnostics,
    assemble_overlap_contigs,
    audit_two_tier_redundancy,
    write_overlap_correction_report,
)
from anvaya.overlap_graph import (
    audit_master_overlap_graph,
    audit_raw_confirmed_master_overlap_graph,
    audit_strain_safe_containment,
)
from anvaya.overlap_reclustering import (
    audit_iterative_reclustering,
    write_iterative_reclustering_report,
)
from anvaya.overlap_progressive import (
    ProgressiveClusteringDiagnostics,
    ProgressiveIterationDiagnostics,
    ProgressiveRawExtensionDiagnostics,
    ProgressiveSequencePool,
    SequenceState,
    discover_progressive_raw_clusters,
    extend_progressive_raw_clusters,
    iterate_progressive_raw_extension,
)
from anvaya.overlap_adaptive_rescue import (
    AdaptiveRescueDiagnostics,
    SelectiveRescueDiagnostics,
    SupportTwoRescueDiagnostics,
    audit_adaptive_support_rescue,
    project_high_confidence_support_two_rescue,
    project_selective_support_rescue,
)
from anvaya.overlap_progressive_links import (
    ProgressiveLinkDiagnostics,
    audit_raw_supported_progressive_links,
)
from anvaya.overlap_scaffolding import (
    PairedScaffoldDiagnostics,
    scaffold_progressive_contigs,
)
from anvaya.paired_reads import PairedMergeDiagnostics, merge_overlapping_pairs
from anvaya.reads import load_reads
from anvaya.damage_string_graph import assemble as assemble_damage_string_graph
from anvaya.carpedeam_backend import run_carpedeam_safe


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
        description="Damage-aware overlap assembly research prototype",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    damage_parser = subparsers.add_parser(
        "damage-assemble",
        help="assemble merged ancient-DNA fragments with the damage-aware string graph",
    )
    damage_parser.add_argument("--input", "-i", required=True, type=Path)
    damage_parser.add_argument("--profile-prefix", required=True, type=Path)
    damage_parser.add_argument("--output-before", required=True, type=Path)
    damage_parser.add_argument("--output-consensus", required=True, type=Path)
    damage_parser.add_argument(
        "--output-assembled", type=Path,
        help="optional consensus FASTA with at least two uniquely placed molecules",
    )
    damage_parser.add_argument(
        "--output-unresolved", type=Path,
        help="optional consensus FASTA with fewer than two uniquely placed molecules",
    )
    damage_parser.add_argument("--diagnostics", required=True, type=Path)
    damage_parser.add_argument("--placements-report", type=Path)
    damage_parser.add_argument("--consensus-report", type=Path)
    damage_parser.add_argument(
        "--max-reads", type=_minimum_count, default=1000,
        help="maximum input fragments for this experimental command (default: 1000)",
    )

    carpedeam_parser = subparsers.add_parser(
        "carpedeam-assemble",
        help="run CarpeDeam safe mode as a pinned external benchmark comparator",
    )
    carpedeam_parser.add_argument("--input", "-i", required=True, type=Path)
    carpedeam_parser.add_argument("--profile-prefix", required=True, type=Path)
    carpedeam_parser.add_argument("--output", "-o", required=True, type=Path)
    carpedeam_parser.add_argument("--temporary-directory", required=True, type=Path)
    carpedeam_parser.add_argument("--diagnostics", required=True, type=Path)
    carpedeam_parser.add_argument("--executable", default="carpedeam")
    carpedeam_parser.add_argument("--threads", type=_minimum_count, default=2)
    carpedeam_parser.add_argument(
        "--min-contig-length", type=_minimum_count, default=31,
    )
    damage_parser.add_argument(
        "--min-output-length", type=_minimum_count, default=31,
        help="minimum emitted contig length (default: 31)",
    )

    overlap_parser = subparsers.add_parser(
        "overlap-assemble",
        help="assemble fragments directly by conservative overlaps",
    )
    overlap_input_group = overlap_parser.add_mutually_exclusive_group(required=True)
    overlap_input_group.add_argument(
        "--input",
        "-i",
        type=Path,
        help="merged-fragment FASTA/FASTQ input, optionally gzipped",
    )
    overlap_input_group.add_argument(
        "--left",
        "-1",
        type=Path,
        help="left paired-end FASTA/FASTQ input, optionally gzipped",
    )
    overlap_parser.add_argument(
        "--right",
        "-2",
        type=Path,
        help="right paired-end FASTA/FASTQ input, optionally gzipped",
    )
    overlap_parser.add_argument(
        "--merge-overlapping-pairs",
        action="store_true",
        help="merge uniquely supported R1/R2 overlaps before assembly",
    )
    overlap_parser.add_argument(
        "--paired-merge-min-overlap",
        type=_minimum_count,
        default=20,
        help="minimum R1/R2 overlap for pair merging (default: 20)",
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
        "--damage-aware-ranking",
        action="store_true",
        help="rank extensions by damage-aware overlap confidence",
    )
    overlap_parser.add_argument(
        "--min-overlap-confidence-margin",
        type=float,
        default=0.0,
        help="minimum confidence gap between competing extensions (default: 0)",
    )
    overlap_parser.add_argument(
        "--damage-mismatch-penalty",
        type=float,
        default=0.25,
        help="fractional error penalty for terminal C/T and G/A mismatches (default: 0.25)",
    )
    overlap_parser.add_argument(
        "--ranking-damage-end-window",
        type=_end_window,
        default=5,
        help="terminal bases eligible for damage-aware ranking (default: 5)",
    )
    overlap_parser.add_argument(
        "--cross-cluster-recruitment-audit",
        action="store_true",
        help="audit cross-cluster fragments against frozen contig ends without changing output",
    )
    overlap_parser.add_argument(
        "--read-supported-contig-link-audit",
        action="store_true",
        help="audit reciprocal contig links crossed by independent reads",
    )
    overlap_parser.add_argument(
        "--min-contig-link-read-support",
        type=_minimum_count,
        default=2,
        help="minimum spanning molecules for a projected contig link (default: 2)",
    )
    overlap_parser.add_argument(
        "--two-tier-redundancy-audit",
        action="store_true",
        help="audit support-three rescue after 97%%/99%% containment filtering",
    )
    overlap_parser.add_argument(
        "--rescue-min-cluster-size",
        type=_minimum_count,
        default=3,
        help="minimum cluster size for two-tier rescue candidates (default: 3)",
    )
    overlap_parser.add_argument(
        "--iterative-reclustering-audit",
        action="store_true",
        help="project fixed-pool global reclustering without changing output",
    )
    overlap_parser.add_argument(
        "--max-recluster-iterations",
        type=_minimum_count,
        default=4,
        help="maximum report-only global reclustering rounds (default: 4)",
    )
    overlap_parser.add_argument(
        "--recluster-derived-min-identity",
        type=float,
        default=0.99,
        help="minimum identity for overlaps involving derived sequences (default: 0.99)",
    )
    overlap_parser.add_argument(
        "--recluster-derived-min-raw-support",
        type=_minimum_count,
        default=2,
        help="minimum pristine source molecules supporting derived extension (default: 2)",
    )
    overlap_parser.add_argument(
        "--recluster-mismatch-min-raw-support",
        type=_minimum_count,
        default=3,
        help="minimum Q20 raw molecules confirming a derived mismatch (default: 3)",
    )
    overlap_parser.add_argument(
        "--recluster-mismatch-min-raw-margin",
        type=_minimum_count,
        default=2,
        help="minimum molecule advantage over a competing mismatch allele (default: 2)",
    )
    overlap_parser.add_argument(
        "--recluster-mismatch-min-base-quality",
        type=_end_window,
        default=20,
        help="minimum Phred score for mismatch evidence (default: 20)",
    )
    overlap_parser.add_argument(
        "--iterative-reclustering-report",
        type=Path,
        help="optional per-round iterative reclustering TSV",
    )
    overlap_parser.add_argument(
        "--iterative-reclustering-projection",
        type=Path,
        help="optional projected FASTA; the primary output remains unchanged",
    )
    overlap_parser.add_argument(
        "--progressive-raw-phase-audit",
        action="store_true",
        help="audit lifecycle-aware raw clustering and extension without changing output",
    )
    overlap_parser.add_argument(
        "--progressive-raw-phase-projection",
        type=Path,
        help="optional first-phase progressive FASTA; primary output remains unchanged",
    )
    overlap_parser.add_argument(
        "--raw-consensus-profile-prefix", type=Path,
        help="opt-in fixed-layout consensus for merged double-stranded fragments; reads PREFIX5p.prof and PREFIX3p.prof",
    )
    overlap_parser.add_argument("--raw-consensus-projection", type=Path, help="damage-aware consensus FASTA from the progressive layout")
    overlap_parser.add_argument("--raw-consensus-quality-control", type=Path, help="same layout and quality model with damage rates set to zero")
    overlap_parser.add_argument("--raw-consensus-report", type=Path, help="TSV of consensus decisions at 0-based positions")
    overlap_parser.add_argument("--raw-consensus-placements", type=Path, help="TSV of original raw-read intervals and contig offsets")
    overlap_parser.add_argument("--raw-consensus-mixture-report", type=Path, help="diagnostic single versus two-allele evidence TSV; does not change calls")
    overlap_parser.add_argument("--raw-consensus-linkage-report", type=Path, help="diagnostic molecule-level neighbouring allele counts TSV")
    overlap_parser.add_argument("--raw-consensus-linked-allele-guard", action="store_true", help="experimental guard against substitutions contradicting linked allele groups")
    overlap_parser.add_argument("--progressive-low-quality-ry-rescue", action="store_true", help="experimental recruitment allowance for one RY mismatch with raw base quality <=15; requires progressive raw audit")
    overlap_parser.add_argument(
        "--max-progressive-raw-iterations",
        type=_minimum_count,
        default=3,
        help="maximum persistent-center raw-extension rounds (default: 3)",
    )
    overlap_parser.add_argument(
        "--evidence-priority-progressive-extension",
        action="store_true",
        help="give conflicting unused reads to the best-supported active contig",
    )
    overlap_parser.add_argument(
        "--progressive-extension-failure-audit",
        action="store_true",
        help="summarize rejected extension-side support without changing output",
    )
    overlap_parser.add_argument(
        "--support-three-progressive-extension-audit",
        action="store_true",
        help="audit unanimous Q20 support-three extensions in a shadow pool",
    )
    overlap_parser.add_argument(
        "--support-three-progressive-projection",
        type=Path,
        help="optional support-three shadow FASTA; trusted output remains unchanged",
    )
    overlap_parser.add_argument(
        "--adaptive-support-rescue-audit",
        action="store_true",
        help="audit support-three raw clusters that exactly attach to trusted contigs",
    )
    overlap_parser.add_argument(
        "--adaptive-support-rescue-projection",
        type=Path,
        help="optional adaptive rescue FASTA; trusted outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--adaptive-rescue-min-support",
        type=_minimum_count,
        default=3,
        help="minimum independent molecules for adaptive rescue (default: 3)",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-adaptive-rescue-links",
        action="store_true",
        help="allow Q20 raw-confirmed near-exact adaptive-rescue links",
    )
    overlap_parser.add_argument(
        "--adaptive-rescue-chains",
        action="store_true",
        help="allow rescue contigs in primary-anchored branchless paths",
    )
    overlap_parser.add_argument(
        "--selective-support-rescue-audit",
        action="store_true",
        help="project novel support-three clusters from unused raw molecules",
    )
    overlap_parser.add_argument(
        "--selective-support-rescue-projection",
        type=Path,
        help="optional selective rescue FASTA; trusted outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--high-confidence-support-two-rescue-audit",
        action="store_true",
        help="project Q20, nonconflicting two-molecule rescue clusters",
    )
    overlap_parser.add_argument(
        "--high-confidence-support-two-rescue-projection",
        type=Path,
        help="optional support-two rescue FASTA; prior outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--support-two-seed-quality-filter",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="experimental whole-cluster rejection for unsupported low-quality seed bases (default: disabled)",
    )
    overlap_parser.add_argument(
        "--support-two-trim-seed-tails",
        action="store_true",
        help="trim exposed unsupported low-quality seed tails after extension; incompatible with seed rejection",
    )
    overlap_parser.add_argument(
        "--support-two-fixed-membership-before",
        type=Path,
        help="write matched untrimmed FASTA and hold rescue membership fixed; requires seed trimming",
    )
    overlap_parser.add_argument(
        "--support-two-recruit-quality-filter",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="experimental quality gate on later support-two recruit overhangs (default: enabled)",
    )
    overlap_parser.add_argument(
        "--support-two-master-graph-audit",
        action="store_true",
        help="project exact branchless paths from support-two rescue contigs",
    )
    overlap_parser.add_argument(
        "--support-two-master-graph-projection",
        type=Path,
        help="optional exact support-two graph FASTA; prior outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--selective-rescue-master-graph-audit",
        action="store_true",
        help="project exact branchless paths from selective rescue contigs",
    )
    overlap_parser.add_argument(
        "--selective-rescue-master-graph-projection",
        type=Path,
        help="optional exact selective-rescue graph FASTA; prior outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-selective-master-graph-audit",
        action="store_true",
        help="add raw-confirmed near-exact edges to the selective master graph",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-selective-master-graph-projection",
        type=Path,
        help="optional raw-confirmed selective graph FASTA; prior outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--raw-supported-progressive-link-audit",
        action="store_true",
        help="audit exact progressive-contig links uniquely spanned by raw molecules",
    )
    overlap_parser.add_argument(
        "--raw-supported-progressive-link-projection",
        type=Path,
        help="optional raw-supported progressive-link FASTA; prior outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--progressive-link-min-read-support",
        type=_minimum_count,
        default=2,
        help="minimum unique raw molecules spanning a progressive link (default: 2)",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-progressive-links",
        action="store_true",
        help="add Q20 raw-confirmed near-exact links after exact progressive links",
    )
    overlap_parser.add_argument(
        "--paired-progressive-scaffold-audit",
        action="store_true",
        help="audit reciprocal paired-end links between progressive contigs",
    )
    overlap_parser.add_argument(
        "--paired-progressive-scaffold-projection",
        type=Path,
        help="optional N-gapped paired scaffold FASTA; prior outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--paired-scaffold-min-support",
        type=_minimum_count,
        default=3,
        help="minimum independent read pairs supporting a scaffold link (default: 3)",
    )
    overlap_parser.add_argument(
        "--paired-scaffold-dominance-ratio",
        type=_dominance_ratio,
        default=3.0,
        help="best-to-runner-up scaffold-link support ratio (default: 3.0)",
    )
    overlap_parser.add_argument(
        "--paired-fragment-mean",
        type=float,
        help="known mean outer fragment length; otherwise estimate from mappings",
    )
    overlap_parser.add_argument(
        "--paired-fragment-sd",
        type=float,
        help="known fragment-length standard deviation; required with mean",
    )
    overlap_parser.add_argument(
        "--paired-scaffold-min-calibration-pairs",
        type=_minimum_count,
        default=20,
        help="same-contig pairs required for automatic calibration (default: 20)",
    )
    overlap_parser.add_argument(
        "--master-overlap-graph-audit",
        action="store_true",
        help="project exact branchless paths from iterative reclustering contigs",
    )
    overlap_parser.add_argument(
        "--master-overlap-graph-projection",
        type=Path,
        help="optional exact string-graph FASTA; accepted outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-master-graph-audit",
        action="store_true",
        help="supplement exact master paths with raw-confirmed near-exact overlaps",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-master-graph-projection",
        type=Path,
        help="optional raw-confirmed master-graph FASTA; prior outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-master-min-identity",
        type=_probability,
        default=0.98,
        help="minimum near-exact overlap identity (default: 0.98)",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-master-max-mismatches",
        type=_minimum_count,
        default=1,
        help="maximum mismatches in a raw-confirmed overlap (default: 1)",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-master-primary-support",
        type=_minimum_count,
        default=3,
        help="minimum Q20 raw molecules supporting the chosen allele (default: 3)",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-master-alternate-support",
        type=_minimum_count,
        default=2,
        help="raw molecules sufficient to identify a strain conflict (default: 2)",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-master-support-margin",
        type=_minimum_count,
        default=2,
        help="minimum chosen-over-competing allele support margin (default: 2)",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-master-min-base-quality",
        type=_end_window,
        default=20,
        help="minimum Phred score for raw mismatch evidence (default: 20)",
    )
    overlap_parser.add_argument(
        "--raw-confirmed-master-damage-end-window",
        type=_end_window,
        default=5,
        help="ignore damage-like evidence this close to read ends (default: 5)",
    )
    overlap_parser.add_argument(
        "--strain-safe-containment-audit",
        action="store_true",
        help="remove contained copies only when raw alleles do not support them",
    )
    overlap_parser.add_argument(
        "--strain-safe-containment-projection",
        type=Path,
        help="optional strain-safe containment FASTA; prior outputs remain unchanged",
    )
    overlap_parser.add_argument(
        "--containment-min-identity",
        type=_probability,
        default=0.99,
        help="minimum contained-sequence identity (default: 0.99)",
    )
    overlap_parser.add_argument(
        "--containment-min-coverage",
        type=_probability,
        default=0.99,
        help="minimum fraction of a candidate covered (default: 0.99)",
    )
    overlap_parser.add_argument(
        "--containment-primary-allele-support",
        type=_minimum_count,
        default=3,
        help="minimum Q20 raw molecules supporting the retained allele (default: 3)",
    )
    overlap_parser.add_argument(
        "--containment-alternate-allele-support",
        type=_minimum_count,
        default=2,
        help="raw molecules sufficient to protect an alternate allele (default: 2)",
    )
    overlap_parser.add_argument(
        "--containment-support-margin",
        type=_minimum_count,
        default=2,
        help="minimum retained-over-alternate allele support margin (default: 2)",
    )
    overlap_parser.add_argument(
        "--containment-min-base-quality",
        type=_end_window,
        default=20,
        help="minimum Phred score for containment allele evidence (default: 20)",
    )
    overlap_parser.add_argument(
        "--containment-damage-end-window",
        type=_end_window,
        default=5,
        help="ignore damage-like alleles this close to raw-read ends (default: 5)",
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


def main(argv: Sequence[str] | None = None) -> int:
    """Run the Anvaya command-line interface."""
    parser = build_parser()
    arguments = parser.parse_args(argv)

    try:
        if arguments.command == "carpedeam-assemble":
            _progress("Running audited CarpeDeam safe-mode comparator")
            result = run_carpedeam_safe(
                arguments.input,
                arguments.profile_prefix,
                arguments.output,
                arguments.temporary_directory,
                arguments.diagnostics,
                executable=arguments.executable,
                threads=arguments.threads,
                minimum_contig_length=arguments.min_contig_length,
            )
            _progress(f"Completed CarpeDeam safe mode in {result['elapsed_seconds']:.2f}s")
            print(f"reads={result['input_scan']['reads']}")
            print(f"contigs={result['output_contigs']}")
            print(f"output={arguments.output}")
            print(f"diagnostics={arguments.diagnostics}")
            return 0

        if arguments.command == "damage-assemble":
            if (arguments.output_assembled is None) != (arguments.output_unresolved is None):
                parser.error("assembled and unresolved outputs must be requested together")
            paths = [
                arguments.input,
                arguments.output_before,
                arguments.output_consensus,
                arguments.output_assembled,
                arguments.output_unresolved,
                arguments.diagnostics,
                arguments.placements_report,
                arguments.consensus_report,
                Path(f"{arguments.profile_prefix}5p.prof"),
                Path(f"{arguments.profile_prefix}3p.prof"),
            ]
            resolved = [path.resolve() for path in paths if path is not None]
            if len(resolved) != len(set(resolved)):
                parser.error("damage assembly inputs, profiles and outputs must be distinct")
            started = time.perf_counter()
            _progress(f"Loading fragments from {arguments.input}")
            reads = load_reads(arguments.input, maximum_reads=arguments.max_reads)
            if any(read.qualities is None for read in reads):
                parser.error("damage assembly requires FASTQ base qualities")
            profile = DamageProfile.from_prefix(arguments.profile_prefix)
            _progress(f"Building damage-aware graph from {len(reads)} fragments")
            graph_timings = {}
            pool, graph = assemble_damage_string_graph(
                reads,
                damage_profile=profile,
                maximum_reads=arguments.max_reads,
                stage_timings=graph_timings,
            )
            layout_contigs = len(pool.active_derived)
            retained_pool = ProgressiveSequencePool(tuple(
                record.consumed()
                if (
                    record.state in {
                        SequenceState.CORRECTED_CENTER,
                        SequenceState.EXTENDED_CONTIG,
                    }
                    and len(record.current.sequence) < arguments.min_output_length
                )
                else record
                for record in pool.records
            ))
            with ExitStack() as stack:
                report = (
                    stack.enter_context(arguments.consensus_report.open("w", encoding="utf-8"))
                    if arguments.consensus_report else None
                )
                placements = (
                    stack.enter_context(arguments.placements_report.open("w", encoding="utf-8"))
                    if arguments.placements_report else None
                )
                consensus_started = time.perf_counter()
                polished, consensus = project_raw_consensus(
                    retained_pool, profile, report=report, placements_report=placements
                )
                consensus_seconds = time.perf_counter() - consensus_started
            before = [
                record.current.sequence for record in retained_pool.active_derived
            ]
            after = [read.sequence for read in polished]
            write_fasta(before, arguments.output_before)
            write_fasta(after, arguments.output_consensus)
            assembled = []
            unresolved = []
            if arguments.output_assembled is not None:
                for index, (record, read) in enumerate(zip(
                    retained_pool.active_derived, polished, strict=True,
                ), start=1):
                    mappings = {}
                    for placement in record.raw_placements or ():
                        mappings.setdefault(placement.read_index, set()).add(
                            (placement.offset, placement.reverse)
                        )
                    ambiguous_molecules = {
                        retained_pool.records[read_index].molecule_id
                        for read_index, positions in mappings.items()
                        if len(positions) > 1
                    }
                    eligible_molecules = {
                        retained_pool.records[placement.read_index].molecule_id
                        for placement in record.raw_placements or ()
                        if retained_pool.records[placement.read_index].molecule_id
                        not in ambiguous_molecules
                    }
                    destination = assembled if len(eligible_molecules) >= 2 else unresolved
                    destination.append((f"unitig_{index}", read.sequence))
                write_named_fasta(assembled, arguments.output_assembled)
                write_named_fasta(unresolved, arguments.output_unresolved)
            elapsed = time.perf_counter() - started
            arguments.diagnostics.write_text(
                json.dumps(
                    {
                        "input_reads": len(reads),
                        "input_bases": sum(len(read.sequence) for read in reads),
                        "layout_contigs": layout_contigs,
                        "filtered_short_contigs": layout_contigs - len(before),
                        "emitted_contigs": len(before),
                        "emitted_bases_before": sum(map(len, before)),
                        "emitted_bases_consensus": sum(map(len, after)),
                        "assembled_minimum_distinct_molecules": 2,
                        "assembled_contigs": len(assembled),
                        "assembled_bases": sum(len(sequence) for _, sequence in assembled),
                        "unresolved_contigs": len(unresolved),
                        "unresolved_bases": sum(len(sequence) for _, sequence in unresolved),
                        "minimum_output_length": arguments.min_output_length,
                        "graph": graph,
                        "consensus": asdict(consensus),
                        "stage_seconds": {
                            **graph_timings,
                            "consensus": consensus_seconds,
                        },
                        "elapsed_seconds": elapsed,
                    },
                    indent=2,
                ) + "\n",
                encoding="utf-8",
            )
            _progress(f"Completed damage-aware graph in {elapsed:.2f}s")
            print(f"reads={len(reads)}")
            print(f"contigs={len(before)}")
            print(f"output_before={arguments.output_before}")
            print(f"output_consensus={arguments.output_consensus}")
            if arguments.output_assembled is not None:
                print(f"output_assembled={arguments.output_assembled}")
                print(f"output_unresolved={arguments.output_unresolved}")
            print(f"diagnostics={arguments.diagnostics}")
            return 0

        if arguments.command == "overlap-assemble":
            if arguments.progressive_low_quality_ry_rescue and not arguments.progressive_raw_phase_audit:
                parser.error("low-quality RY rescue requires --progressive-raw-phase-audit")
            consensus_profile = None
            consensus_outputs = [arguments.raw_consensus_projection, arguments.raw_consensus_quality_control,
                                 arguments.raw_consensus_report, arguments.raw_consensus_placements,
                                 arguments.raw_consensus_mixture_report, arguments.raw_consensus_linkage_report]
            if arguments.raw_consensus_profile_prefix is not None:
                if not arguments.progressive_raw_phase_audit or arguments.raw_consensus_projection is None:
                    parser.error("raw consensus requires progressive raw audit and a consensus projection")
                if arguments.input is None:
                    parser.error("raw consensus currently requires a single input of merged double-stranded fragments")
                consensus_profile = DamageProfile.from_prefix(arguments.raw_consensus_profile_prefix)
                paths = [path.resolve() for path in consensus_outputs if path is not None]
                reserved = {path.resolve() for path in (arguments.input, arguments.output,
                    arguments.progressive_raw_phase_projection, arguments.selective_support_rescue_projection,
                    arguments.high_confidence_support_two_rescue_projection,
                    arguments.support_two_fixed_membership_before) if path is not None}
                reserved.update(Path(f"{arguments.raw_consensus_profile_prefix}{end}.prof").resolve() for end in ("5p", "3p"))
                if len(set(paths)) != len(paths) or set(paths) & reserved:
                    parser.error("raw consensus outputs must be distinct from each other and existing inputs/projections")
            elif any(path is not None for path in consensus_outputs) or arguments.raw_consensus_linked_allele_guard:
                parser.error("raw consensus outputs require --raw-consensus-profile-prefix")
            if arguments.support_two_fixed_membership_before is not None:
                if not arguments.support_two_trim_seed_tails:
                    parser.error("fixed-membership comparison requires --support-two-trim-seed-tails")
                after_path = arguments.high_confidence_support_two_rescue_projection
                if after_path is None:
                    parser.error("fixed-membership comparison requires a support-two projection path")
                if arguments.support_two_fixed_membership_before.resolve() == after_path.resolve():
                    parser.error("fixed-membership before and after paths must differ")
            if arguments.support_two_trim_seed_tails:
                if arguments.support_two_seed_quality_filter:
                    parser.error("support-two seed rejection and trimming cannot be combined")
                if not arguments.high_confidence_support_two_rescue_audit:
                    parser.error("seed trimming requires --high-confidence-support-two-rescue-audit")
            started = time.perf_counter()
            input_paths = _resolve_input_paths(parser, arguments)
            _progress(f"Loading fragments from {', '.join(map(str, input_paths))}")
            read_groups = [load_reads(path) for path in input_paths]
            if len(read_groups) == 2 and len(read_groups[0]) != len(read_groups[1]):
                parser.error(
                    "paired-end input files must contain the same number of reads"
                )
            paired_merge = PairedMergeDiagnostics()
            if arguments.merge_overlapping_pairs:
                if len(read_groups) != 2:
                    parser.error("--merge-overlapping-pairs requires paired-end input")
                reads, molecule_ids, paired_merge = merge_overlapping_pairs(
                    read_groups[0],
                    read_groups[1],
                    minimum_overlap=arguments.paired_merge_min_overlap,
                )
            else:
                reads = [read for group in read_groups for read in group]
                molecule_ids = (
                    list(range(len(read_groups[0]))) * 2
                    if len(read_groups) == 2
                    else list(range(len(reads)))
                )
            _progress(f"Loaded {len(reads)} fragments")
            correction_events: list[OverlapCorrectionEvent] = []
            contigs, summary = assemble_overlap_contigs(
                reads,
                molecule_ids=molecule_ids,
                anchor_k=arguments.anchor_k,
                anchors_per_read=arguments.anchors_per_read,
                maximum_anchor_occurrences=arguments.max_anchor_occurrences,
                minimum_anchor_matches=arguments.min_anchor_matches,
                minimum_overlap=arguments.min_overlap,
                ranked_extension=arguments.ranked_extension,
                extension_consensus=arguments.extension_consensus,
                minimum_ranked_extension_support=arguments.min_ranked_extension_support,
                reciprocal_best_extension=arguments.reciprocal_best_extension,
                damage_aware_ranking=arguments.damage_aware_ranking,
                minimum_confidence_margin=arguments.min_overlap_confidence_margin,
                damage_mismatch_penalty=arguments.damage_mismatch_penalty,
                ranking_damage_end_window=arguments.ranking_damage_end_window,
                cross_cluster_recruitment_audit=(
                    arguments.cross_cluster_recruitment_audit
                ),
                read_supported_contig_link_audit=(
                    arguments.read_supported_contig_link_audit
                ),
                minimum_contig_link_read_support=(
                    arguments.min_contig_link_read_support
                ),
                maximum_rounds=arguments.max_rounds,
                maximum_contig_iterations=arguments.max_contig_iterations,
                minimum_cluster_size=arguments.min_cluster_size,
                minimum_output_length=arguments.min_output_length,
                damage_end_window=arguments.damage_end_window,
                correction_events=correction_events,
            )
            rescue_audit = TwoTierRedundancyDiagnostics()
            if arguments.two_tier_redundancy_audit:
                if arguments.rescue_min_cluster_size >= arguments.min_cluster_size:
                    parser.error(
                        "--rescue-min-cluster-size must be smaller than "
                        "--min-cluster-size"
                    )
                rescue_contigs, _ = assemble_overlap_contigs(
                    reads,
                    molecule_ids=molecule_ids,
                    anchor_k=arguments.anchor_k,
                    anchors_per_read=arguments.anchors_per_read,
                    maximum_anchor_occurrences=arguments.max_anchor_occurrences,
                    minimum_anchor_matches=arguments.min_anchor_matches,
                    minimum_overlap=arguments.min_overlap,
                    ranked_extension=arguments.ranked_extension,
                    extension_consensus=arguments.extension_consensus,
                    minimum_ranked_extension_support=(
                        arguments.min_ranked_extension_support
                    ),
                    reciprocal_best_extension=(
                        arguments.reciprocal_best_extension
                    ),
                    damage_aware_ranking=arguments.damage_aware_ranking,
                    minimum_confidence_margin=(
                        arguments.min_overlap_confidence_margin
                    ),
                    damage_mismatch_penalty=arguments.damage_mismatch_penalty,
                    ranking_damage_end_window=(
                        arguments.ranking_damage_end_window
                    ),
                    maximum_rounds=arguments.max_rounds,
                    maximum_contig_iterations=(
                        arguments.max_contig_iterations
                    ),
                    minimum_cluster_size=arguments.rescue_min_cluster_size,
                    minimum_output_length=arguments.min_output_length,
                    damage_end_window=arguments.damage_end_window,
                )
                rescue_audit = audit_two_tier_redundancy(
                    contigs,
                    rescue_contigs,
                    anchor_k=arguments.anchor_k,
                    anchors_per_read=arguments.anchors_per_read,
                    maximum_anchor_occurrences=(
                        arguments.max_anchor_occurrences
                    ),
                    minimum_anchor_matches=arguments.min_anchor_matches,
                )
            reclustering_audit = IterativeReclusteringDiagnostics()
            progressive_clustering = ProgressiveClusteringDiagnostics()
            progressive_extension = ProgressiveRawExtensionDiagnostics()
            progressive_iterations = ProgressiveIterationDiagnostics()
            support_three_iterations = ProgressiveIterationDiagnostics()
            adaptive_rescue = AdaptiveRescueDiagnostics()
            selective_rescue = SelectiveRescueDiagnostics()
            support_two_rescue = SupportTwoRescueDiagnostics()
            support_two_master_graph = MasterOverlapGraphDiagnostics()
            selective_rescue_master_graph = MasterOverlapGraphDiagnostics()
            raw_selective_master_graph = RawConfirmedMasterGraphDiagnostics()
            progressive_links = ProgressiveLinkDiagnostics()
            paired_scaffold = PairedScaffoldDiagnostics()
            master_graph_audit = MasterOverlapGraphDiagnostics()
            raw_master_graph_audit = RawConfirmedMasterGraphDiagnostics()
            containment_audit = StrainSafeContainmentDiagnostics()
            if arguments.progressive_raw_phase_audit:
                if not (
                    arguments.ranked_extension
                    and arguments.reciprocal_best_extension
                    and arguments.damage_aware_ranking
                ):
                    parser.error(
                        "--progressive-raw-phase-audit requires ranked, "
                        "reciprocal-best, and damage-aware extension"
                    )
                progressive_pool = ProgressiveSequencePool.from_reads(
                    reads, molecule_ids
                )
                progressive_clusters, progressive_clustering = (
                    discover_progressive_raw_clusters(
                        progressive_pool,
                        anchor_k=arguments.anchor_k,
                        anchors_per_read=arguments.anchors_per_read,
                        maximum_anchor_occurrences=(
                            arguments.max_anchor_occurrences
                        ),
                        minimum_anchor_matches=arguments.min_anchor_matches,
                        minimum_overlap=arguments.min_overlap,
                        minimum_cluster_size=arguments.min_cluster_size,
                    )
                )
                progressive_pool, progressive_extension = (
                    extend_progressive_raw_clusters(
                        progressive_pool,
                        progressive_clusters,
                        anchor_k=arguments.anchor_k,
                        anchors_per_read=arguments.anchors_per_read,
                        maximum_anchor_occurrences=(
                            arguments.max_anchor_occurrences
                        ),
                        minimum_anchor_matches=arguments.min_anchor_matches,
                        minimum_overlap=arguments.min_overlap,
                        minimum_extension_support=(
                            arguments.min_ranked_extension_support
                        ),
                        minimum_confidence_margin=(
                            arguments.min_overlap_confidence_margin
                        ),
                        damage_mismatch_penalty=arguments.damage_mismatch_penalty,
                        damage_end_window=arguments.ranking_damage_end_window,
                        reciprocal_best_extension=True,
                        extension_consensus=arguments.extension_consensus,
                        maximum_rounds=arguments.max_rounds,
                        track_raw_placements=consensus_profile is not None,
                        low_quality_ry_rescue=arguments.progressive_low_quality_ry_rescue,
                    )
                )
                progressive_pool, progressive_iterations = (
                    iterate_progressive_raw_extension(
                        progressive_pool,
                        anchor_k=arguments.anchor_k,
                        anchors_per_read=arguments.anchors_per_read,
                        maximum_anchor_occurrences=(
                            arguments.max_anchor_occurrences
                        ),
                        minimum_anchor_matches=arguments.min_anchor_matches,
                        minimum_overlap=arguments.min_overlap,
                        maximum_iterations=(
                            arguments.max_progressive_raw_iterations
                        ),
                        minimum_confidence_margin=(
                            arguments.min_overlap_confidence_margin
                        ),
                        damage_mismatch_penalty=arguments.damage_mismatch_penalty,
                        damage_end_window=arguments.ranking_damage_end_window,
                        evidence_priority=(
                            arguments.evidence_priority_progressive_extension
                        ),
                        low_quality_ry_rescue=arguments.progressive_low_quality_ry_rescue,
                        audit_rejected_extensions=(
                            arguments.progressive_extension_failure_audit
                        ),
                    )
                )
                progressive_projection_contigs = [
                    record.current for record in progressive_pool.active_derived
                ]
                if consensus_profile is not None:
                    # Match the written progressive FASTA's length filter and IDs.
                    consensus_pool = ProgressiveSequencePool(tuple(
                        record.consumed() if record.state in (SequenceState.CORRECTED_CENTER, SequenceState.EXTENDED_CONTIG)
                        and len(record.current.sequence) < arguments.min_output_length else record
                        for record in progressive_pool.records
                    ))
                    with ExitStack() as stack:
                        report = stack.enter_context(arguments.raw_consensus_report.open("w", encoding="utf-8", newline="")) if arguments.raw_consensus_report else None
                        placements = stack.enter_context(arguments.raw_consensus_placements.open("w", encoding="utf-8", newline="")) if arguments.raw_consensus_placements else None
                        mixture = stack.enter_context(arguments.raw_consensus_mixture_report.open("w", encoding="utf-8", newline="")) if arguments.raw_consensus_mixture_report else None
                        linkage = stack.enter_context(arguments.raw_consensus_linkage_report.open("w", encoding="utf-8", newline="")) if arguments.raw_consensus_linkage_report else None
                        consensus_contigs, consensus_diagnostics = project_raw_consensus(
                            consensus_pool, consensus_profile, report=report, placements_report=placements,
                            mixture_report=mixture, linkage_report=linkage,
                            linked_allele_guard=arguments.raw_consensus_linked_allele_guard,
                        )
                    write_fasta([read.sequence for read in consensus_contigs], arguments.raw_consensus_projection)
                    print(f"raw_consensus_profile_prefix={arguments.raw_consensus_profile_prefix}")
                    print(f"raw_consensus_linked_allele_guard={str(arguments.raw_consensus_linked_allele_guard).lower()}")
                    for label, value in asdict(consensus_diagnostics).items():
                        print(f"raw_consensus_{label}={value}")
                    if arguments.raw_consensus_quality_control is not None:
                        quality_contigs, quality_diagnostics = project_raw_consensus(consensus_pool, DamageProfile((), ()))
                        write_fasta([read.sequence for read in quality_contigs], arguments.raw_consensus_quality_control)
                        for label, value in asdict(quality_diagnostics).items():
                            print(f"raw_quality_control_{label}={value}")
                if (
                    arguments.support_three_progressive_projection is not None
                    and not arguments.support_three_progressive_extension_audit
                ):
                    parser.error(
                        "support-three projection requires "
                        "--support-three-progressive-extension-audit"
                    )
                if arguments.support_three_progressive_extension_audit:
                    support_three_pool, support_three_iterations = (
                        iterate_progressive_raw_extension(
                            progressive_pool,
                            anchor_k=arguments.anchor_k,
                            anchors_per_read=arguments.anchors_per_read,
                            maximum_anchor_occurrences=(
                                arguments.max_anchor_occurrences
                            ),
                            minimum_anchor_matches=(
                                arguments.min_anchor_matches
                            ),
                            minimum_overlap=arguments.min_overlap,
                            minimum_consensus_support=3,
                            maximum_iterations=(
                                arguments.max_progressive_raw_iterations
                            ),
                            minimum_confidence_margin=(
                                arguments.min_overlap_confidence_margin
                            ),
                            damage_mismatch_penalty=(
                                arguments.damage_mismatch_penalty
                            ),
                            damage_end_window=(
                                arguments.ranking_damage_end_window
                            ),
                            audit_minimum_base_quality=20,
                            require_strict_boundary_evidence=True,
                        )
                    )
                    if arguments.support_three_progressive_projection is not None:
                        write_fasta(
                            [
                                record.current.sequence
                                for record in support_three_pool.active_derived
                                if len(record.current.sequence)
                                >= arguments.min_output_length
                            ],
                            arguments.support_three_progressive_projection,
                        )
                if (
                    arguments.adaptive_support_rescue_projection is not None
                    and not arguments.adaptive_support_rescue_audit
                ):
                    parser.error(
                        "adaptive rescue projection requires "
                        "--adaptive-support-rescue-audit"
                    )
                if (
                    arguments.raw_confirmed_adaptive_rescue_links
                    and not arguments.adaptive_support_rescue_audit
                ):
                    parser.error(
                        "--raw-confirmed-adaptive-rescue-links requires "
                        "--adaptive-support-rescue-audit"
                    )
                if (
                    arguments.adaptive_rescue_chains
                    and not arguments.adaptive_support_rescue_audit
                ):
                    parser.error(
                        "--adaptive-rescue-chains requires "
                        "--adaptive-support-rescue-audit"
                    )
                if arguments.adaptive_support_rescue_audit:
                    adaptive_contigs, adaptive_rescue = (
                        audit_adaptive_support_rescue(
                            progressive_pool,
                            minimum_rescue_support=(
                                arguments.adaptive_rescue_min_support
                            ),
                            anchor_k=arguments.anchor_k,
                            anchors_per_read=arguments.anchors_per_read,
                            maximum_anchor_occurrences=(
                                arguments.max_anchor_occurrences
                            ),
                            minimum_anchor_matches=(
                                arguments.min_anchor_matches
                            ),
                            minimum_overlap=arguments.min_overlap,
                            allow_rescue_chains=arguments.adaptive_rescue_chains,
                            allow_near_exact=(
                                arguments.raw_confirmed_adaptive_rescue_links
                            ),
                            near_exact_minimum_identity=(
                                arguments.raw_confirmed_master_min_identity
                            ),
                            near_exact_maximum_mismatches=(
                                arguments.raw_confirmed_master_max_mismatches
                            ),
                            minimum_primary_allele_support=(
                                arguments.raw_confirmed_master_primary_support
                            ),
                            minimum_alternate_allele_support=(
                                arguments.raw_confirmed_master_alternate_support
                            ),
                            minimum_support_margin=(
                                arguments.raw_confirmed_master_support_margin
                            ),
                            minimum_base_quality=(
                                arguments.raw_confirmed_master_min_base_quality
                            ),
                            damage_end_window=(
                                arguments.raw_confirmed_master_damage_end_window
                            ),
                        )
                    )
                    if arguments.adaptive_support_rescue_projection is not None:
                        write_fasta(
                            [
                                contig.sequence
                                for contig in adaptive_contigs
                                if len(contig.sequence)
                                >= arguments.min_output_length
                            ],
                            arguments.adaptive_support_rescue_projection,
                        )
                if (
                    arguments.selective_support_rescue_projection is not None
                    and not arguments.selective_support_rescue_audit
                ):
                    parser.error(
                        "selective rescue projection requires "
                        "--selective-support-rescue-audit"
                    )
                if (
                    arguments.selective_rescue_master_graph_audit
                    and not arguments.selective_support_rescue_audit
                ):
                    parser.error(
                        "--selective-rescue-master-graph-audit requires "
                        "--selective-support-rescue-audit"
                    )
                if (
                    arguments.high_confidence_support_two_rescue_audit
                    and not arguments.selective_support_rescue_audit
                ):
                    parser.error(
                        "--high-confidence-support-two-rescue-audit requires "
                        "--selective-support-rescue-audit"
                    )
                if (
                    arguments.high_confidence_support_two_rescue_projection
                    is not None
                    and not arguments.high_confidence_support_two_rescue_audit
                ):
                    parser.error(
                        "support-two rescue projection requires "
                        "--high-confidence-support-two-rescue-audit"
                    )
                if (
                    arguments.support_two_master_graph_audit
                    and not arguments.high_confidence_support_two_rescue_audit
                ):
                    parser.error(
                        "--support-two-master-graph-audit requires "
                        "--high-confidence-support-two-rescue-audit"
                    )
                if (
                    arguments.support_two_master_graph_projection is not None
                    and not arguments.support_two_master_graph_audit
                ):
                    parser.error(
                        "support-two master-graph projection requires "
                        "--support-two-master-graph-audit"
                    )
                if (
                    arguments.selective_rescue_master_graph_projection is not None
                    and not arguments.selective_rescue_master_graph_audit
                ):
                    parser.error(
                        "selective rescue master-graph projection requires "
                        "--selective-rescue-master-graph-audit"
                    )
                if (
                    arguments.raw_confirmed_selective_master_graph_audit
                    and not arguments.selective_rescue_master_graph_audit
                ):
                    parser.error(
                        "--raw-confirmed-selective-master-graph-audit requires "
                        "--selective-rescue-master-graph-audit"
                    )
                if (
                    arguments.raw_confirmed_selective_master_graph_projection
                    is not None
                    and not arguments.raw_confirmed_selective_master_graph_audit
                ):
                    parser.error(
                        "raw-confirmed selective master-graph projection requires "
                        "--raw-confirmed-selective-master-graph-audit"
                    )
                if arguments.selective_support_rescue_audit:
                    selective_contigs, selective_rescue = (
                        project_selective_support_rescue(
                            progressive_pool,
                            minimum_rescue_support=(
                                arguments.adaptive_rescue_min_support
                            ),
                            anchor_k=arguments.anchor_k,
                            anchors_per_read=arguments.anchors_per_read,
                            maximum_anchor_occurrences=(
                                arguments.max_anchor_occurrences
                            ),
                            minimum_anchor_matches=(
                                arguments.min_anchor_matches
                            ),
                            minimum_overlap=arguments.min_overlap,
                        )
                    )
                    if (
                        arguments.selective_support_rescue_projection
                        is not None
                    ):
                        write_fasta(
                            [
                                contig.sequence
                                for contig in selective_contigs
                                if len(contig.sequence)
                                >= arguments.min_output_length
                            ],
                            arguments.selective_support_rescue_projection,
                        )
                    if arguments.high_confidence_support_two_rescue_audit:
                        fixed_before = [] if arguments.support_two_fixed_membership_before is not None else None
                        support_two_contigs, support_two_rescue = (
                            project_high_confidence_support_two_rescue(
                                progressive_pool,
                                selective_contigs,
                                filter_seed_quality=arguments.support_two_seed_quality_filter,
                                filter_recruit_quality=arguments.support_two_recruit_quality_filter,
                                trim_seed_tails=arguments.support_two_trim_seed_tails,
                                minimum_output_length=arguments.min_output_length,
                                fixed_membership_before=fixed_before,
                                anchor_k=arguments.anchor_k,
                                anchors_per_read=arguments.anchors_per_read,
                                maximum_anchor_occurrences=min(
                                    20, arguments.max_anchor_occurrences
                                ),
                                minimum_anchor_matches=max(
                                    2, arguments.min_anchor_matches
                                ),
                                minimum_overlap=arguments.min_overlap,
                                minimum_base_quality=(
                                    arguments.raw_confirmed_master_min_base_quality
                                ),
                                damage_end_window=(
                                    arguments.raw_confirmed_master_damage_end_window
                                ),
                            )
                        )
                        if (
                            arguments.high_confidence_support_two_rescue_projection
                            is not None
                        ):
                            write_fasta(
                                [
                                    contig.sequence
                                    for contig in support_two_contigs
                                    if len(contig.sequence)
                                    >= arguments.min_output_length
                                ],
                                arguments.high_confidence_support_two_rescue_projection,
                            )
                        if fixed_before is not None:
                            write_fasta(
                                [contig.sequence for contig in fixed_before if len(contig.sequence) >= arguments.min_output_length],
                                arguments.support_two_fixed_membership_before,
                            )
                        if arguments.support_two_master_graph_audit:
                            (
                                support_two_master_contigs,
                                support_two_master_graph,
                            ) = audit_master_overlap_graph(
                                support_two_contigs,
                                anchor_k=arguments.anchor_k,
                                anchors_per_read=arguments.anchors_per_read,
                                maximum_anchor_occurrences=(
                                    arguments.max_anchor_occurrences
                                ),
                                minimum_anchor_matches=(
                                    arguments.min_anchor_matches
                                ),
                                minimum_overlap=arguments.min_overlap,
                            )
                            if (
                                arguments.support_two_master_graph_projection
                                is not None
                            ):
                                write_fasta(
                                    [
                                        contig.sequence
                                        for contig in support_two_master_contigs
                                        if len(contig.sequence)
                                        >= arguments.min_output_length
                                    ],
                                    arguments.support_two_master_graph_projection,
                                )
                    if arguments.selective_rescue_master_graph_audit:
                        (
                            selective_master_contigs,
                            selective_rescue_master_graph,
                        ) = audit_master_overlap_graph(
                            selective_contigs,
                            anchor_k=arguments.anchor_k,
                            anchors_per_read=arguments.anchors_per_read,
                            maximum_anchor_occurrences=(
                                arguments.max_anchor_occurrences
                            ),
                            minimum_anchor_matches=(
                                arguments.min_anchor_matches
                            ),
                            minimum_overlap=arguments.min_overlap,
                        )
                        if (
                            arguments.selective_rescue_master_graph_projection
                            is not None
                        ):
                            write_fasta(
                                [
                                    contig.sequence
                                    for contig in selective_master_contigs
                                    if len(contig.sequence)
                                    >= arguments.min_output_length
                                ],
                                arguments.selective_rescue_master_graph_projection,
                            )
                        if arguments.raw_confirmed_selective_master_graph_audit:
                            (
                                raw_selective_master_contigs,
                                raw_selective_master_graph,
                            ) = audit_raw_confirmed_master_overlap_graph(
                                selective_master_contigs,
                                reads,
                                molecule_ids=molecule_ids,
                                anchor_k=arguments.anchor_k,
                                anchors_per_read=arguments.anchors_per_read,
                                maximum_anchor_occurrences=(
                                    arguments.max_anchor_occurrences
                                ),
                                minimum_anchor_matches=(
                                    arguments.min_anchor_matches
                                ),
                                minimum_overlap=arguments.min_overlap,
                                minimum_identity=(
                                    arguments.raw_confirmed_master_min_identity
                                ),
                                maximum_mismatches=(
                                    arguments.raw_confirmed_master_max_mismatches
                                ),
                                minimum_primary_allele_support=(
                                    arguments.raw_confirmed_master_primary_support
                                ),
                                minimum_alternate_allele_support=(
                                    arguments.raw_confirmed_master_alternate_support
                                ),
                                minimum_support_margin=(
                                    arguments.raw_confirmed_master_support_margin
                                ),
                                minimum_base_quality=(
                                    arguments.raw_confirmed_master_min_base_quality
                                ),
                                damage_end_window=(
                                    arguments.raw_confirmed_master_damage_end_window
                                ),
                            )
                            if (
                                arguments.raw_confirmed_selective_master_graph_projection
                                is not None
                            ):
                                write_fasta(
                                    [
                                        contig.sequence
                                        for contig in raw_selective_master_contigs
                                        if len(contig.sequence)
                                        >= arguments.min_output_length
                                    ],
                                    arguments.raw_confirmed_selective_master_graph_projection,
                                )
                if (
                    arguments.raw_supported_progressive_link_projection
                    is not None
                    and not arguments.raw_supported_progressive_link_audit
                ):
                    parser.error(
                        "raw-supported progressive-link projection requires "
                        "--raw-supported-progressive-link-audit"
                    )
                if (
                    arguments.raw_confirmed_progressive_links
                    and not arguments.raw_supported_progressive_link_audit
                ):
                    parser.error(
                        "--raw-confirmed-progressive-links requires "
                        "--raw-supported-progressive-link-audit"
                    )
                if (
                    arguments.paired_progressive_scaffold_projection is not None
                    and not arguments.paired_progressive_scaffold_audit
                ):
                    parser.error(
                        "paired progressive scaffold projection requires "
                        "--paired-progressive-scaffold-audit"
                    )
                if (
                    arguments.paired_progressive_scaffold_audit
                    and len(read_groups) != 2
                ):
                    parser.error(
                        "--paired-progressive-scaffold-audit requires paired-end "
                        "--left/-1 and --right/-2 input"
                    )
                scaffold_input = progressive_projection_contigs
                if arguments.raw_supported_progressive_link_audit:
                    progressive_link_contigs, progressive_links = (
                        audit_raw_supported_progressive_links(
                            progressive_pool,
                            anchor_k=arguments.anchor_k,
                            anchors_per_read=arguments.anchors_per_read,
                            maximum_anchor_occurrences=(
                                arguments.max_anchor_occurrences
                            ),
                            minimum_anchor_matches=(
                                arguments.min_anchor_matches
                            ),
                            minimum_overlap=arguments.min_overlap,
                            minimum_read_support=(
                                arguments.progressive_link_min_read_support
                            ),
                            allow_near_exact=(
                                arguments.raw_confirmed_progressive_links
                            ),
                            near_exact_minimum_identity=(
                                arguments.raw_confirmed_master_min_identity
                            ),
                            near_exact_maximum_mismatches=(
                                arguments.raw_confirmed_master_max_mismatches
                            ),
                            minimum_primary_allele_support=(
                                arguments.raw_confirmed_master_primary_support
                            ),
                            minimum_alternate_allele_support=(
                                arguments.raw_confirmed_master_alternate_support
                            ),
                            minimum_support_margin=(
                                arguments.raw_confirmed_master_support_margin
                            ),
                            minimum_base_quality=(
                                arguments.raw_confirmed_master_min_base_quality
                            ),
                            damage_end_window=(
                                arguments.raw_confirmed_master_damage_end_window
                            ),
                        )
                    )
                    if (
                        arguments.raw_supported_progressive_link_projection
                        is not None
                    ):
                        write_fasta(
                            [
                                contig.sequence
                                for contig in progressive_link_contigs
                                if len(contig.sequence)
                                >= arguments.min_output_length
                            ],
                            arguments.raw_supported_progressive_link_projection,
                        )
                    scaffold_input = progressive_link_contigs
                if arguments.paired_progressive_scaffold_audit:
                    scaffold_contigs, paired_scaffold = (
                        scaffold_progressive_contigs(
                            scaffold_input,
                            read_groups[0],
                            read_groups[1],
                            anchor_k=arguments.anchor_k,
                            anchors_per_read=arguments.anchors_per_read,
                            maximum_anchor_occurrences=(
                                arguments.max_anchor_occurrences
                            ),
                            minimum_anchor_matches=(
                                arguments.min_anchor_matches
                            ),
                            minimum_link_support=(
                                arguments.paired_scaffold_min_support
                            ),
                            dominance_ratio=(
                                arguments.paired_scaffold_dominance_ratio
                            ),
                            fragment_mean=arguments.paired_fragment_mean,
                            fragment_sd=arguments.paired_fragment_sd,
                            minimum_calibration_pairs=(
                                arguments.paired_scaffold_min_calibration_pairs
                            ),
                        )
                    )
                    if (
                        arguments.paired_progressive_scaffold_projection
                        is not None
                    ):
                        write_fasta(
                            [
                                contig.sequence
                                for contig in scaffold_contigs
                                if len(contig.sequence)
                                >= arguments.min_output_length
                            ],
                            arguments.paired_progressive_scaffold_projection,
                        )
                if arguments.progressive_raw_phase_projection is not None:
                    write_fasta(
                        [
                            contig.sequence
                            for contig in progressive_projection_contigs
                            if len(contig.sequence)
                            >= arguments.min_output_length
                        ],
                        arguments.progressive_raw_phase_projection,
                    )
            elif (
                arguments.progressive_raw_phase_projection is not None
                or arguments.evidence_priority_progressive_extension
                or arguments.progressive_extension_failure_audit
                or arguments.support_three_progressive_extension_audit
                or arguments.support_three_progressive_projection is not None
                or arguments.raw_supported_progressive_link_audit
                or arguments.raw_supported_progressive_link_projection is not None
                or arguments.raw_confirmed_progressive_links
                or arguments.paired_progressive_scaffold_audit
                or arguments.paired_progressive_scaffold_projection is not None
                or arguments.adaptive_support_rescue_audit
                or arguments.adaptive_support_rescue_projection is not None
                or arguments.raw_confirmed_adaptive_rescue_links
                or arguments.adaptive_rescue_chains
                or arguments.selective_support_rescue_audit
                or arguments.selective_support_rescue_projection is not None
                or arguments.high_confidence_support_two_rescue_audit
                or arguments.high_confidence_support_two_rescue_projection
                is not None
                or arguments.support_two_master_graph_audit
                or arguments.support_two_master_graph_projection is not None
                or arguments.selective_rescue_master_graph_audit
                or arguments.selective_rescue_master_graph_projection is not None
                or arguments.raw_confirmed_selective_master_graph_audit
                or arguments.raw_confirmed_selective_master_graph_projection
                is not None
            ):
                parser.error(
                    "progressive projections and link audits require "
                    "--progressive-raw-phase-audit"
                )
            if (
                arguments.raw_confirmed_master_graph_audit
                and not arguments.master_overlap_graph_audit
            ):
                parser.error(
                    "--raw-confirmed-master-graph-audit requires "
                    "--master-overlap-graph-audit"
                )
            if (
                arguments.raw_confirmed_master_graph_projection is not None
                and not arguments.raw_confirmed_master_graph_audit
            ):
                parser.error(
                    "raw-confirmed master graph projection requires "
                    "--raw-confirmed-master-graph-audit"
                )
            if (
                arguments.strain_safe_containment_audit
                and not arguments.master_overlap_graph_audit
            ):
                parser.error(
                    "--strain-safe-containment-audit requires "
                    "--master-overlap-graph-audit"
                )
            if (
                arguments.strain_safe_containment_projection is not None
                and not arguments.strain_safe_containment_audit
            ):
                parser.error(
                    "strain-safe containment projection requires "
                    "--strain-safe-containment-audit"
                )
            if arguments.iterative_reclustering_audit:
                if arguments.max_contig_iterations:
                    parser.error(
                        "--iterative-reclustering-audit requires "
                        "--max-contig-iterations 0"
                    )
                if not (
                    arguments.ranked_extension
                    and arguments.reciprocal_best_extension
                    and arguments.damage_aware_ranking
                ):
                    parser.error(
                        "--iterative-reclustering-audit requires ranked, "
                        "reciprocal-best, and damage-aware extension"
                    )
                (
                    reclustered_contigs,
                    reclustering_audit,
                ) = audit_iterative_reclustering(
                    reads,
                    anchor_k=arguments.anchor_k,
                    anchors_per_read=arguments.anchors_per_read,
                    maximum_anchor_occurrences=arguments.max_anchor_occurrences,
                    minimum_anchor_matches=arguments.min_anchor_matches,
                    minimum_overlap=arguments.min_overlap,
                    minimum_cluster_size=arguments.min_cluster_size,
                    minimum_ranked_extension_support=(
                        arguments.min_ranked_extension_support
                    ),
                    reciprocal_best_extension=True,
                    damage_aware_ranking=True,
                    minimum_confidence_margin=(
                        arguments.min_overlap_confidence_margin
                    ),
                    damage_mismatch_penalty=arguments.damage_mismatch_penalty,
                    ranking_damage_end_window=(
                        arguments.ranking_damage_end_window
                    ),
                    extension_consensus=arguments.extension_consensus,
                    maximum_iterations=arguments.max_recluster_iterations,
                    minimum_output_length=arguments.min_output_length,
                    derived_minimum_identity=(
                        arguments.recluster_derived_min_identity
                    ),
                    derived_minimum_raw_support=(
                        arguments.recluster_derived_min_raw_support
                    ),
                    mismatch_minimum_raw_support=(
                        arguments.recluster_mismatch_min_raw_support
                    ),
                    mismatch_minimum_raw_margin=(
                        arguments.recluster_mismatch_min_raw_margin
                    ),
                    mismatch_minimum_base_quality=(
                        arguments.recluster_mismatch_min_base_quality
                    ),
                )
                if arguments.iterative_reclustering_report is not None:
                    write_iterative_reclustering_report(
                        reclustering_audit,
                        arguments.iterative_reclustering_report,
                    )
                if arguments.iterative_reclustering_projection is not None:
                    write_fasta(
                        [contig.sequence for contig in reclustered_contigs],
                        arguments.iterative_reclustering_projection,
                    )
                if arguments.master_overlap_graph_audit:
                    master_contigs, master_graph_audit = audit_master_overlap_graph(
                        reclustered_contigs,
                        anchor_k=arguments.anchor_k,
                        anchors_per_read=arguments.anchors_per_read,
                        maximum_anchor_occurrences=(
                            arguments.max_anchor_occurrences
                        ),
                        minimum_anchor_matches=arguments.min_anchor_matches,
                        minimum_overlap=arguments.min_overlap,
                    )
                    if arguments.master_overlap_graph_projection is not None:
                        write_fasta(
                            [contig.sequence for contig in master_contigs],
                            arguments.master_overlap_graph_projection,
                        )
                    if arguments.raw_confirmed_master_graph_audit:
                        (
                            raw_master_contigs,
                            raw_master_graph_audit,
                        ) = audit_raw_confirmed_master_overlap_graph(
                            master_contigs,
                            reads,
                            anchor_k=arguments.anchor_k,
                            anchors_per_read=arguments.anchors_per_read,
                            maximum_anchor_occurrences=(
                                arguments.max_anchor_occurrences
                            ),
                            minimum_anchor_matches=(
                                arguments.min_anchor_matches
                            ),
                            minimum_overlap=arguments.min_overlap,
                            minimum_identity=(
                                arguments.raw_confirmed_master_min_identity
                            ),
                            maximum_mismatches=(
                                arguments.raw_confirmed_master_max_mismatches
                            ),
                            minimum_primary_allele_support=(
                                arguments.raw_confirmed_master_primary_support
                            ),
                            minimum_alternate_allele_support=(
                                arguments.raw_confirmed_master_alternate_support
                            ),
                            minimum_support_margin=(
                                arguments.raw_confirmed_master_support_margin
                            ),
                            minimum_base_quality=(
                                arguments.raw_confirmed_master_min_base_quality
                            ),
                            damage_end_window=(
                                arguments.raw_confirmed_master_damage_end_window
                            ),
                        )
                        if (
                            arguments.raw_confirmed_master_graph_projection
                            is not None
                        ):
                            write_fasta(
                                [
                                    contig.sequence
                                    for contig in raw_master_contigs
                                ],
                                arguments.raw_confirmed_master_graph_projection,
                            )
                    if arguments.strain_safe_containment_audit:
                        (
                            containment_contigs,
                            containment_audit,
                        ) = audit_strain_safe_containment(
                            master_contigs,
                            reads,
                            anchor_k=arguments.anchor_k,
                            anchors_per_read=arguments.anchors_per_read,
                            maximum_anchor_occurrences=(
                                arguments.max_anchor_occurrences
                            ),
                            minimum_anchor_matches=arguments.min_anchor_matches,
                            minimum_identity=arguments.containment_min_identity,
                            minimum_coverage=arguments.containment_min_coverage,
                            minimum_primary_allele_support=(
                                arguments.containment_primary_allele_support
                            ),
                            minimum_alternate_allele_support=(
                                arguments.containment_alternate_allele_support
                            ),
                            minimum_support_margin=(
                                arguments.containment_support_margin
                            ),
                            minimum_base_quality=(
                                arguments.containment_min_base_quality
                            ),
                            damage_end_window=(
                                arguments.containment_damage_end_window
                            ),
                        )
                        if (
                            arguments.strain_safe_containment_projection
                            is not None
                        ):
                            write_fasta(
                                [
                                    contig.sequence
                                    for contig in containment_contigs
                                ],
                                arguments.strain_safe_containment_projection,
                            )
            elif (
                arguments.iterative_reclustering_report is not None
                or arguments.iterative_reclustering_projection is not None
                or arguments.master_overlap_graph_audit
                or arguments.master_overlap_graph_projection is not None
                or arguments.raw_confirmed_master_graph_audit
                or arguments.raw_confirmed_master_graph_projection is not None
                or arguments.strain_safe_containment_audit
                or arguments.strain_safe_containment_projection is not None
            ):
                parser.error(
                    "iterative reclustering and master-graph projections require "
                    "--iterative-reclustering-audit"
                )
            if (
                arguments.master_overlap_graph_projection is not None
                and not arguments.master_overlap_graph_audit
            ):
                parser.error(
                    "master overlap graph projection requires "
                    "--master-overlap-graph-audit"
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
            print(f"overlap_confidence_ranked_sides={summary.confidence_ranked_sides}")
            print(
                "overlap_confidence_changed_winners="
                f"{summary.confidence_changed_winners}"
            )
            print(
                "overlap_confidence_ambiguous_extensions="
                f"{summary.confidence_ambiguous_extensions}"
            )
            print(f"overlap_deferred_reads={summary.deferred_reads}")
            print(f"overlap_cross_cluster_reads={summary.cross_cluster_reads}")
            print(
                "overlap_deferred_candidate_reads="
                f"{summary.deferred_candidate_reads}"
            )
            print(
                "overlap_cross_cluster_candidate_reads="
                f"{summary.cross_cluster_candidate_reads}"
            )
            print(
                "overlap_deferred_extension_candidates="
                f"{summary.deferred_extension_candidates}"
            )
            print(
                "overlap_cross_cluster_extension_candidates="
                f"{summary.cross_cluster_extension_candidates}"
            )
            print(
                "overlap_deferred_ambiguous_assignments="
                f"{summary.deferred_ambiguous_assignments}"
            )
            print(
                "overlap_cross_cluster_ambiguous_assignments="
                f"{summary.cross_cluster_ambiguous_assignments}"
            )
            print(
                "overlap_deferred_unique_assignments="
                f"{summary.deferred_unique_assignments}"
            )
            print(
                "overlap_cross_cluster_unique_assignments="
                f"{summary.cross_cluster_unique_assignments}"
            )
            print(
                "overlap_deferred_reciprocal_assignments="
                f"{summary.deferred_reciprocal_assignments}"
            )
            print(
                "overlap_cross_cluster_reciprocal_assignments="
                f"{summary.cross_cluster_reciprocal_assignments}"
            )
            print(
                "overlap_recruitment_supported_contigs="
                f"{summary.recruitment_supported_contigs}"
            )
            print(
                "overlap_recruitment_supported_bases="
                f"{summary.recruitment_supported_bases}"
            )
            print(
                "overlap_deferred_supported_contigs="
                f"{summary.deferred_supported_contigs}"
            )
            print(
                "overlap_deferred_supported_bases="
                f"{summary.deferred_supported_bases}"
            )
            print(
                "overlap_cross_cluster_supported_contigs="
                f"{summary.cross_cluster_supported_contigs}"
            )
            print(
                "overlap_cross_cluster_supported_bases="
                f"{summary.cross_cluster_supported_bases}"
            )
            print(
                "overlap_contig_link_candidate_overlaps="
                f"{summary.contig_link_candidate_overlaps}"
            )
            print(
                "overlap_contig_link_unique_overlaps="
                f"{summary.contig_link_unique_overlaps}"
            )
            print(
                "overlap_contig_link_reciprocal_overlaps="
                f"{summary.contig_link_reciprocal_overlaps}"
            )
            print(
                "overlap_contig_link_read_supported_overlaps="
                f"{summary.contig_link_read_supported_overlaps}"
            )
            print(
                "overlap_contig_link_support_at_least_1="
                f"{summary.contig_link_support_at_least_1}"
            )
            print(
                "overlap_contig_link_support_at_least_2="
                f"{summary.contig_link_support_at_least_2}"
            )
            print(
                "overlap_contig_link_support_at_least_3="
                f"{summary.contig_link_support_at_least_3}"
            )
            print(
                "overlap_contig_link_support_at_least_5="
                f"{summary.contig_link_support_at_least_5}"
            )
            print(
                "overlap_contig_link_max_read_support="
                f"{summary.contig_link_max_read_support}"
            )
            print(
                "overlap_contig_link_ambiguous_ends="
                f"{summary.contig_link_ambiguous_ends}"
            )
            print(
                "overlap_contig_link_cyclic_components="
                f"{summary.contig_link_cyclic_components}"
            )
            print(
                "overlap_contig_link_linear_chains="
                f"{summary.contig_link_linear_chains}"
            )
            print(
                "overlap_contig_link_projected_joins="
                f"{summary.contig_link_projected_joins}"
            )
            print(
                "overlap_contig_link_projected_merged_bases="
                f"{summary.contig_link_projected_merged_bases}"
            )
            print(
                "overlap_contig_link_projected_longest_contig="
                f"{summary.contig_link_projected_longest_contig}"
            )
            print(
                "overlap_two_tier_redundancy_audit="
                f"{str(arguments.two_tier_redundancy_audit).lower()}"
            )
            print(
                "overlap_rescue_min_cluster_size="
                f"{arguments.rescue_min_cluster_size}"
            )
            print("overlap_rescue_min_identity=0.97")
            print("overlap_rescue_min_ry_identity=0.99")
            print("overlap_rescue_min_coverage=0.99")
            print(f"overlap_rescue_contigs={rescue_audit.rescue_contigs}")
            print(f"overlap_rescue_bases={rescue_audit.rescue_bases}")
            print(
                "overlap_rescue_contained_by_primary="
                f"{rescue_audit.contained_by_primary}"
            )
            print(
                "overlap_rescue_redundant_with_rescue="
                f"{rescue_audit.redundant_with_rescue}"
            )
            print(
                "overlap_rescue_extends_primary="
                f"{rescue_audit.extends_primary}"
            )
            print(
                "overlap_rescue_ambiguous_primary_extensions="
                f"{rescue_audit.ambiguous_primary_extensions}"
            )
            print(
                "overlap_rescue_replacement_contigs="
                f"{rescue_audit.replacement_contigs}"
            )
            print(
                "overlap_rescue_replacement_added_bases="
                f"{rescue_audit.replacement_added_bases}"
            )
            print(
                "overlap_rescue_replacement_projected_bases="
                f"{rescue_audit.replacement_projected_bases}"
            )
            print(
                "overlap_rescue_replacement_projected_n50="
                f"{rescue_audit.replacement_projected_n50}"
            )
            print(
                "overlap_rescue_replacement_projected_longest_contig="
                f"{rescue_audit.replacement_projected_longest_contig}"
            )
            print(
                f"overlap_rescue_novel_contigs={rescue_audit.novel_contigs}"
            )
            print(f"overlap_rescue_novel_bases={rescue_audit.novel_bases}")
            print(
                "overlap_rescue_projected_contigs="
                f"{rescue_audit.projected_contigs}"
            )
            print(
                "overlap_rescue_projected_bases="
                f"{rescue_audit.projected_bases}"
            )
            print(
                f"overlap_rescue_projected_n50={rescue_audit.projected_n50}"
            )
            print(
                "overlap_rescue_projected_longest_contig="
                f"{rescue_audit.projected_longest_contig}"
            )
            print(
                "overlap_iterative_reclustering_audit="
                f"{str(arguments.iterative_reclustering_audit).lower()}"
            )
            print(
                "overlap_recluster_report="
                f"{arguments.iterative_reclustering_report or ''}"
            )
            print(
                "overlap_recluster_projection="
                f"{arguments.iterative_reclustering_projection or ''}"
            )
            print(f"overlap_recluster_iterations={reclustering_audit.iterations}")
            print(
                "overlap_recluster_converged="
                f"{str(reclustering_audit.converged).lower()}"
            )
            print(
                "overlap_recluster_rollback_rejected="
                f"{str(reclustering_audit.rollback_rejected).lower()}"
            )
            print(
                "overlap_recluster_extended_sequences="
                f"{reclustering_audit.extended_sequences}"
            )
            print(
                f"overlap_recluster_added_bases={reclustering_audit.added_bases}"
            )
            print(
                "overlap_recluster_candidate_extensions="
                f"{reclustering_audit.candidate_extensions}"
            )
            print(
                "overlap_recluster_reused_candidates="
                f"{reclustering_audit.reused_candidates}"
            )
            print(
                "overlap_recluster_conflicting_candidates="
                f"{reclustering_audit.conflicting_candidates}"
            )
            print(
                "overlap_recluster_reciprocal_checks="
                f"{reclustering_audit.reciprocal_checks}"
            )
            print(
                "overlap_recluster_reciprocal_rejections="
                f"{reclustering_audit.reciprocal_rejections}"
            )
            print(
                "overlap_recluster_derived_overlap_rejections="
                f"{reclustering_audit.derived_overlap_rejections}"
            )
            print(
                "overlap_recluster_raw_confirmed_mismatch_overlaps="
                f"{reclustering_audit.raw_confirmed_mismatch_overlaps}"
            )
            print(
                "overlap_recluster_raw_support_rejections="
                f"{reclustering_audit.raw_support_rejections}"
            )
            print(
                "overlap_recluster_assembled_sequences="
                f"{reclustering_audit.assembled_sequences}"
            )
            print(
                "overlap_recluster_redundant_sequences="
                f"{reclustering_audit.redundant_sequences}"
            )
            print(
                "overlap_recluster_projected_contigs="
                f"{reclustering_audit.projected_contigs}"
            )
            print(
                "overlap_recluster_projected_bases="
                f"{reclustering_audit.projected_bases}"
            )
            print(
                "overlap_recluster_projected_n50="
                f"{reclustering_audit.projected_n50}"
            )
            print(
                "overlap_recluster_projected_longest_contig="
                f"{reclustering_audit.projected_longest_contig}"
            )
            print(
                "overlap_progressive_raw_phase_audit="
                f"{str(arguments.progressive_raw_phase_audit).lower()}"
            )
            print(
                "overlap_progressive_raw_phase_projection="
                f"{arguments.progressive_raw_phase_projection or ''}"
            )
            print(
                "overlap_progressive_candidate_centers="
                f"{progressive_clustering.candidate_centers}"
            )
            print(
                "overlap_progressive_clusters="
                f"{progressive_clustering.clusters}"
            )
            print(
                "overlap_progressive_clusters_below_minimum_size="
                f"{progressive_clustering.clusters_below_minimum_size}"
            )
            print(
                "overlap_progressive_clustered_reads="
                f"{progressive_clustering.clustered_reads}"
            )
            print(
                "overlap_progressive_retained_unclustered_reads="
                f"{progressive_clustering.retained_unclustered_reads}"
            )
            print(
                "overlap_progressive_extended_centers="
                f"{progressive_extension.extended_centers}"
            )
            print(
                "overlap_progressive_extension_rounds="
                f"{progressive_extension.extension_rounds}"
            )
            print(
                "overlap_progressive_recruited_reads="
                f"{progressive_extension.recruited_reads}"
            )
            print(
                "overlap_progressive_added_bases="
                f"{progressive_extension.added_bases}"
            )
            print(
                "overlap_progressive_consumed_reads="
                f"{progressive_extension.consumed_reads}"
            )
            print(
                "overlap_progressive_ambiguous_extensions="
                f"{progressive_extension.ambiguous_extensions}"
            )
            print(
                "overlap_progressive_reciprocal_checks="
                f"{progressive_extension.reciprocal_checks}"
            )
            print(
                "overlap_progressive_reciprocal_rejections="
                f"{progressive_extension.reciprocal_rejections}"
            )
            for reason, count in progressive_extension.round_outcomes:
                print(f"overlap_progressive_round_{reason}={count}")
            for reason, count in progressive_extension.recruit_search:
                print(f"overlap_progressive_recruit_search_{reason}={count}")
            print(f"overlap_progressive_low_quality_ry_rescue={str(arguments.progressive_low_quality_ry_rescue).lower()}")
            for reason, count in progressive_iterations.recruit_search:
                print(f"overlap_progressive_iteration_search_{reason}={count}")
            print(
                "overlap_progressive_iterations="
                f"{progressive_iterations.iterations}"
            )
            print(
                "overlap_progressive_converged="
                f"{str(progressive_iterations.converged).lower()}"
            )
            print(
                "overlap_progressive_iteration_candidates="
                f"{progressive_iterations.candidate_alignments}"
            )
            print(
                "overlap_progressive_iteration_extended_centers="
                f"{progressive_iterations.extended_centers}"
            )
            print(
                "overlap_progressive_iteration_added_bases="
                f"{progressive_iterations.added_bases}"
            )
            print(
                "overlap_progressive_iteration_consumed_reads="
                f"{progressive_iterations.consumed_reads}"
            )
            print(
                "overlap_progressive_iteration_ambiguous_extensions="
                f"{progressive_iterations.ambiguous_extensions}"
            )
            print(
                "overlap_progressive_iteration_insufficient_consensus_support="
                f"{progressive_iterations.insufficient_consensus_support}"
            )
            print(
                "overlap_progressive_evidence_priority="
                f"{str(arguments.evidence_priority_progressive_extension).lower()}"
            )
            print(
                "overlap_progressive_evidence_priority_sweeps="
                f"{progressive_iterations.evidence_priority_sweeps}"
            )
            print(
                "overlap_progressive_evidence_priority_reordered_centers="
                f"{progressive_iterations.evidence_priority_reordered_centers}"
            )
            print(
                "overlap_progressive_evidence_priority_claim_conflicts="
                f"{progressive_iterations.evidence_priority_claim_conflicts}"
            )
            print(
                "overlap_progressive_failure_audit="
                f"{str(arguments.progressive_extension_failure_audit).lower()}"
            )
            for label, value in (
                ("left_sides", progressive_iterations.rejected_left_sides),
                ("right_sides", progressive_iterations.rejected_right_sides),
                ("support_1", progressive_iterations.rejected_support_1),
                ("support_2", progressive_iterations.rejected_support_2),
                ("support_3", progressive_iterations.rejected_support_3),
                ("support_4", progressive_iterations.rejected_support_4),
                ("support_5_plus", progressive_iterations.rejected_support_5_plus),
                ("agreeing_sides", progressive_iterations.rejected_agreeing_sides),
                (
                    "conflicting_sides",
                    progressive_iterations.rejected_conflicting_sides,
                ),
                ("extension_1_5", progressive_iterations.rejected_extension_1_5),
                ("extension_6_10", progressive_iterations.rejected_extension_6_10),
                (
                    "extension_11_20",
                    progressive_iterations.rejected_extension_11_20,
                ),
                (
                    "extension_21_plus",
                    progressive_iterations.rejected_extension_21_plus,
                ),
                (
                    "high_quality_boundary_observations",
                    progressive_iterations.rejected_high_quality_boundary_observations,
                ),
                (
                    "low_quality_boundary_observations",
                    progressive_iterations.rejected_low_quality_boundary_observations,
                ),
                (
                    "missing_quality_boundary_observations",
                    progressive_iterations.rejected_missing_quality_boundary_observations,
                ),
                (
                    "damage_compatible_mismatches",
                    progressive_iterations.rejected_damage_compatible_mismatches,
                ),
                (
                    "ordinary_mismatches",
                    progressive_iterations.rejected_ordinary_mismatches,
                ),
            ):
                print(f"overlap_progressive_rejected_{label}={value}")
            print(
                "overlap_support_three_progressive_audit="
                f"{str(arguments.support_three_progressive_extension_audit).lower()}"
            )
            print(
                "overlap_support_three_progressive_projection="
                f"{arguments.support_three_progressive_projection or ''}"
            )
            for label, value in (
                ("iterations", support_three_iterations.iterations),
                ("extended_centers", support_three_iterations.extended_centers),
                ("added_bases", support_three_iterations.added_bases),
                (
                    "accepted_sides",
                    support_three_iterations.strict_boundary_accepted_sides,
                ),
                (
                    "conflict_rejections",
                    support_three_iterations.strict_boundary_conflict_rejections,
                ),
                (
                    "quality_rejections",
                    support_three_iterations.strict_boundary_quality_rejections,
                ),
                ("projected_contigs", support_three_iterations.projected_contigs),
                ("projected_bases", support_three_iterations.projected_bases),
                ("projected_n50", support_three_iterations.projected_n50),
                (
                    "projected_longest_contig",
                    support_three_iterations.projected_longest_contig,
                ),
            ):
                print(f"overlap_support_three_progressive_{label}={value}")
            print(
                "overlap_adaptive_rescue_audit="
                f"{str(arguments.adaptive_support_rescue_audit).lower()}"
            )
            print(
                "overlap_adaptive_rescue_projection="
                f"{arguments.adaptive_support_rescue_projection or ''}"
            )
            print(
                "overlap_adaptive_rescue_primary_contigs="
                f"{adaptive_rescue.primary_contigs}"
            )
            print(
                "overlap_adaptive_rescue_eligible_raw_reads="
                f"{adaptive_rescue.eligible_raw_reads}"
            )
            print(
                "overlap_adaptive_rescue_clusters="
                f"{adaptive_rescue.rescue_clusters}"
            )
            print(
                "overlap_adaptive_rescue_clustered_reads="
                f"{adaptive_rescue.rescue_clustered_reads}"
            )
            print(
                "overlap_adaptive_rescue_contigs="
                f"{adaptive_rescue.rescue_contigs}"
            )
            print(
                "overlap_adaptive_rescue_candidate_primary_links="
                f"{adaptive_rescue.candidate_primary_links}"
            )
            print(
                "overlap_adaptive_rescue_exact_primary_links="
                f"{adaptive_rescue.exact_primary_links}"
            )
            print(
                "overlap_adaptive_rescue_chains_enabled="
                f"{str(arguments.adaptive_rescue_chains).lower()}"
            )
            print(
                "overlap_adaptive_rescue_candidate_rescue_links="
                f"{adaptive_rescue.candidate_rescue_links}"
            )
            print(
                "overlap_adaptive_rescue_exact_rescue_links="
                f"{adaptive_rescue.exact_rescue_links}"
            )
            print(
                "overlap_adaptive_rescue_near_exact_enabled="
                f"{str(arguments.raw_confirmed_adaptive_rescue_links).lower()}"
            )
            print(
                "overlap_adaptive_rescue_near_exact_candidates="
                f"{adaptive_rescue.near_exact_candidates}"
            )
            print(
                "overlap_adaptive_rescue_near_exact_raw_confirmed="
                f"{adaptive_rescue.near_exact_raw_confirmed}"
            )
            print(
                "overlap_adaptive_rescue_near_exact_accepted="
                f"{adaptive_rescue.near_exact_accepted}"
            )
            print(
                "overlap_adaptive_rescue_mismatch_positions="
                f"{adaptive_rescue.mismatch_positions}"
            )
            print(
                "overlap_adaptive_rescue_insufficient_raw_support="
                f"{adaptive_rescue.insufficient_raw_support}"
            )
            print(
                "overlap_adaptive_rescue_strain_conflicts="
                f"{adaptive_rescue.strain_conflicts}"
            )
            print(
                "overlap_adaptive_rescue_ambiguous_near_molecules="
                f"{adaptive_rescue.ambiguous_near_molecules}"
            )
            print(
                "overlap_adaptive_rescue_near_unique_support_rejections="
                f"{adaptive_rescue.near_exact_unique_support_rejections}"
            )
            print(
                "overlap_adaptive_rescue_exact_preferred="
                f"{adaptive_rescue.exact_preferred}"
            )
            print(
                "overlap_adaptive_rescue_anchored_components="
                f"{adaptive_rescue.anchored_components}"
            )
            print(
                "overlap_adaptive_rescue_rejected_unanchored_components="
                f"{adaptive_rescue.rejected_unanchored_components}"
            )
            print(
                "overlap_adaptive_rescue_promoted_contigs="
                f"{adaptive_rescue.promoted_rescue_contigs}"
            )
            print(
                "overlap_adaptive_rescue_rejected_unattached="
                f"{adaptive_rescue.rejected_unattached_contigs}"
            )
            print(
                "overlap_adaptive_rescue_ambiguous_ends="
                f"{adaptive_rescue.ambiguous_ends}"
            )
            print(
                "overlap_adaptive_rescue_reciprocal_edges="
                f"{adaptive_rescue.reciprocal_edges}"
            )
            print(
                "overlap_adaptive_rescue_corrected_overlap_bases="
                f"{adaptive_rescue.corrected_overlap_bases}"
            )
            print(
                "overlap_adaptive_rescue_projected_contigs="
                f"{adaptive_rescue.projected_contigs}"
            )
            print(
                "overlap_adaptive_rescue_projected_bases="
                f"{adaptive_rescue.projected_bases}"
            )
            print(
                "overlap_adaptive_rescue_projected_n50="
                f"{adaptive_rescue.projected_n50}"
            )
            print(
                "overlap_adaptive_rescue_projected_longest_contig="
                f"{adaptive_rescue.projected_longest_contig}"
            )
            print(
                "overlap_selective_rescue_audit="
                f"{str(arguments.selective_support_rescue_audit).lower()}"
            )
            print(
                "overlap_selective_rescue_projection="
                f"{arguments.selective_support_rescue_projection or ''}"
            )
            for label, value in (
                ("primary_contigs", selective_rescue.primary_contigs),
                ("eligible_raw_reads", selective_rescue.eligible_raw_reads),
                ("clusters", selective_rescue.rescue_clusters),
                ("clustered_reads", selective_rescue.rescue_clustered_reads),
                ("candidate_contigs", selective_rescue.rescue_contigs),
                (
                    "insufficient_molecule_support",
                    selective_rescue.insufficient_molecule_support,
                ),
                ("contained_by_primary", selective_rescue.contained_by_primary),
                (
                    "redundant_with_rescue",
                    selective_rescue.redundant_with_rescue,
                ),
                ("primary_extensions", selective_rescue.primary_extensions),
                (
                    "ambiguous_primary_extensions",
                    selective_rescue.ambiguous_primary_extensions,
                ),
                ("novel_contigs", selective_rescue.novel_contigs),
                ("novel_bases", selective_rescue.novel_bases),
                ("projected_contigs", selective_rescue.projected_contigs),
                ("projected_bases", selective_rescue.projected_bases),
                ("projected_n50", selective_rescue.projected_n50),
                (
                    "projected_longest_contig",
                    selective_rescue.projected_longest_contig,
                ),
            ):
                print(f"overlap_selective_rescue_{label}={value}")
            print(
                "overlap_support_two_rescue_audit="
                f"{str(arguments.high_confidence_support_two_rescue_audit).lower()}"
            )
            print(
                "overlap_support_two_rescue_projection="
                f"{arguments.high_confidence_support_two_rescue_projection or ''}"
            )
            for label, value in (
                ("eligible_raw_reads", support_two_rescue.eligible_raw_reads),
                ("seed_quality_filter", str(arguments.support_two_seed_quality_filter).lower()),
                ("trim_seed_tails", str(arguments.support_two_trim_seed_tails).lower()),
                ("seed_trimmed_contigs", support_two_rescue.seed_trimmed_contigs),
                ("seed_trimmed_left_bases", support_two_rescue.seed_trimmed_left_bases),
                ("seed_trimmed_right_bases", support_two_rescue.seed_trimmed_right_bases),
                ("trimmed_short_contigs", support_two_rescue.trimmed_short_contigs),
                ("fixed_membership", str(arguments.support_two_fixed_membership_before is not None).lower()),
                ("fixed_membership_short_exclusions", support_two_rescue.fixed_membership_short_exclusions),
                ("primary_extension_exclusions", support_two_rescue.primary_extension_exclusions),
                ("ambiguous_primary_extension_exclusions", support_two_rescue.ambiguous_primary_extension_exclusions),
                ("recruit_quality_filter", str(arguments.support_two_recruit_quality_filter).lower()),
                ("seed_quality_rejections", support_two_rescue.seed_quality_rejections),
                ("seed_quality_rejected_bases", support_two_rescue.seed_quality_rejected_bases),
                ("recruit_quality_rejections", support_two_rescue.recruit_quality_rejections),
                ("recruit_quality_rejected_overhang_bases", support_two_rescue.recruit_quality_rejected_overhang_bases),
                ("candidate_clusters", support_two_rescue.candidate_clusters),
                ("admitted_clusters", support_two_rescue.admitted_clusters),
                (
                    "ordinary_mismatch_rejections",
                    support_two_rescue.ordinary_mismatch_rejections,
                ),
                (
                    "internal_damage_rejections",
                    support_two_rescue.internal_damage_rejections,
                ),
                (
                    "boundary_quality_rejections",
                    support_two_rescue.boundary_quality_rejections,
                ),
                ("rescue_contigs", support_two_rescue.rescue_contigs),
                ("contained_by_primary", support_two_rescue.contained_by_primary),
                (
                    "redundant_with_rescue",
                    support_two_rescue.redundant_with_rescue,
                ),
                ("novel_contigs", support_two_rescue.novel_contigs),
                ("novel_bases", support_two_rescue.novel_bases),
                ("projected_contigs", support_two_rescue.projected_contigs),
                ("projected_bases", support_two_rescue.projected_bases),
                ("projected_n50", support_two_rescue.projected_n50),
                (
                    "projected_longest_contig",
                    support_two_rescue.projected_longest_contig,
                ),
            ):
                print(f"overlap_support_two_rescue_{label}={value}")
            print(
                "overlap_support_two_master_graph_audit="
                f"{str(arguments.support_two_master_graph_audit).lower()}"
            )
            print(
                "overlap_support_two_master_graph_projection="
                f"{arguments.support_two_master_graph_projection or ''}"
            )
            for label, value in (
                ("input_contigs", support_two_master_graph.input_contigs),
                ("candidate_dovetails", support_two_master_graph.candidate_dovetails),
                ("exact_dovetails", support_two_master_graph.exact_dovetails),
                (
                    "transitive_edges_removed",
                    support_two_master_graph.transitive_edges_removed,
                ),
                ("ambiguous_ends", support_two_master_graph.ambiguous_ends),
                ("reciprocal_edges", support_two_master_graph.reciprocal_edges),
                ("linear_paths", support_two_master_graph.linear_paths),
                ("cyclic_components", support_two_master_graph.cyclic_components),
                ("merged_contigs", support_two_master_graph.merged_contigs),
                ("added_bases", support_two_master_graph.added_bases),
                ("projected_contigs", support_two_master_graph.projected_contigs),
                ("projected_bases", support_two_master_graph.projected_bases),
                ("projected_n50", support_two_master_graph.projected_n50),
                (
                    "projected_longest_contig",
                    support_two_master_graph.projected_longest_contig,
                ),
            ):
                print(f"overlap_support_two_master_graph_{label}={value}")
            print(
                "overlap_selective_rescue_master_graph_audit="
                f"{str(arguments.selective_rescue_master_graph_audit).lower()}"
            )
            print(
                "overlap_selective_rescue_master_graph_projection="
                f"{arguments.selective_rescue_master_graph_projection or ''}"
            )
            for label, value in (
                ("input_contigs", selective_rescue_master_graph.input_contigs),
                (
                    "candidate_dovetails",
                    selective_rescue_master_graph.candidate_dovetails,
                ),
                ("exact_dovetails", selective_rescue_master_graph.exact_dovetails),
                (
                    "transitive_edges_removed",
                    selective_rescue_master_graph.transitive_edges_removed,
                ),
                ("ambiguous_ends", selective_rescue_master_graph.ambiguous_ends),
                ("reciprocal_edges", selective_rescue_master_graph.reciprocal_edges),
                ("linear_paths", selective_rescue_master_graph.linear_paths),
                (
                    "cyclic_components",
                    selective_rescue_master_graph.cyclic_components,
                ),
                ("merged_contigs", selective_rescue_master_graph.merged_contigs),
                ("added_bases", selective_rescue_master_graph.added_bases),
                ("projected_contigs", selective_rescue_master_graph.projected_contigs),
                ("projected_bases", selective_rescue_master_graph.projected_bases),
                ("projected_n50", selective_rescue_master_graph.projected_n50),
                (
                    "projected_longest_contig",
                    selective_rescue_master_graph.projected_longest_contig,
                ),
            ):
                print(f"overlap_selective_rescue_master_graph_{label}={value}")
            print(
                "overlap_raw_confirmed_selective_master_graph_audit="
                f"{str(arguments.raw_confirmed_selective_master_graph_audit).lower()}"
            )
            print(
                "overlap_raw_confirmed_selective_master_graph_projection="
                f"{arguments.raw_confirmed_selective_master_graph_projection or ''}"
            )
            for label, value in (
                ("input_contigs", raw_selective_master_graph.input_contigs),
                ("candidate_dovetails", raw_selective_master_graph.candidate_dovetails),
                ("exact_dovetails", raw_selective_master_graph.exact_dovetails),
                (
                    "near_exact_candidates",
                    raw_selective_master_graph.near_exact_candidates,
                ),
                ("near_exact_accepted", raw_selective_master_graph.near_exact_accepted),
                ("mismatch_positions", raw_selective_master_graph.mismatch_positions),
                (
                    "insufficient_raw_support",
                    raw_selective_master_graph.insufficient_raw_support,
                ),
                ("strain_conflicts", raw_selective_master_graph.strain_conflicts),
                ("exact_preferred", raw_selective_master_graph.exact_preferred),
                (
                    "transitive_edges_removed",
                    raw_selective_master_graph.transitive_edges_removed,
                ),
                ("ambiguous_ends", raw_selective_master_graph.ambiguous_ends),
                ("reciprocal_edges", raw_selective_master_graph.reciprocal_edges),
                ("linear_paths", raw_selective_master_graph.linear_paths),
                ("cyclic_components", raw_selective_master_graph.cyclic_components),
                ("merged_contigs", raw_selective_master_graph.merged_contigs),
                (
                    "corrected_overlap_bases",
                    raw_selective_master_graph.corrected_overlap_bases,
                ),
                ("added_bases", raw_selective_master_graph.added_bases),
                ("projected_contigs", raw_selective_master_graph.projected_contigs),
                ("projected_bases", raw_selective_master_graph.projected_bases),
                ("projected_n50", raw_selective_master_graph.projected_n50),
                (
                    "projected_longest_contig",
                    raw_selective_master_graph.projected_longest_contig,
                ),
            ):
                print(
                    "overlap_raw_confirmed_selective_master_graph_"
                    f"{label}={value}"
                )
            print(
                "overlap_progressive_link_audit="
                f"{str(arguments.raw_supported_progressive_link_audit).lower()}"
            )
            print(
                "overlap_progressive_link_projection="
                f"{arguments.raw_supported_progressive_link_projection or ''}"
            )
            print(
                "overlap_progressive_link_input_contigs="
                f"{progressive_links.input_contigs}"
            )
            print(
                "overlap_progressive_link_candidate_dovetails="
                f"{progressive_links.candidate_dovetails}"
            )
            print(
                "overlap_progressive_link_exact_dovetails="
                f"{progressive_links.exact_dovetails}"
            )
            print(
                "overlap_progressive_link_near_exact_enabled="
                f"{str(arguments.raw_confirmed_progressive_links).lower()}"
            )
            print(
                "overlap_progressive_link_near_exact_candidates="
                f"{progressive_links.near_exact_candidates}"
            )
            print(
                "overlap_progressive_link_near_exact_raw_confirmed="
                f"{progressive_links.near_exact_raw_confirmed}"
            )
            print(
                "overlap_progressive_link_near_exact_accepted="
                f"{progressive_links.near_exact_accepted}"
            )
            print(
                "overlap_progressive_link_mismatch_positions="
                f"{progressive_links.mismatch_positions}"
            )
            print(
                "overlap_progressive_link_insufficient_raw_support="
                f"{progressive_links.insufficient_raw_support}"
            )
            print(
                "overlap_progressive_link_strain_conflicts="
                f"{progressive_links.strain_conflicts}"
            )
            print(
                "overlap_progressive_link_ambiguous_near_molecules="
                f"{progressive_links.ambiguous_near_molecules}"
            )
            print(
                "overlap_progressive_link_near_unique_support_rejections="
                f"{progressive_links.near_exact_unique_support_rejections}"
            )
            print(
                "overlap_progressive_link_exact_preferred="
                f"{progressive_links.exact_preferred}"
            )
            print(
                "overlap_progressive_link_bridge_candidate_molecules="
                f"{progressive_links.bridge_candidate_molecules}"
            )
            print(
                "overlap_progressive_link_ambiguous_bridge_molecules="
                f"{progressive_links.ambiguous_bridge_molecules}"
            )
            print(
                "overlap_progressive_link_support_at_least_1="
                f"{progressive_links.support_at_least_1}"
            )
            print(
                "overlap_progressive_link_support_at_least_2="
                f"{progressive_links.support_at_least_2}"
            )
            print(
                "overlap_progressive_link_support_at_least_3="
                f"{progressive_links.support_at_least_3}"
            )
            print(
                "overlap_progressive_link_support_at_least_5="
                f"{progressive_links.support_at_least_5}"
            )
            print(
                "overlap_progressive_link_max_read_support="
                f"{progressive_links.max_read_support}"
            )
            print(
                "overlap_progressive_link_supported_dovetails="
                f"{progressive_links.supported_dovetails}"
            )
            print(
                "overlap_progressive_link_transitive_edges_removed="
                f"{progressive_links.transitive_edges_removed}"
            )
            print(
                "overlap_progressive_link_ambiguous_ends="
                f"{progressive_links.ambiguous_ends}"
            )
            print(
                "overlap_progressive_link_reciprocal_edges="
                f"{progressive_links.reciprocal_edges}"
            )
            print(
                "overlap_progressive_link_linear_paths="
                f"{progressive_links.linear_paths}"
            )
            print(
                "overlap_progressive_link_cyclic_components="
                f"{progressive_links.cyclic_components}"
            )
            print(
                "overlap_progressive_link_merged_contigs="
                f"{progressive_links.merged_contigs}"
            )
            print(
                "overlap_progressive_link_corrected_overlap_bases="
                f"{progressive_links.corrected_overlap_bases}"
            )
            print(
                "overlap_progressive_link_added_bases="
                f"{progressive_links.added_bases}"
            )
            print(
                "overlap_progressive_link_projected_contigs="
                f"{progressive_links.projected_contigs}"
            )
            print(
                "overlap_progressive_link_projected_bases="
                f"{progressive_links.projected_bases}"
            )
            print(
                "overlap_progressive_link_projected_n50="
                f"{progressive_links.projected_n50}"
            )
            print(
                "overlap_progressive_link_projected_longest_contig="
                f"{progressive_links.projected_longest_contig}"
            )
            print(
                "overlap_paired_scaffold_audit="
                f"{str(arguments.paired_progressive_scaffold_audit).lower()}"
            )
            print(
                "overlap_paired_scaffold_projection="
                f"{arguments.paired_progressive_scaffold_projection or ''}"
            )
            print(
                "overlap_paired_scaffold_input_contigs="
                f"{paired_scaffold.input_contigs}"
            )
            print(
                "overlap_paired_scaffold_molecules="
                f"{paired_scaffold.paired_molecules}"
            )
            print(
                "overlap_paired_scaffold_uniquely_mapped_pairs="
                f"{paired_scaffold.uniquely_mapped_pairs}"
            )
            print(
                "overlap_paired_scaffold_ambiguous_reads="
                f"{paired_scaffold.ambiguous_reads}"
            )
            print(
                "overlap_paired_scaffold_same_contig_pairs="
                f"{paired_scaffold.same_contig_pairs}"
            )
            print(
                "overlap_paired_scaffold_calibration_pairs="
                f"{paired_scaffold.calibration_pairs}"
            )
            print(
                "overlap_paired_scaffold_calibration_valid="
                f"{str(paired_scaffold.calibration_valid).lower()}"
            )
            print(
                "overlap_paired_scaffold_fragment_mean="
                f"{paired_scaffold.fragment_mean:.3f}"
            )
            print(
                "overlap_paired_scaffold_fragment_sd="
                f"{paired_scaffold.fragment_sd:.3f}"
            )
            print(
                "overlap_paired_scaffold_cross_contig_pairs="
                f"{paired_scaffold.cross_contig_pairs}"
            )
            print(
                "overlap_paired_scaffold_negative_gap_pairs="
                f"{paired_scaffold.negative_gap_pairs}"
            )
            print(
                "overlap_paired_scaffold_candidate_links="
                f"{paired_scaffold.candidate_links}"
            )
            print(
                "overlap_paired_scaffold_supported_links="
                f"{paired_scaffold.supported_links}"
            )
            print(
                "overlap_paired_scaffold_ambiguous_ends="
                f"{paired_scaffold.ambiguous_ends}"
            )
            print(
                "overlap_paired_scaffold_reciprocal_links="
                f"{paired_scaffold.reciprocal_links}"
            )
            print(
                "overlap_paired_scaffold_cyclic_components="
                f"{paired_scaffold.cyclic_components}"
            )
            print(
                "overlap_paired_scaffold_linear_paths="
                f"{paired_scaffold.linear_paths}"
            )
            print(
                "overlap_paired_scaffold_joined_contigs="
                f"{paired_scaffold.joined_contigs}"
            )
            print(
                "overlap_paired_scaffold_inserted_gap_bases="
                f"{paired_scaffold.inserted_gap_bases}"
            )
            print(
                "overlap_paired_scaffold_projected_contigs="
                f"{paired_scaffold.projected_contigs}"
            )
            print(
                "overlap_paired_scaffold_projected_bases="
                f"{paired_scaffold.projected_bases}"
            )
            print(
                "overlap_paired_scaffold_projected_n50="
                f"{paired_scaffold.projected_n50}"
            )
            print(
                "overlap_paired_scaffold_projected_longest_contig="
                f"{paired_scaffold.projected_longest_contig}"
            )
            print(
                "overlap_paired_merge_enabled="
                f"{str(arguments.merge_overlapping_pairs).lower()}"
            )
            print(f"overlap_paired_merge_input_pairs={paired_merge.input_pairs}")
            print(f"overlap_paired_merge_merged_pairs={paired_merge.merged_pairs}")
            print(f"overlap_paired_merge_unmerged_pairs={paired_merge.unmerged_pairs}")
            print(f"overlap_paired_merge_ambiguous_pairs={paired_merge.ambiguous_pairs}")
            print(f"overlap_paired_merge_merged_bases={paired_merge.merged_bases}")
            print(
                "overlap_master_graph_audit="
                f"{str(arguments.master_overlap_graph_audit).lower()}"
            )
            print(
                "overlap_master_graph_projection="
                f"{arguments.master_overlap_graph_projection or ''}"
            )
            print(
                "overlap_master_graph_input_contigs="
                f"{master_graph_audit.input_contigs}"
            )
            print(
                "overlap_master_graph_candidate_dovetails="
                f"{master_graph_audit.candidate_dovetails}"
            )
            print(
                "overlap_master_graph_exact_dovetails="
                f"{master_graph_audit.exact_dovetails}"
            )
            print(
                "overlap_master_graph_transitive_edges_removed="
                f"{master_graph_audit.transitive_edges_removed}"
            )
            print(
                "overlap_master_graph_ambiguous_ends="
                f"{master_graph_audit.ambiguous_ends}"
            )
            print(
                "overlap_master_graph_reciprocal_edges="
                f"{master_graph_audit.reciprocal_edges}"
            )
            print(
                "overlap_master_graph_linear_paths="
                f"{master_graph_audit.linear_paths}"
            )
            print(
                "overlap_master_graph_cyclic_components="
                f"{master_graph_audit.cyclic_components}"
            )
            print(
                "overlap_master_graph_merged_contigs="
                f"{master_graph_audit.merged_contigs}"
            )
            print(
                "overlap_master_graph_added_bases="
                f"{master_graph_audit.added_bases}"
            )
            print(
                "overlap_master_graph_projected_contigs="
                f"{master_graph_audit.projected_contigs}"
            )
            print(
                "overlap_master_graph_projected_bases="
                f"{master_graph_audit.projected_bases}"
            )
            print(
                "overlap_master_graph_projected_n50="
                f"{master_graph_audit.projected_n50}"
            )
            print(
                "overlap_master_graph_projected_longest_contig="
                f"{master_graph_audit.projected_longest_contig}"
            )
            print(
                "overlap_raw_confirmed_master_graph_audit="
                f"{str(arguments.raw_confirmed_master_graph_audit).lower()}"
            )
            print(
                "overlap_raw_confirmed_master_graph_projection="
                f"{arguments.raw_confirmed_master_graph_projection or ''}"
            )
            print(
                "overlap_raw_confirmed_master_input_contigs="
                f"{raw_master_graph_audit.input_contigs}"
            )
            print(
                "overlap_raw_confirmed_master_candidate_dovetails="
                f"{raw_master_graph_audit.candidate_dovetails}"
            )
            print(
                "overlap_raw_confirmed_master_exact_dovetails="
                f"{raw_master_graph_audit.exact_dovetails}"
            )
            print(
                "overlap_raw_confirmed_master_near_exact_candidates="
                f"{raw_master_graph_audit.near_exact_candidates}"
            )
            print(
                "overlap_raw_confirmed_master_near_exact_accepted="
                f"{raw_master_graph_audit.near_exact_accepted}"
            )
            print(
                "overlap_raw_confirmed_master_mismatch_positions="
                f"{raw_master_graph_audit.mismatch_positions}"
            )
            print(
                "overlap_raw_confirmed_master_insufficient_raw_support="
                f"{raw_master_graph_audit.insufficient_raw_support}"
            )
            print(
                "overlap_raw_confirmed_master_strain_conflicts="
                f"{raw_master_graph_audit.strain_conflicts}"
            )
            print(
                "overlap_raw_confirmed_master_exact_preferred="
                f"{raw_master_graph_audit.exact_preferred}"
            )
            print(
                "overlap_raw_confirmed_master_transitive_edges_removed="
                f"{raw_master_graph_audit.transitive_edges_removed}"
            )
            print(
                "overlap_raw_confirmed_master_ambiguous_ends="
                f"{raw_master_graph_audit.ambiguous_ends}"
            )
            print(
                "overlap_raw_confirmed_master_reciprocal_edges="
                f"{raw_master_graph_audit.reciprocal_edges}"
            )
            print(
                "overlap_raw_confirmed_master_linear_paths="
                f"{raw_master_graph_audit.linear_paths}"
            )
            print(
                "overlap_raw_confirmed_master_cyclic_components="
                f"{raw_master_graph_audit.cyclic_components}"
            )
            print(
                "overlap_raw_confirmed_master_merged_contigs="
                f"{raw_master_graph_audit.merged_contigs}"
            )
            print(
                "overlap_raw_confirmed_master_corrected_overlap_bases="
                f"{raw_master_graph_audit.corrected_overlap_bases}"
            )
            print(
                "overlap_raw_confirmed_master_added_bases="
                f"{raw_master_graph_audit.added_bases}"
            )
            print(
                "overlap_raw_confirmed_master_projected_contigs="
                f"{raw_master_graph_audit.projected_contigs}"
            )
            print(
                "overlap_raw_confirmed_master_projected_bases="
                f"{raw_master_graph_audit.projected_bases}"
            )
            print(
                "overlap_raw_confirmed_master_projected_n50="
                f"{raw_master_graph_audit.projected_n50}"
            )
            print(
                "overlap_raw_confirmed_master_projected_longest_contig="
                f"{raw_master_graph_audit.projected_longest_contig}"
            )
            print(
                "overlap_strain_safe_containment_audit="
                f"{str(arguments.strain_safe_containment_audit).lower()}"
            )
            print(
                "overlap_containment_projection="
                f"{arguments.strain_safe_containment_projection or ''}"
            )
            print(
                "overlap_containment_input_contigs="
                f"{containment_audit.input_contigs}"
            )
            print(
                "overlap_containment_candidate_containments="
                f"{containment_audit.candidate_containments}"
            )
            print(
                "overlap_containment_exact_duplicates="
                f"{containment_audit.exact_duplicates}"
            )
            print(
                "overlap_containment_near_duplicate_candidates="
                f"{containment_audit.near_duplicate_candidates}"
            )
            print(
                "overlap_containment_strain_protected="
                f"{containment_audit.strain_protected}"
            )
            print(
                "overlap_containment_insufficient_evidence="
                f"{containment_audit.insufficient_evidence}"
            )
            print(
                "overlap_containment_unsupported_variants="
                f"{containment_audit.unsupported_variants}"
            )
            print(
                "overlap_containment_removed_contigs="
                f"{containment_audit.removed_contigs}"
            )
            print(
                "overlap_containment_removed_bases="
                f"{containment_audit.removed_bases}"
            )
            print(
                "overlap_containment_projected_contigs="
                f"{containment_audit.projected_contigs}"
            )
            print(
                "overlap_containment_projected_bases="
                f"{containment_audit.projected_bases}"
            )
            print(
                "overlap_containment_projected_n50="
                f"{containment_audit.projected_n50}"
            )
            print(
                "overlap_containment_projected_longest_contig="
                f"{containment_audit.projected_longest_contig}"
            )
            print(f"output={arguments.output}")
            return 0
    except (OSError, ValueError) as error:
        parser.error(str(error))

    parser.error(f"unknown command: {arguments.command}")
