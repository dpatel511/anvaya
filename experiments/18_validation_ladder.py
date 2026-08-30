"""Run the frozen, truth-labelled Anvaya validation ladder."""

import argparse
import csv
import hashlib
import json
import platform
import random
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.dead_end_attribution import write_dead_end_attribution_report
from anvaya.damage_consensus import apply_damage_aware_consensus
from anvaya.read_correction import write_correction_report
from anvaya.reads import Read
from anvaya.sequences import canonical_sequence, reverse_complement
from anvaya.unitig_graph import CompactedUnitigGraph, build_compacted_unitig_graph
from helpers.simulation import (
    SimulatedRead,
    add_sequencing_errors,
    add_terminal_damage,
    simulate_fragments,
)


DEFAULT_CONFIG = Path(__file__).with_name("benchmark_ladder.json")
SCHEMA_VERSION = 1


@dataclass(frozen=True)
class BenchmarkRecord:
    seed: int
    scenario: str
    unitigs: int
    total_bases: int
    largest_contig: int
    n50: int
    reference_bases_recovered: int
    reference_fraction: float
    false_contigs: int
    rare_strain_bases_recovered: int
    damage_events: int
    dead_ends: int
    read_boundaries: int
    coverage_gaps: int
    filtered_unique: int
    filtered_conflicts: int
    low_quality_dead_ends: int
    terminal_only_dead_ends: int
    correction_candidates: int
    damage_protected_corrections: int
    consensus_corrected_bases: int
    consensus_true_corrections: int
    consensus_false_corrections: int
    consensus_precision: float
    consensus_unitigs: int
    consensus_largest_contig: int
    consensus_n50: int
    consensus_reference_bases_recovered: int
    consensus_reference_fraction: float
    consensus_false_contigs: int
    consensus_rare_strain_bases_recovered: int
    consensus_n50_delta: int
    consensus_reference_bases_delta: int
    consensus_false_contigs_delta: int
    oracle_unique_rescued_kmers: int
    oracle_unique_unitigs: int
    oracle_unique_total_bases: int
    oracle_unique_largest_contig: int
    oracle_unique_n50: int
    oracle_unique_reference_bases_recovered: int
    oracle_unique_reference_fraction: float
    oracle_unique_false_contigs: int
    oracle_unique_rare_strain_bases_recovered: int
    oracle_unique_n50_delta: int
    oracle_unique_largest_contig_delta: int
    oracle_unique_reference_bases_delta: int
    oracle_unique_false_contigs_delta: int
    oracle_all_truth_rescued_kmers: int
    oracle_all_truth_unitigs: int
    oracle_all_truth_total_bases: int
    oracle_all_truth_largest_contig: int
    oracle_all_truth_n50: int
    oracle_all_truth_reference_bases_recovered: int
    oracle_all_truth_reference_fraction: float
    oracle_all_truth_false_contigs: int
    oracle_all_truth_rare_strain_bases_recovered: int
    oracle_all_truth_n50_delta: int
    oracle_all_truth_largest_contig_delta: int
    oracle_all_truth_reference_bases_delta: int
    oracle_all_truth_false_contigs_delta: int


def _load_config(path: Path) -> dict[str, object]:
    config = json.loads(path.read_text(encoding="utf-8"))
    if config.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("unsupported benchmark ladder schema")
    required = {
        "reference_length",
        "read_length",
        "coverage",
        "k",
        "min_count",
        "seeds",
        "damage_profile",
        "sequencing_error_rate",
        "rare_strain_coverage",
        "contamination_coverage",
        "scenarios",
        "comparators",
        "metrics",
    }
    missing = sorted(required - set(config))
    if missing:
        raise ValueError(f"benchmark config missing fields: {', '.join(missing)}")
    scenarios = config["scenarios"]
    if not isinstance(scenarios, list) or not scenarios:
        raise ValueError("benchmark scenarios must be a non-empty list")
    supported = {
        "clean",
        "damage",
        "sequencing_error",
        "damage_error",
        "rare_strain",
        "contamination",
    }
    unknown = sorted(set(scenarios) - supported)
    if unknown:
        raise ValueError(f"unsupported benchmark scenarios: {', '.join(unknown)}")
    seeds = config["seeds"]
    if not isinstance(seeds, list) or not seeds or any(
        not isinstance(seed, int) or isinstance(seed, bool) for seed in seeds
    ):
        raise ValueError("benchmark seeds must be a non-empty list of integers")
    return config


def _random_reference(length: int, seed: int) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def _mutate_reference(reference: str, seed: int, count: int = 50) -> str:
    generator = random.Random(seed)
    sequence = list(reference)
    positions = generator.sample(range(len(sequence)), min(count, len(sequence)))
    for position in positions:
        sequence[position] = generator.choice(
            [base for base in "ACGT" if base != sequence[position]]
        )
    return "".join(sequence)


def _n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    total = sum(lengths)
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running * 2 >= total:
            return length
    return 0


def _covered_bases(contigs: list[str], references: tuple[str, ...]) -> int:
    intervals: list[tuple[int, int]] = []
    for contig in contigs:
        for reference in references:
            position = reference.find(contig)
            if position >= 0:
                intervals.append((position, position + len(contig)))
                break
            position = reverse_complement(reference).find(contig)
            if position >= 0:
                intervals.append((position, position + len(contig)))
                break
    covered = 0
    end = -1
    for start, stop in sorted(intervals):
        if start > end:
            covered += stop - start
        elif stop > end:
            covered += stop - end
        end = max(end, stop)
    return covered


def _truth_kmers(references: tuple[str, ...], k: int) -> set[str]:
    return {
        canonical_sequence(reference[start : start + k])
        for reference in references
        for start in range(len(reference) - k + 1)
    }


def _oracle_supplements(
    reads: list[SimulatedRead],
    references: tuple[str, ...],
    k: int,
    min_count: int,
) -> tuple[list[str], int]:
    """Add only enough support to retain observed, truth-valid weak k-mers."""
    counts: dict[str, int] = {}
    for read in reads:
        for start in range(len(read.sequence) - k + 1):
            kmer = canonical_sequence(read.sequence[start : start + k])
            counts[kmer] = counts.get(kmer, 0) + 1

    truth = _truth_kmers(references, k)
    rescued = {
        kmer: count
        for kmer, count in counts.items()
        if count < min_count and kmer in truth
    }
    supplements = [
        kmer
        for kmer, count in sorted(rescued.items())
        for _ in range(min_count - count)
    ]
    return supplements, len(rescued)


def _oracle_unique_kmers(
    graph: CompactedUnitigGraph,
    reads: list[SimulatedRead],
    references: tuple[str, ...],
    min_count: int,
) -> set[str]:
    """Return truth-valid weak k-mers that uniquely continue a dead end."""
    context_length = graph.k - 1
    contexts = {
        graph.oriented_sequence(oriented_handle)[-context_length:]
        for unitig_id in range(graph.unitig_count)
        for oriented_handle in ((unitig_id << 1) | 1, unitig_id << 1)
        if not graph.outgoing(oriented_handle)
    }
    extensions = {context: [0, 0, 0, 0] for context in contexts}
    base_index = {base: index for index, base in enumerate("ACGT")}
    for read in reads:
        for sequence in (read.sequence, reverse_complement(read.sequence)):
            for start in range(len(sequence) - graph.k + 1):
                context = sequence[start : start + context_length]
                counts = extensions.get(context)
                if counts is None:
                    continue
                base = sequence[start + context_length]
                if base in base_index:
                    counts[base_index[base]] += 1

    truth = _truth_kmers(references, graph.k)
    global_counts: dict[str, int] = {}
    for read in reads:
        for start in range(len(read.sequence) - graph.k + 1):
            kmer = canonical_sequence(read.sequence[start : start + graph.k])
            global_counts[kmer] = global_counts.get(kmer, 0) + 1
    candidates: set[str] = set()
    for context, counts in extensions.items():
        observed = [index for index, count in enumerate(counts) if count]
        if len(observed) != 1 or counts[observed[0]] >= min_count:
            continue
        kmer = canonical_sequence(context + "ACGT"[observed[0]])
        if kmer in truth and global_counts.get(kmer, 0) < min_count:
            candidates.add(kmer)
    return candidates


def _restricted_oracle_supplements(
    reads: list[SimulatedRead],
    candidates: set[str],
    k: int,
    min_count: int,
) -> list[str]:
    counts: dict[str, int] = {}
    for read in reads:
        for start in range(len(read.sequence) - k + 1):
            kmer = canonical_sequence(read.sequence[start : start + k])
            counts[kmer] = counts.get(kmer, 0) + 1
    return [
        kmer
        for kmer in sorted(candidates)
        for _ in range(max(0, min_count - counts.get(kmer, 0)))
    ]


def _scenario_reads(
    reference: str,
    seed: int,
    config: dict[str, object],
    scenario: str,
) -> tuple[list[SimulatedRead], tuple[str, ...], tuple[str, ...], int]:
    coverage = float(config["coverage"])
    reads = simulate_fragments(
        reference, int(config["read_length"]), coverage, seed + 100
    )
    truth = (reference,)
    rare_truth: tuple[str, ...] = ()
    damage_events = 0
    profile = config["damage_profile"]
    if scenario in {"damage", "damage_error", "rare_strain", "contamination"}:
        reads = add_terminal_damage(reads, profile, profile, seed + 200)
        damage_events += sum(len(read.damage_positions) for read in reads)
    if scenario in {"sequencing_error", "damage_error"}:
        reads = add_sequencing_errors(
            reads, float(config["sequencing_error_rate"]), seed + 300
        )
    if scenario == "rare_strain":
        rare = _mutate_reference(reference, seed + 400)
        rare_reads = simulate_fragments(
            rare,
            int(config["read_length"]),
            float(config["rare_strain_coverage"]),
            seed + 500,
        )
        reads = reads + add_terminal_damage(rare_reads, profile, profile, seed + 600)
        truth = (reference, rare)
        rare_truth = (rare,)
    elif scenario == "contamination":
        contaminant = _random_reference(len(reference), seed + 700)
        contaminant_reads = simulate_fragments(
            contaminant,
            int(config["read_length"]),
            float(config["contamination_coverage"]),
            seed + 800,
        )
        reads = reads + contaminant_reads
        truth = (reference, contaminant)
    return reads, truth, rare_truth, damage_events


def _run_record(
    reference: str,
    seed: int,
    config: dict[str, object],
    scenario: str,
) -> BenchmarkRecord:
    reads, truth, rare_truth, damage_events = _scenario_reads(
        reference, seed, config, scenario
    )
    graph = build_bidirected_dbg(
        [read.sequence for read in reads],
        int(config["k"]),
        int(config["min_count"]),
    )
    unitig_graph = build_compacted_unitig_graph(graph)
    contigs = list(unitig_graph.sequences)
    audited_reads = [
        Read(
            read.name,
            read.sequence,
            tuple(
                8 if position in read.error_positions else 35
                for position in range(len(read.sequence))
            ),
        )
        for read in reads
    ]
    with tempfile.TemporaryDirectory() as directory:
        audit_root = Path(directory)
        dead_end_summary = write_dead_end_attribution_report(
            unitig_graph,
            audited_reads,
            audit_root / "dead-ends.tsv",
            min_count=int(config["min_count"]),
        )
        correction_summary = write_correction_report(
            audited_reads,
            audit_root / "corrections.tsv",
            k=int(config["k"]),
            min_count=int(config["min_count"]),
        )
    lengths = [len(contig) for contig in contigs]
    exact_contigs = [
        contig
        for contig in contigs
        if any(
            contig in reference or contig in reverse_complement(reference)
            for reference in truth
        )
    ]
    reference_bases = _covered_bases(exact_contigs, (truth[0],))
    rare_bases = _covered_bases(exact_contigs, rare_truth) if rare_truth else 0
    consensus_reads, consensus_summary = apply_damage_aware_consensus(
        audited_reads,
        end_window=len(config["damage_profile"]),
    )
    true_corrections = 0
    false_corrections = 0
    for simulated, corrected in zip(reads, consensus_reads, strict=True):
        for position, (before, after) in enumerate(
            zip(simulated.sequence, corrected.sequence, strict=True)
        ):
            if before == after:
                continue
            if after == simulated.source_sequence[position]:
                true_corrections += 1
            else:
                false_corrections += 1
    consensus_graph = build_bidirected_dbg(
        [read.sequence for read in consensus_reads],
        int(config["k"]),
        int(config["min_count"]),
    )
    consensus_contigs = list(
        build_compacted_unitig_graph(consensus_graph).sequences
    )
    consensus_lengths = [len(contig) for contig in consensus_contigs]
    consensus_exact_contigs = [
        contig
        for contig in consensus_contigs
        if any(
            contig in truth_reference
            or contig in reverse_complement(truth_reference)
            for truth_reference in truth
        )
    ]
    consensus_reference_bases = _covered_bases(
        consensus_exact_contigs, (truth[0],)
    )
    consensus_rare_bases = (
        _covered_bases(consensus_exact_contigs, rare_truth)
        if rare_truth
        else 0
    )
    consensus_false_contigs = (
        len(consensus_contigs) - len(consensus_exact_contigs)
    )
    unique_kmers = _oracle_unique_kmers(
        unitig_graph,
        reads,
        truth,
        int(config["min_count"]),
    )
    unique_supplements = _restricted_oracle_supplements(
        reads,
        unique_kmers,
        int(config["k"]),
        int(config["min_count"]),
    )
    unique_graph = build_bidirected_dbg(
        [read.sequence for read in reads] + unique_supplements,
        int(config["k"]),
        int(config["min_count"]),
    )
    unique_contigs = list(build_compacted_unitig_graph(unique_graph).sequences)
    unique_lengths = [len(contig) for contig in unique_contigs]
    unique_exact_contigs = [
        contig
        for contig in unique_contigs
        if any(
            contig in truth_reference
            or contig in reverse_complement(truth_reference)
            for truth_reference in truth
        )
    ]
    unique_reference_bases = _covered_bases(
        unique_exact_contigs, (truth[0],)
    )
    unique_rare_bases = (
        _covered_bases(unique_exact_contigs, rare_truth) if rare_truth else 0
    )
    unique_false_contigs = len(unique_contigs) - len(unique_exact_contigs)

    all_truth_supplements, all_truth_rescued_kmers = _oracle_supplements(
        reads,
        truth,
        int(config["k"]),
        int(config["min_count"]),
    )
    all_truth_graph = build_bidirected_dbg(
        [read.sequence for read in reads] + all_truth_supplements,
        int(config["k"]),
        int(config["min_count"]),
    )
    all_truth_contigs = list(
        build_compacted_unitig_graph(all_truth_graph).sequences
    )
    all_truth_lengths = [len(contig) for contig in all_truth_contigs]
    all_truth_exact_contigs = [
        contig
        for contig in all_truth_contigs
        if any(
            contig in truth_reference
            or contig in reverse_complement(truth_reference)
            for truth_reference in truth
        )
    ]
    all_truth_reference_bases = _covered_bases(
        all_truth_exact_contigs, (truth[0],)
    )
    all_truth_rare_bases = (
        _covered_bases(all_truth_exact_contigs, rare_truth)
        if rare_truth
        else 0
    )
    all_truth_false_contigs = (
        len(all_truth_contigs) - len(all_truth_exact_contigs)
    )
    return BenchmarkRecord(
        seed=seed,
        scenario=scenario,
        unitigs=len(contigs),
        total_bases=sum(lengths),
        largest_contig=max(lengths, default=0),
        n50=_n50(lengths),
        reference_bases_recovered=reference_bases,
        reference_fraction=reference_bases / len(truth[0]),
        false_contigs=len(contigs) - len(exact_contigs),
        rare_strain_bases_recovered=rare_bases,
        damage_events=damage_events,
        dead_ends=dead_end_summary.dead_ends,
        read_boundaries=dead_end_summary.read_boundaries,
        coverage_gaps=dead_end_summary.coverage_gaps,
        filtered_unique=dead_end_summary.filtered_unique,
        filtered_conflicts=dead_end_summary.filtered_conflicts,
        low_quality_dead_ends=dead_end_summary.low_quality,
        terminal_only_dead_ends=dead_end_summary.terminal_only,
        correction_candidates=correction_summary.would_correct,
        damage_protected_corrections=correction_summary.protected_damage,
        consensus_corrected_bases=consensus_summary.corrected_bases,
        consensus_true_corrections=true_corrections,
        consensus_false_corrections=false_corrections,
        consensus_precision=(
            true_corrections / (true_corrections + false_corrections)
            if true_corrections + false_corrections
            else 1.0
        ),
        consensus_unitigs=len(consensus_contigs),
        consensus_largest_contig=max(consensus_lengths, default=0),
        consensus_n50=_n50(consensus_lengths),
        consensus_reference_bases_recovered=consensus_reference_bases,
        consensus_reference_fraction=consensus_reference_bases / len(truth[0]),
        consensus_false_contigs=consensus_false_contigs,
        consensus_rare_strain_bases_recovered=consensus_rare_bases,
        consensus_n50_delta=_n50(consensus_lengths) - _n50(lengths),
        consensus_reference_bases_delta=(
            consensus_reference_bases - reference_bases
        ),
        consensus_false_contigs_delta=(
            consensus_false_contigs - (len(contigs) - len(exact_contigs))
        ),
        oracle_unique_rescued_kmers=len(unique_kmers),
        oracle_unique_unitigs=len(unique_contigs),
        oracle_unique_total_bases=sum(unique_lengths),
        oracle_unique_largest_contig=max(unique_lengths, default=0),
        oracle_unique_n50=_n50(unique_lengths),
        oracle_unique_reference_bases_recovered=unique_reference_bases,
        oracle_unique_reference_fraction=unique_reference_bases / len(truth[0]),
        oracle_unique_false_contigs=unique_false_contigs,
        oracle_unique_rare_strain_bases_recovered=unique_rare_bases,
        oracle_unique_n50_delta=_n50(unique_lengths) - _n50(lengths),
        oracle_unique_largest_contig_delta=(
            max(unique_lengths, default=0) - max(lengths, default=0)
        ),
        oracle_unique_reference_bases_delta=(
            unique_reference_bases - reference_bases
        ),
        oracle_unique_false_contigs_delta=(
            unique_false_contigs - (len(contigs) - len(exact_contigs))
        ),
        oracle_all_truth_rescued_kmers=all_truth_rescued_kmers,
        oracle_all_truth_unitigs=len(all_truth_contigs),
        oracle_all_truth_total_bases=sum(all_truth_lengths),
        oracle_all_truth_largest_contig=max(all_truth_lengths, default=0),
        oracle_all_truth_n50=_n50(all_truth_lengths),
        oracle_all_truth_reference_bases_recovered=all_truth_reference_bases,
        oracle_all_truth_reference_fraction=(
            all_truth_reference_bases / len(truth[0])
        ),
        oracle_all_truth_false_contigs=all_truth_false_contigs,
        oracle_all_truth_rare_strain_bases_recovered=all_truth_rare_bases,
        oracle_all_truth_n50_delta=_n50(all_truth_lengths) - _n50(lengths),
        oracle_all_truth_largest_contig_delta=(
            max(all_truth_lengths, default=0) - max(lengths, default=0)
        ),
        oracle_all_truth_reference_bases_delta=(
            all_truth_reference_bases - reference_bases
        ),
        oracle_all_truth_false_contigs_delta=(
            all_truth_false_contigs - (len(contigs) - len(exact_contigs))
        ),
    )


def _git_state(project_root: Path) -> dict[str, object]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=project_root,
        check=False,
        capture_output=True,
        text=True,
    )
    status = subprocess.run(
        ["git", "status", "--short"],
        cwd=project_root,
        check=False,
        capture_output=True,
        text=True,
    )
    return {
        "commit": commit.stdout.strip() if commit.returncode == 0 else None,
        "dirty": bool(status.stdout.strip()) if status.returncode == 0 else None,
    }


def _write_outputs(output_dir: Path, config_path: Path, config: dict[str, object], records: list[BenchmarkRecord]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(records[0])))
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)
    (output_dir / "summary.json").write_text(
        json.dumps([asdict(record) for record in records], indent=2) + "\n",
        encoding="utf-8",
    )
    project_root = Path(__file__).resolve().parents[1]
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project": _git_state(project_root),
        "runtime": {"python": sys.version, "platform": platform.platform()},
        "config": {
            "path": str(config_path.resolve()),
            "sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
            "values": config,
        },
        "outputs": {
            name: {
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for name in ("summary.tsv", "summary.json")
            for path in [output_dir / name]
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def run_benchmark(config_path: Path, output_dir: Path) -> list[BenchmarkRecord]:
    config = _load_config(config_path)
    records: list[BenchmarkRecord] = []
    for seed in config["seeds"]:
        reference = _random_reference(int(config["reference_length"]), seed)
        for scenario in config["scenarios"]:
            records.append(_run_record(reference, seed, config, scenario))
    _write_outputs(output_dir, config_path, config, records)
    return records


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main() -> int:
    arguments = build_parser().parse_args()
    config = _load_config(arguments.config)
    run_count = len(config["seeds"]) * len(config["scenarios"])
    if arguments.dry_run:
        print(f"benchmark_schema={SCHEMA_VERSION}")
        print(f"runs={run_count}")
        print(f"k={config['k']}")
        print(f"min_count={config['min_count']}")
        print(f"comparators={','.join(config['comparators'])}")
        return 0
    records = run_benchmark(arguments.config, arguments.output_dir)
    print(f"runs={len(records)}")
    print(f"output={arguments.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
