"""Run the frozen, truth-labelled Anvaya validation ladder."""

import argparse
import csv
import hashlib
import json
import platform
import random
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.sequences import reverse_complement
from anvaya.unitig_graph import build_compacted_unitig_graph
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
