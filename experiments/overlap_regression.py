"""Run deterministic overlap CLI fixtures; this is not an accuracy benchmark."""

import argparse
import contextlib
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))

from anvaya.cli import main as assemble_cli
from anvaya.reads import load_reads
from helpers.simulation import (
    add_sequencing_errors,
    add_terminal_damage,
    simulate_fragments,
)


def _reference(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _reads(config: dict, seed: int, scenario: str):
    reference = _reference(config["reference_length"], seed)
    reads = simulate_fragments(
        reference, config["read_length"], config["coverage"], seed + 100,
        reverse_fraction=0.5,
    )
    profile = config["damage_profile"]
    if scenario in {"damage", "damage_error", "rare_strain", "contamination"}:
        reads = add_terminal_damage(reads, profile, profile, seed + 200)
    if scenario in {"sequencing_error", "damage_error"}:
        reads = add_sequencing_errors(reads, config["sequencing_error_rate"], seed + 300)
    if scenario in {"rare_strain", "contamination"}:
        other = _reference(len(reference), seed + 400)
        coverage = config["contamination_coverage"]
        if scenario == "rare_strain":
            rng = random.Random(seed + 400)
            bases = list(reference)
            for position in rng.sample(range(len(bases)), min(50, len(bases))):
                bases[position] = rng.choice(
                    [base for base in "ACGT" if base != bases[position]]
                )
            other = "".join(bases)
            coverage = config["rare_strain_coverage"]
        extra = simulate_fragments(
            other, config["read_length"], coverage, seed + 500,
            reverse_fraction=0.5,
        )
        if scenario == "rare_strain":
            extra = add_terminal_damage(extra, profile, profile, seed + 600)
        reads += extra
    return reads


def run_regression(config_path: Path, output_dir: Path) -> dict:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    scenarios = {
        "clean", "damage", "sequencing_error", "damage_error",
        "rare_strain", "contamination",
    }
    if (
        not config["seeds"] or not config["scenarios"]
        or not set(config["scenarios"]) <= scenarios
    ):
        raise ValueError("regression requires seeds and supported scenarios")
    output_dir.mkdir(parents=True, exist_ok=False)
    records = []
    for seed in config["seeds"]:
        for scenario in config["scenarios"]:
            run = output_dir / f"{seed}-{scenario}"
            run.mkdir()
            reads = _reads(config, seed, scenario)
            input_path = run / "reads.fastq"
            with input_path.open("w", encoding="utf-8", newline="\n") as handle:
                for index, read in enumerate(reads):
                    # Uniform Q30 does not encode the simulated error labels.
                    handle.write(
                        f"@fragment_{index}\n{read.sequence}\n+\n"
                        f"{'?' * len(read.sequence)}\n"
                    )
            args = [
                "overlap-assemble", "-i", str(input_path),
                *config["assembly_args"],
                "--progressive-raw-phase-projection", str(run / "progressive.fasta"),
                "--selective-support-rescue-projection", str(run / "selective.fasta"),
                "--high-confidence-support-two-rescue-projection",
                str(run / "support-two.fasta"),
                "-o", str(run / "contigs.fasta"),
            ]
            with (run / "assembly.log").open("w", encoding="utf-8") as log:
                with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                    if assemble_cli(args) != 0:
                        raise RuntimeError("overlap assembly failed")
            outputs = {}
            for name in (
                "contigs.fasta", "progressive.fasta",
                "selective.fasta", "support-two.fasta",
            ):
                path = run / name
                content = path.read_bytes()
                contigs = load_reads(path) if content else []
                outputs[name] = {
                    "sha256": hashlib.sha256(content).hexdigest(),
                    "contigs": len(contigs),
                    "bases": sum(len(read.sequence) for read in contigs),
                }
            records.append({
                "seed": seed,
                "scenario": scenario,
                "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
                "outputs": outputs,
            })
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT,
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    dirty = bool(subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=normal"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout.strip())
    result = {
        "schema_version": 1,
        "purpose": "structural_regression_not_biological_acceptance",
        "commit": commit,
        "dirty": dirty,
        "python": sys.version,
        "config": config,
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "records": records,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8",
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=Path(__file__).with_name("overlap_regression.json"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--compare", type=Path,
        help="baseline manifest whose input/output hashes must match",
    )
    args = parser.parse_args()
    result = run_regression(args.config, args.output_dir)
    if args.compare:
        baseline = json.loads(args.compare.read_text(encoding="utf-8"))
        if (
            result["config"] != baseline["config"]
            or result["records"] != baseline["records"]
        ):
            parser.exit(
                1, "Regression mismatch: configuration, inputs or FASTA outputs changed\n"
            )
        print("All input and FASTA checksums match the baseline")
    print(f"runs={len(result['records'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
