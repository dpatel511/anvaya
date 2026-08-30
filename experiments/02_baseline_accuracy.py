"""Run the initial reference-based assembly accuracy benchmark."""

import argparse
import hashlib
import json
import platform
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


TOOL_NAMES = {
    "art": "art_illumina",
    "anvaya": "anvaya",
    "megahit": "megahit",
    "spades": "spades.py",
    "quast": "quast.py",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _write_manifest(
    path: Path,
    *,
    arguments: argparse.Namespace,
    tools: dict[str, str],
    commands: list[tuple[str, list[str]]],
    artifacts: list[Path],
) -> None:
    project_root = Path(__file__).resolve().parents[1]
    existing_artifacts = {
        str(artifact.resolve()): {
            "bytes": artifact.stat().st_size,
            "sha256": _sha256(artifact),
        }
        for artifact in artifacts
        if artifact.is_file()
    }
    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project": _git_state(project_root),
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
        },
        "parameters": {
            "reference": str(arguments.reference.resolve()),
            "k": arguments.k,
            "threads": arguments.threads,
            "coverage": arguments.coverage,
            "read_length": arguments.read_length,
            "fragment_mean": arguments.fragment_mean,
            "fragment_sd": arguments.fragment_sd,
            "seed": arguments.seed,
        },
        "tools": tools,
        "commands": [
            {"name": name, "argv": command} for name, command in commands
        ],
        "artifacts": existing_artifacts,
    }
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Simulate paired reads and compare Anvaya with established assemblers"
    )
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--coverage", type=float, default=20.0)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--fragment-mean", type=int, default=300)
    parser.add_argument("--fragment-sd", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--dry-run", action="store_true", help="print commands without running them"
    )
    return parser


def _find_tools() -> dict[str, str]:
    tools = {name: shutil.which(command) for name, command in TOOL_NAMES.items()}
    missing = [TOOL_NAMES[name] for name, path in tools.items() if path is None]
    if not Path("/usr/bin/time").is_file():
        missing.append("/usr/bin/time")
    if missing:
        raise SystemExit(f"missing required tools: {', '.join(missing)}")
    return {name: path for name, path in tools.items() if path is not None}


def _commands(
    reference: Path,
    output_dir: Path,
    tools: dict[str, str],
    k: int,
    threads: int,
    coverage: float,
    read_length: int,
    fragment_mean: int,
    fragment_sd: int,
    seed: int,
) -> list[tuple[str, list[str]]]:
    reads_prefix = output_dir / "reads" / "simulated"
    left = Path(f"{reads_prefix}1.fq")
    right = Path(f"{reads_prefix}2.fq")
    anvaya_contigs = output_dir / "anvaya" / "contigs.fasta"
    megahit_dir = output_dir / "megahit"
    spades_dir = output_dir / "metaspades"

    return [
        (
            "art",
            [
                tools["art"],
                "-ss",
                "HS25",
                "-i",
                str(reference),
                "-p",
                "-l",
                str(read_length),
                "-f",
                str(coverage),
                "-m",
                str(fragment_mean),
                "-s",
                str(fragment_sd),
                "-rs",
                str(seed),
                "-na",
                "-o",
                str(reads_prefix),
            ],
        ),
        (
            "anvaya",
            [
                tools["anvaya"],
                "assemble",
                "-1",
                str(left),
                "-2",
                str(right),
                "--k",
                str(k),
                "-o",
                str(anvaya_contigs),
            ],
        ),
        (
            "megahit",
            [
                tools["megahit"],
                "-1",
                str(left),
                "-2",
                str(right),
                "-t",
                str(threads),
                "-o",
                str(megahit_dir),
            ],
        ),
        (
            "metaspades",
            [
                tools["spades"],
                "--meta",
                "--only-assembler",
                "-1",
                str(left),
                "-2",
                str(right),
                "-t",
                str(threads),
                "-o",
                str(spades_dir),
            ],
        ),
        (
            "quast",
            [
                tools["quast"],
                str(anvaya_contigs),
                str(megahit_dir / "final.contigs.fa"),
                str(spades_dir / "contigs.fasta"),
                "-r",
                str(reference),
                "--labels",
                "Anvaya,MEGAHIT,metaSPAdes",
                "--min-contig",
                "200",
                "-t",
                str(threads),
                "-o",
                str(output_dir / "quast"),
            ],
        ),
    ]


def main() -> int:
    arguments = build_parser().parse_args()
    if arguments.k < 2 or arguments.threads < 1:
        raise SystemExit("--k must be at least 2 and --threads must be positive")
    if not arguments.reference.is_file():
        raise SystemExit(f"reference not found: {arguments.reference}")

    tools = TOOL_NAMES if arguments.dry_run else _find_tools()
    commands = _commands(
        arguments.reference.resolve(),
        arguments.output_dir.resolve(),
        tools,
        arguments.k,
        arguments.threads,
        arguments.coverage,
        arguments.read_length,
        arguments.fragment_mean,
        arguments.fragment_sd,
        arguments.seed,
    )

    if arguments.dry_run:
        for name, command in commands:
            print(f"[{name}] {shlex.join(command)}")
        return 0

    if arguments.output_dir.exists():
        raise SystemExit(f"output directory already exists: {arguments.output_dir}")
    (arguments.output_dir / "reads").mkdir(parents=True)
    (arguments.output_dir / "anvaya").mkdir()
    (arguments.output_dir / "timing").mkdir()

    command_log = arguments.output_dir / "commands.txt"
    with command_log.open("w", encoding="utf-8") as commands_handle:
        for name, command in commands:
            print(f"[baseline] Running {name}", flush=True)
            commands_handle.write(shlex.join(command) + "\n")
            commands_handle.flush()
            with (arguments.output_dir / f"{name}.log").open(
                "w", encoding="utf-8"
            ) as log_handle:
                subprocess.run(
                    [
                        "/usr/bin/time",
                        "-v",
                        "-o",
                        str(arguments.output_dir / "timing" / f"{name}.txt"),
                        *command,
                    ],
                    check=True,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                )

    reads_prefix = arguments.output_dir / "reads" / "simulated"
    _write_manifest(
        arguments.output_dir / "manifest.json",
        arguments=arguments,
        tools=tools,
        commands=commands,
        artifacts=[
            arguments.reference,
            Path(f"{reads_prefix}1.fq"),
            Path(f"{reads_prefix}2.fq"),
            arguments.output_dir / "anvaya" / "contigs.fasta",
            arguments.output_dir / "megahit" / "final.contigs.fa",
            arguments.output_dir / "metaspades" / "contigs.fasta",
            arguments.output_dir / "quast" / "report.txt",
        ],
    )

    print(f"[baseline] Results: {arguments.output_dir / 'quast' / 'report.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
