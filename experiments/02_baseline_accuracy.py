"""Run the initial reference-based assembly accuracy benchmark."""

import argparse
import shlex
import shutil
import subprocess
from pathlib import Path


TOOL_NAMES = {
    "art": "art_illumina",
    "anvaya": "anvaya",
    "megahit": "megahit",
    "spades": "spades.py",
    "quast": "quast.py",
}


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

    print(f"[baseline] Results: {arguments.output_dir / 'quast' / 'report.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
