"""Validate conservative tip cleaning across clean simulated datasets."""

import argparse
import csv
import re
import shlex
import shutil
import subprocess
from pathlib import Path


TOOL_NAMES = {
    "art": "art_illumina",
    "anvaya": "anvaya",
    "quast": "quast.py",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare Anvaya with and without tip cleaning"
    )
    parser.add_argument("--reference", action="append", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seeds", nargs="+", type=int, default=[42, 43])
    parser.add_argument("--coverages", nargs="+", type=float, default=[5.0, 20.0])
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--resume", action="store_true", help="skip stages with verified outputs"
    )
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


def _condition_commands(
    reference: Path,
    condition_dir: Path,
    tools: dict[str, str],
    seed: int,
    coverage: float,
    k: int,
    threads: int,
) -> list[tuple[str, list[str]]]:
    reads_prefix = condition_dir / "reads" / "simulated"
    left = Path(f"{reads_prefix}1.fq")
    right = Path(f"{reads_prefix}2.fq")
    baseline = condition_dir / "baseline" / "contigs.fasta"
    cleaned = condition_dir / "cleaned" / "contigs.fasta"

    common_assembly = [
        tools["anvaya"],
        "assemble",
        "-1",
        str(left),
        "-2",
        str(right),
        "--k",
        str(k),
        "--min-count",
        "2",
        "--orientation-aware",
    ]
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
                "150",
                "-f",
                str(coverage),
                "-m",
                "300",
                "-s",
                "30",
                "-rs",
                str(seed),
                "-na",
                "-o",
                str(reads_prefix),
            ],
        ),
        ("baseline", [*common_assembly, "-o", str(baseline)]),
        (
            "cleaned",
            [*common_assembly, "--clean-tips", "-o", str(cleaned)],
        ),
        (
            "quast",
            [
                tools["quast"],
                str(baseline),
                str(cleaned),
                "-r",
                str(reference),
                "--labels",
                "Baseline,Tip_cleaning",
                "--min-contig",
                "200",
                "-t",
                str(threads),
                "-o",
                str(condition_dir / "quast"),
            ],
        ),
    ]


def _conditions(arguments: argparse.Namespace, tools: dict[str, str]):
    for reference in arguments.reference:
        reference_name = reference.stem
        for coverage in arguments.coverages:
            coverage_name = f"cov_{coverage:g}x"
            for seed in arguments.seeds:
                condition_dir = (
                    arguments.output_dir
                    / reference_name
                    / coverage_name
                    / f"seed_{seed}"
                )
                yield condition_dir, _condition_commands(
                    reference.resolve(),
                    condition_dir.resolve(),
                    tools,
                    seed,
                    coverage,
                    arguments.k,
                    arguments.threads,
                )


def _completed_stage(condition_dir: Path, name: str) -> bool:
    timing = condition_dir / "timing" / f"{name}.txt"
    if not timing.is_file() or "Exit status: 0" not in timing.read_text():
        return False
    if name == "art":
        outputs = [
            condition_dir / "reads" / "simulated1.fq",
            condition_dir / "reads" / "simulated2.fq",
        ]
    elif name in {"baseline", "cleaned"}:
        outputs = [condition_dir / name / "contigs.fasta"]
    else:
        outputs = [condition_dir / "quast" / "report.tsv"]
    return all(path.is_file() and path.stat().st_size > 0 for path in outputs)


def _timing_values(path: Path) -> tuple[str, str]:
    text = path.read_text()
    wall = re.search(r"Elapsed \(wall clock\) time .*: (.+)", text)
    memory = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if wall is None or memory is None:
        raise ValueError(f"incomplete timing report: {path}")
    return wall.group(1), memory.group(1)


def _write_summary(output_dir: Path) -> None:
    columns = [
        "reference",
        "coverage",
        "seed",
        "tips_removed",
        "baseline_n50",
        "cleaned_n50",
        "baseline_largest_contig",
        "cleaned_largest_contig",
        "baseline_genome_fraction",
        "cleaned_genome_fraction",
        "baseline_misassemblies",
        "cleaned_misassemblies",
        "baseline_duplication_ratio",
        "cleaned_duplication_ratio",
        "baseline_wall_time",
        "cleaned_wall_time",
        "baseline_peak_memory_kib",
        "cleaned_peak_memory_kib",
    ]
    with (output_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        for report in sorted(output_dir.glob("*/*/seed_*/quast/report.tsv")):
            condition_dir = report.parents[1]
            rows = {
                row[0]: row[1:]
                for row in csv.reader(report.open(encoding="utf-8"), delimiter="\t")
            }
            tip_log = (condition_dir / "cleaned.log").read_text()
            tips = re.search(r"tips_removed=(\d+)", tip_log)
            if tips is None:
                raise ValueError(f"tip count missing from {condition_dir / 'cleaned.log'}")
            baseline_time = _timing_values(condition_dir / "timing" / "baseline.txt")
            cleaned_time = _timing_values(condition_dir / "timing" / "cleaned.txt")
            relative = condition_dir.relative_to(output_dir)
            writer.writerow(
                [
                    relative.parts[0],
                    relative.parts[1],
                    relative.parts[2].removeprefix("seed_"),
                    tips.group(1),
                    *rows["N50"],
                    *rows["Largest contig"],
                    *rows["Genome fraction (%)"],
                    *rows["# misassemblies"],
                    *rows["Duplication ratio"],
                    baseline_time[0],
                    cleaned_time[0],
                    baseline_time[1],
                    cleaned_time[1],
                ]
            )
def main() -> int:
    arguments = build_parser().parse_args()
    if arguments.k < 2 or arguments.threads < 1:
        raise SystemExit("--k must be at least 2 and --threads must be positive")
    if any(coverage <= 0 for coverage in arguments.coverages):
        raise SystemExit("coverage values must be positive")
    missing_references = [path for path in arguments.reference if not path.is_file()]
    if missing_references:
        raise SystemExit(f"reference not found: {missing_references[0]}")

    tools = TOOL_NAMES if arguments.dry_run else _find_tools()
    conditions = list(_conditions(arguments, tools))
    if arguments.dry_run:
        for condition_dir, commands in conditions:
            print(f"[condition] {condition_dir}")
            for name, command in commands:
                print(f"[{name}] {shlex.join(command)}")
        return 0

    if arguments.output_dir.exists() and not arguments.resume:
        raise SystemExit(f"output directory already exists: {arguments.output_dir}")

    for condition_dir, commands in conditions:
        (condition_dir / "reads").mkdir(parents=True, exist_ok=True)
        (condition_dir / "baseline").mkdir(exist_ok=True)
        (condition_dir / "cleaned").mkdir(exist_ok=True)
        (condition_dir / "timing").mkdir(exist_ok=True)
        with (condition_dir / "commands.txt").open(
            "a" if arguments.resume else "w", encoding="utf-8"
        ) as commands_handle:
            for name, command in commands:
                if arguments.resume and _completed_stage(condition_dir, name):
                    print(f"[validation] {condition_dir.name}: skipping {name}")
                    continue
                print(f"[validation] {condition_dir.name}: {name}", flush=True)
                commands_handle.write(shlex.join(command) + "\n")
                commands_handle.flush()
                with (condition_dir / f"{name}.log").open(
                    "w", encoding="utf-8"
                ) as log_handle:
                    subprocess.run(
                        [
                            "/usr/bin/time",
                            "-v",
                            "-o",
                            str(condition_dir / "timing" / f"{name}.txt"),
                            *command,
                        ],
                        check=True,
                        stdout=log_handle,
                        stderr=subprocess.STDOUT,
                    )

    _write_summary(arguments.output_dir)
    print(f"[validation] Results: {arguments.output_dir / 'summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
