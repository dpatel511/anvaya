"""Reproducible safe-mode execution of the external CarpeDeam assembler."""

import gzip
import hashlib
import json
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.raw_consensus import DamageProfile
from anvaya.sequences import normalize_dna


@dataclass(frozen=True, slots=True)
class InputScan:
    """Streaming input summary; no reads are retained in memory."""

    format: str
    reads: int
    bases: int
    minimum_length: int
    maximum_length: int


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _open_text(path: Path):
    return (
        gzip.open(path, mode="rt", encoding="utf-8")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8")
    )


def scan_sequence_input(path: Path) -> InputScan:
    """Validate FASTA/FASTQ records while retaining only length statistics."""
    reads = bases = 0
    minimum = None
    maximum = 0

    def observe(sequence: str) -> None:
        nonlocal reads, bases, minimum, maximum
        length = len(normalize_dna(sequence))
        if length == 0:
            raise ValueError("sequence records must not be empty")
        reads += 1
        bases += length
        minimum = length if minimum is None else min(minimum, length)
        maximum = max(maximum, length)

    with _open_text(path) as handle:
        first = handle.readline().rstrip("\r\n")
        if not first:
            raise ValueError("sequence input must not be empty")
        if first.startswith("@"):
            header = first
            while header:
                if not header.startswith("@") or not header[1:].strip():
                    raise ValueError("FASTQ records require non-empty '@' headers")
                sequence = handle.readline().strip()
                separator = handle.readline().rstrip("\r\n")
                qualities = handle.readline().rstrip("\r\n")
                if not sequence or not separator.startswith("+") or not qualities:
                    raise ValueError("incomplete FASTQ record")
                if len(sequence) != len(qualities):
                    raise ValueError("FASTQ sequence and quality lengths must match")
                observe(sequence)
                header = handle.readline().rstrip("\r\n")
                while header == "":
                    next_line = handle.readline()
                    if not next_line:
                        break
                    header = next_line.rstrip("\r\n")
            format_name = "FASTQ"
        elif first.startswith(">"):
            if not first[1:].strip():
                raise ValueError("FASTA records require non-empty '>' headers")
            sequence_parts = []
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if not line[1:].strip():
                        raise ValueError("FASTA records require non-empty '>' headers")
                    observe("".join(sequence_parts))
                    sequence_parts = []
                else:
                    sequence_parts.append(line)
            observe("".join(sequence_parts))
            format_name = "FASTA"
        else:
            raise ValueError("input must be FASTA or FASTQ")

    return InputScan(format_name, reads, bases, minimum or 0, maximum)


def _fasta_lengths(path: Path) -> list[int]:
    lengths = []
    current = 0
    saw_header = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not line[1:].strip():
                    raise ValueError("CarpeDeam emitted a FASTA record without a name")
                if current:
                    lengths.append(current)
                elif saw_header:
                    raise ValueError("CarpeDeam emitted an empty FASTA record")
                saw_header = True
                current = 0
            else:
                if not saw_header:
                    raise ValueError("CarpeDeam output is not FASTA")
                current += len(normalize_dna(line))
    if current:
        lengths.append(current)
    if not lengths:
        raise ValueError("CarpeDeam did not emit any FASTA contigs")
    return lengths


def _n50(lengths: list[int]) -> int:
    threshold = (sum(lengths) + 1) // 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            return length
    return 0


def _resolve_executable(value: str | Path) -> Path:
    candidate = Path(value)
    resolved = (
        candidate.resolve()
        if candidate.parent != Path(".") or candidate.exists()
        else Path(shutil.which(str(value)) or "")
    )
    if not str(resolved) or not resolved.is_file():
        raise ValueError(f"CarpeDeam executable not found: {value}")
    return resolved.resolve()


def run_carpedeam_safe(
    input_path: Path,
    profile_prefix: Path,
    output_path: Path,
    temporary_directory: Path,
    diagnostics_path: Path,
    *,
    executable: str | Path = "carpedeam",
    threads: int = 2,
    minimum_contig_length: int = 31,
) -> dict:
    """Run a pinned external binary in safe mode and write an audit manifest."""
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if minimum_contig_length < 1:
        raise ValueError("minimum contig length must be at least 1")
    input_path = input_path.resolve()
    output_path = output_path.resolve()
    temporary_directory = temporary_directory.resolve()
    diagnostics_path = diagnostics_path.resolve()
    profile_paths = [Path(f"{profile_prefix}{end}.prof").resolve() for end in ("5p", "3p")]
    if not input_path.is_file():
        raise ValueError(f"input does not exist: {input_path}")
    DamageProfile.from_prefix(profile_prefix)
    reserved = [input_path, output_path, temporary_directory, diagnostics_path, *profile_paths]
    if len(reserved) != len(set(reserved)):
        raise ValueError("input, profiles, output, temporary directory and diagnostics must be distinct")
    for path in (output_path, temporary_directory, diagnostics_path):
        if path.exists():
            raise ValueError(f"refusing to reuse existing path: {path}")

    scan = scan_sequence_input(input_path)
    if scan.minimum_length < 20:
        raise ValueError(
            f"CarpeDeam requires reads of at least 20 bp; observed {scan.minimum_length}"
        )
    executable_path = _resolve_executable(executable)
    help_result = subprocess.run(
        [str(executable_path), "ancient_assemble", "-h"],
        check=False,
        capture_output=True,
        text=True,
    )
    help_text = (help_result.stdout or "") + (help_result.stderr or "")
    required_help = (
        "usage: carpedeam ancient_assemble",
        "--ancient-damage",
        "--unsafe",
        "--min-contig-len",
        "--num-iter-reads-only",
        "--num-iterations",
        "--min-merge-seq-id",
        "--min-cov-safe",
    )
    if any(marker not in help_text for marker in required_help):
        raise ValueError("executable does not expose the expected CarpeDeam safe-mode CLI")

    stdout_path = diagnostics_path.with_name(diagnostics_path.stem + "-stdout.txt")
    stderr_path = diagnostics_path.with_name(diagnostics_path.stem + "-stderr.txt")
    for path in (stdout_path, stderr_path):
        if path.exists():
            raise ValueError(f"refusing to reuse existing log path: {path}")
    command = [
        str(executable_path),
        "ancient_assemble",
        str(input_path),
        str(output_path),
        str(temporary_directory),
        "--ancient-damage",
        str(profile_prefix.resolve()),
        "--unsafe",
        "0",
        "--min-contig-len",
        str(minimum_contig_length),
        "--threads",
        str(threads),
        "--num-iter-reads-only",
        "5",
        "--num-iterations",
        "10",
        "--min-merge-seq-id",
        "0.99",
        "--min-cov-safe",
        "5",
    ]
    started = time.perf_counter()
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        completed = subprocess.run(
            command,
            check=False,
            stdout=stdout,
            stderr=stderr,
            text=True,
        )
    elapsed = time.perf_counter() - started
    diagnostics = {
        "backend": "CarpeDeam",
        "mode": "safe",
        "status": "completed" if completed.returncode == 0 else "failed",
        "returncode": completed.returncode,
        "command": command,
        "elapsed_seconds": elapsed,
        "executable": str(executable_path),
        "executable_sha256": _sha256_file(executable_path),
        "help_sha256": hashlib.sha256(help_text.encode("utf-8")).hexdigest(),
        "input": str(input_path),
        "input_sha256": _sha256_file(input_path),
        "input_scan": asdict(scan),
        "profile_prefix": str(profile_prefix.resolve()),
        "profile_sha256": {
            path.name: _sha256_file(path) for path in profile_paths
        },
        "output": str(output_path),
        "temporary_directory": str(temporary_directory),
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
        "settings": {
            "unsafe": False,
            "minimum_contig_length": minimum_contig_length,
            "threads": threads,
            "raw_read_iterations": 5,
            "total_iterations": 10,
            "minimum_merge_identity": 0.99,
            "minimum_safe_coverage": 5,
        },
        "limitations": [
            "CarpeDeam output does not expose raw-read placements to Anvaya",
            "backend contiguity and sequence accuracy require separate evaluation",
            "executable SHA256 identifies the binary when a version flag is unavailable",
        ],
    }
    output_error = None
    if completed.returncode == 0:
        if not output_path.is_file():
            output_error = "CarpeDeam returned success without creating its output FASTA"
        else:
            try:
                lengths = _fasta_lengths(output_path)
            except (OSError, ValueError) as error:
                output_error = f"invalid CarpeDeam output FASTA: {error}"
            else:
                diagnostics["output_sha256"] = _sha256_file(output_path)
                diagnostics["output_contigs"] = len(lengths)
                diagnostics["output_bases"] = sum(lengths)
                diagnostics["output_n50"] = _n50(lengths)
                diagnostics["output_longest"] = max(lengths)
        if output_error is not None:
            diagnostics["status"] = "invalid_output"
            diagnostics["output_error"] = output_error
    diagnostics_path.write_text(
        json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8"
    )
    if completed.returncode != 0:
        raise ValueError(
            f"CarpeDeam failed with exit status {completed.returncode}; see {stderr_path}"
        )
    if output_error is not None:
        raise ValueError(output_error)
    return diagnostics
