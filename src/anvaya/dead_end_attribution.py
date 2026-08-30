"""Report raw-read evidence explaining compacted-graph dead ends."""

import csv
from dataclasses import dataclass
from pathlib import Path

from anvaya.reads import Read
from anvaya.sequences import reverse_complement
from anvaya.unitig_graph import CompactedUnitigGraph


_BASE_INDEX = {base: index for index, base in enumerate("ACGT")}


@dataclass(slots=True, frozen=True)
class DeadEndAttributionSummary:
    """Counts of dead ends by their strongest observed cause."""

    dead_ends: int = 0
    read_boundaries: int = 0
    coverage_gaps: int = 0
    filtered_unique: int = 0
    filtered_conflicts: int = 0
    retained_context_elsewhere: int = 0
    low_quality: int = 0
    terminal_only: int = 0
    missing_quality: int = 0


@dataclass(slots=True)
class _Evidence:
    counts: list[int]
    terminal_counts: list[int]
    quality_sums: list[int]
    quality_counts: list[int]
    low_quality_counts: list[int]
    missing_quality_counts: list[int]
    boundary_stops: int = 0


def _cause(
    evidence: _Evidence,
    min_count: int,
) -> str:
    observed = [count for count in evidence.counts if count]
    if any(count >= min_count for count in observed):
        return "retained_context_elsewhere"
    if len(observed) > 1:
        return "filtered_conflict"
    if len(observed) == 1:
        index = next(index for index, count in enumerate(evidence.counts) if count)
        count = evidence.counts[index]
        if evidence.missing_quality_counts[index] == count:
            return "filtered_unique_missing_quality"
        if evidence.low_quality_counts[index] * 2 >= count:
            return "filtered_unique_low_quality"
        if evidence.terminal_counts[index] == count:
            return "filtered_unique_terminal_only"
        return "filtered_unique_low_support"
    if evidence.boundary_stops:
        return "read_boundary"
    return "coverage_gap"


def _scan_orientation(
    sequence: str,
    qualities: tuple[int, ...] | None,
    contexts: dict[str, list[int]],
    evidence: list[_Evidence],
    k: int,
    end_window: int,
    low_quality_threshold: int,
) -> None:
    context_length = k - 1
    if len(sequence) < context_length:
        return
    for index in contexts.get(sequence[-context_length:], ()):
        evidence[index].boundary_stops += 1
    if len(sequence) < k:
        return
    for start in range(len(sequence) - k + 1):
        indices = contexts.get(sequence[start : start + context_length])
        if not indices:
            continue
        base = sequence[start + context_length]
        if base not in _BASE_INDEX:
            continue
        base_index = _BASE_INDEX[base]
        quality = (
            None if qualities is None else qualities[start + context_length]
        )
        distance = len(sequence) - 1 - (start + context_length)
        for index in indices:
            item = evidence[index]
            item.counts[base_index] += 1
            if distance < end_window:
                item.terminal_counts[base_index] += 1
            if quality is None:
                item.missing_quality_counts[base_index] += 1
            else:
                item.quality_sums[base_index] += quality
                item.quality_counts[base_index] += 1
                if quality <= low_quality_threshold:
                    item.low_quality_counts[base_index] += 1


def write_dead_end_attribution_report(
    graph: CompactedUnitigGraph,
    reads: list[Read],
    path: Path,
    *,
    min_count: int,
    end_window: int = 5,
    low_quality_threshold: int = 20,
) -> DeadEndAttributionSummary:
    """Explain dead ends from raw continuation and read-boundary evidence.

    This second pass indexes only retained dead-end contexts. It never edits a
    read, graph edge, unitig, or output sequence.
    """
    if min_count < 1:
        raise ValueError("min_count must be at least 1")
    if end_window < 1:
        raise ValueError("end_window must be at least 1")
    dead_ends: list[tuple[int, str, int, str]] = []
    contexts: dict[str, list[int]] = {}
    for unitig_id in range(graph.unitig_count):
        for side, oriented_handle in (
            ("left", (unitig_id << 1) | 1),
            ("right", unitig_id << 1),
        ):
            if graph.outgoing(oriented_handle):
                continue
            context = graph.oriented_sequence(oriented_handle)[-(graph.k - 1) :]
            index = len(dead_ends)
            dead_ends.append((unitig_id, side, oriented_handle, context))
            contexts.setdefault(context, []).append(index)

    evidence = [
        _Evidence([0] * 4, [0] * 4, [0] * 4, [0] * 4, [0] * 4, [0] * 4)
        for _ in dead_ends
    ]
    for read in reads:
        _scan_orientation(
            read.sequence,
            read.qualities,
            contexts,
            evidence,
            graph.k,
            end_window,
            low_quality_threshold,
        )
        _scan_orientation(
            reverse_complement(read.sequence),
            None if read.qualities is None else tuple(reversed(read.qualities)),
            contexts,
            evidence,
            graph.k,
            end_window,
            low_quality_threshold,
        )

    cause_counts: dict[str, int] = {}
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "unitig_id",
        "side",
        "oriented_handle",
        "context",
        "cause",
        "raw_extension_counts",
        "terminal_extension_counts",
        "read_boundary_stops",
        "mean_extension_qualities",
        "low_quality_extension_counts",
        "missing_quality_extension_counts",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for record, item in zip(dead_ends, evidence, strict=True):
            cause = _cause(item, min_count)
            cause_counts[cause] = cause_counts.get(cause, 0) + 1
            means = [
                ""
                if count == 0
                else f"{quality_sum / count:.2f}"
                for quality_sum, count in zip(
                    item.quality_sums, item.quality_counts, strict=True
                )
            ]
            writer.writerow(
                {
                    "unitig_id": record[0] + 1,
                    "side": record[1],
                    "oriented_handle": record[2],
                    "context": record[3],
                    "cause": cause,
                    "raw_extension_counts": ",".join(
                        f"{base}:{count}"
                        for base, count in zip("ACGT", item.counts, strict=True)
                    ),
                    "terminal_extension_counts": ",".join(
                        f"{base}:{count}"
                        for base, count in zip(
                            "ACGT", item.terminal_counts, strict=True
                        )
                    ),
                    "read_boundary_stops": item.boundary_stops,
                    "mean_extension_qualities": ",".join(
                        f"{base}:{mean}"
                        for base, mean in zip("ACGT", means, strict=True)
                    ),
                    "low_quality_extension_counts": ",".join(
                        f"{base}:{count}"
                        for base, count in zip(
                            "ACGT", item.low_quality_counts, strict=True
                        )
                    ),
                    "missing_quality_extension_counts": ",".join(
                        f"{base}:{count}"
                        for base, count in zip(
                            "ACGT", item.missing_quality_counts, strict=True
                        )
                    ),
                }
            )

    return DeadEndAttributionSummary(
        dead_ends=len(dead_ends),
        read_boundaries=cause_counts.get("read_boundary", 0),
        coverage_gaps=cause_counts.get("coverage_gap", 0),
        filtered_unique=sum(
            count
            for cause, count in cause_counts.items()
            if cause.startswith("filtered_unique")
        ),
        filtered_conflicts=cause_counts.get("filtered_conflict", 0),
        retained_context_elsewhere=cause_counts.get(
            "retained_context_elsewhere", 0
        ),
        low_quality=cause_counts.get("filtered_unique_low_quality", 0),
        terminal_only=cause_counts.get("filtered_unique_terminal_only", 0),
        missing_quality=cause_counts.get("filtered_unique_missing_quality", 0),
    )
