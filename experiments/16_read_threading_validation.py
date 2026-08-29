"""Validate reciprocal read-threaded junctions against sequence truth."""

import argparse
import csv
import importlib.util
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

from anvaya.bidirected import build_bidirected_dbg
from anvaya.paired_extension import spell_extended_paths
from anvaya.read_threading import (
    audit_read_threads,
    resolve_read_thread_extensions,
)
from anvaya.reads import Read
from anvaya.sequences import reverse_complement
from anvaya.unitig_graph import build_compacted_unitig_graph


def _load_previous_experiment():
    path = Path(__file__).with_name("06_tip_matching_validation.py")
    spec = importlib.util.spec_from_file_location("tip_validation", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


previous = _load_previous_experiment()


@dataclass(frozen=True)
class RunRecord:
    condition: str
    seed: int
    reads: int
    junctions: int
    candidate_transitions: int
    supported_transitions: int
    strong_winner_junctions: int
    reciprocal_links: int
    joined_links: int
    repeat_guarded_links: int
    true_joined_links: int
    false_joined_links: int
    ry_true_joined_links: int
    ry_false_joined_links: int
    true_links: int
    false_links: int
    ry_true_links: int
    ry_false_links: int
    precision: float
    baseline_contigs: int
    false_baseline_contigs: int
    extended_contigs: int
    false_extended_contigs: int


def _contains(reference: str, sequence: str) -> bool:
    return sequence in reference or reverse_complement(sequence) in reference


_RY_TRANSLATION = str.maketrans({"A": "R", "G": "R", "C": "Y", "T": "Y"})


def _contains_ry(encoded_reference: str, sequence: str) -> bool:
    encoded_sequence = sequence.translate(_RY_TRANSLATION)
    reverse = reverse_complement(sequence).translate(_RY_TRANSLATION)
    return encoded_sequence in encoded_reference or reverse in encoded_reference


def _physical_link(source: int, target: int) -> tuple[int, int]:
    return min((source, target), (target ^ 1, source ^ 1))


def _evaluate(
    condition: str,
    seed: int,
    simulated_reads,
    truth_references: tuple[str, ...],
) -> RunRecord:
    reads = [
        Read(read.name, read.sequence, (30,) * len(read.sequence))
        for read in simulated_reads
    ]
    graph = build_bidirected_dbg(
        [read.sequence for read in reads],
        previous.K,
        previous.MIN_COUNT,
    )
    unitigs = build_compacted_unitig_graph(graph, track_edge_handles=True)
    rows, summary = audit_read_threads(
        graph,
        unitigs,
        ((read,) for read in reads),
        end_window=previous.END_WINDOW,
    )

    accepted = {
        _physical_link(row.source, row.target): row
        for row in rows
        if row.reciprocal
    }
    ry_references = tuple(
        reference.translate(_RY_TRANSLATION) for reference in truth_references
    )
    true_physical_links: set[tuple[int, int]] = set()
    ry_true_physical_links: set[tuple[int, int]] = set()
    for physical, row in accepted.items():
        sequence = unitigs.oriented_sequence(row.source)
        sequence += unitigs.oriented_sequence(row.target)[previous.K - 1 :]
        if any(_contains(reference, sequence) for reference in truth_references):
            true_physical_links.add(physical)
        if any(_contains_ry(reference, sequence) for reference in ry_references):
            ry_true_physical_links.add(physical)
    true_links = len(true_physical_links)
    false_links = len(accepted) - true_links
    successors = resolve_read_thread_extensions(rows)
    joined_physical_links = {
        _physical_link(source, target)
        for source, target in successors.items()
    }
    extended = spell_extended_paths(unitigs, successors)
    false_baseline = sum(
        not any(_contains(reference, sequence) for reference in truth_references)
        for sequence in unitigs.sequences
    )
    false_extended = sum(
        not any(_contains(reference, sequence) for reference in truth_references)
        for sequence in extended
    )
    return RunRecord(
        condition=condition,
        seed=seed,
        reads=len(reads),
        junctions=summary.junctions,
        candidate_transitions=summary.candidate_transitions,
        supported_transitions=summary.supported_transitions,
        strong_winner_junctions=summary.resolvable_junctions,
        reciprocal_links=len(accepted),
        joined_links=len(successors) // 2,
        repeat_guarded_links=len(accepted) - len(successors) // 2,
        true_joined_links=len(joined_physical_links & true_physical_links),
        false_joined_links=len(joined_physical_links - true_physical_links),
        ry_true_joined_links=len(joined_physical_links & ry_true_physical_links),
        ry_false_joined_links=len(joined_physical_links - ry_true_physical_links),
        true_links=true_links,
        false_links=false_links,
        ry_true_links=len(ry_true_physical_links),
        ry_false_links=len(accepted) - len(ry_true_physical_links),
        precision=(true_links / len(accepted) if accepted else 1.0),
        baseline_contigs=unitigs.unitig_count,
        false_baseline_contigs=false_baseline,
        extended_contigs=len(extended),
        false_extended_contigs=false_extended,
    )


def run_matrix(seeds: tuple[int, ...]) -> list[RunRecord]:
    """Run clean, damage, error, and related-strain conditions."""
    records: list[RunRecord] = []
    profile = previous.DAMAGE_PROFILES["standard"]
    for seed in seeds:
        reference = previous.random_reference(seed)
        source = previous.simulate_fragments(
            reference,
            previous.READ_LENGTH,
            previous.TOTAL_COVERAGE,
            seed + 10_000,
            reverse_fraction=0.5,
        )
        damaged = previous.add_terminal_damage(
            source,
            profile,
            profile,
            seed + 20_000,
        )
        errors = previous.add_sequencing_errors(
            source,
            0.005,
            seed + 30_000,
        )
        rare_reference, _ = previous.rare_reference(reference, seed + 40_000)
        major_reads = previous.simulate_fragments(
            reference,
            previous.READ_LENGTH,
            15,
            seed + 50_000,
            reverse_fraction=0.5,
        )
        rare_reads = previous.simulate_fragments(
            rare_reference,
            previous.READ_LENGTH,
            5,
            seed + 60_000,
            reverse_fraction=0.5,
        )
        for condition, reads, truth in (
            ("clean", source, (reference,)),
            ("damage", damaged, (reference,)),
            ("error", errors, (reference,)),
            ("rare_strain", major_reads + rare_reads, (reference, rare_reference)),
        ):
            records.append(_evaluate(condition, seed, reads, truth))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seeds", nargs="+", type=int, default=range(1501, 1506))
    arguments = parser.parse_args()
    records = run_matrix(tuple(arguments.seeds))
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RunRecord.__dataclass_fields__)
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)
    total_links = sum(record.reciprocal_links for record in records)
    total_true = sum(record.true_links for record in records)
    print(f"runs={len(records)}")
    print(f"reciprocal_links={total_links}")
    print(f"true_links={total_true}")
    print(f"false_links={total_links - total_true}")
    print(f"precision={total_true / total_links if total_links else 1.0:.6f}")
    print(f"output={arguments.output}")


if __name__ == "__main__":
    main()
