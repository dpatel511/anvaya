"""Audit mate-pair support across ambiguous reduced-graph transitions.

Inputs must be synchronized R1/R2 FASTQ files with exactly one record per
molecule in matching order. Coordinates are 0-based and intervals half-open.
"""
import argparse
import hashlib
import inspect
import json
from collections import defaultdict
from pathlib import Path

from anvaya.damage_string_graph import assemble
from anvaya.overlap_assembly import _n50
from anvaya.raw_consensus import DamageProfile
from anvaya.reads import load_reads
from anvaya.sequences import reverse_complement

from branch_transition_opportunity import (
    _all_placements,
    _projected_n50,
    _raw_anchor_index,
    _spell_transition,
    _transition_key,
)


def _pair_name(name):
    token = name.split()[0]
    return token[:-2] if token.endswith(("/1", "/2")) else token


def audit_paired_transitions(contigs, reduced, excluded, left, right):
    reads = [read for pair in zip(left, right) for read in pair]
    profile = DamageProfile((), ())
    anchors = _raw_anchor_index(reads, anchors_per_read=8)
    outgoing = defaultdict(list)
    incoming = defaultdict(list)
    for edge in reduced.values():
        outgoing[edge.source].append(edge)
        incoming[edge.target].append(edge)
    ambiguous = {
        node for node in set(outgoing) | set(incoming)
        if len(outgoing[node]) > 1 or len(incoming[node]) > 1
    }
    transition_support = {}
    ambiguous_pair_placements = 0
    for centre in sorted(ambiguous):
        for first in incoming.get(centre, ()):
            for second in outgoing.get(centre, ()):
                key = _transition_key(first, second)
                if key in transition_support:
                    continue
                sequence, centre_start, centre_end = _spell_transition(
                    contigs, first, second,
                )
                placements = _all_placements(
                    sequence, reads, profile, anchors,
                )
                endpoint_sequences = {
                    min(
                        contigs[node[0]].sequence,
                        reverse_complement(contigs[node[0]].sequence),
                    )
                    for node in (first.source, centre, second.target)
                }
                endpoint_molecules = {
                    read_index // 2 for read_index, read in enumerate(reads)
                    if min(read.sequence, reverse_complement(read.sequence))
                    in endpoint_sequences
                }
                by_read = defaultdict(set)
                for read_index, reverse, offset in placements:
                    by_read[read_index].add((reverse, offset))
                supported = set()
                for molecule in range(len(left)):
                    if molecule in endpoint_molecules:
                        continue
                    first_read = by_read.get(2 * molecule, set())
                    second_read = by_read.get(2 * molecule + 1, set())
                    pairings = set()
                    for first_reverse, first_offset in first_read:
                        for second_reverse, second_offset in second_read:
                            if first_reverse == second_reverse:
                                continue
                            forward_offset, forward_length = (
                                (first_offset, len(left[molecule].sequence))
                                if not first_reverse
                                else (second_offset, len(right[molecule].sequence))
                            )
                            reverse_offset, reverse_length = (
                                (first_offset, len(left[molecule].sequence))
                                if first_reverse
                                else (second_offset, len(right[molecule].sequence))
                            )
                            if forward_offset >= reverse_offset:
                                continue
                            if (
                                forward_offset < centre_start
                                and reverse_offset + reverse_length > centre_end
                            ):
                                pairings.add((
                                    first_reverse, first_offset,
                                    second_reverse, second_offset,
                                ))
                    if len(pairings) == 1:
                        supported.add(molecule)
                    elif len(pairings) > 1:
                        ambiguous_pair_placements += 1
                transition_support[key] = supported

    by_molecule = defaultdict(set)
    for key, molecules in transition_support.items():
        for molecule in molecules:
            by_molecule[molecule].add(key)
    conflicting = {molecule for molecule, keys in by_molecule.items() if len(keys) > 1}
    for molecules in transition_support.values():
        molecules.difference_update(conflicting)

    selected = set()
    for branch in {key[1][0] for key in transition_support}:
        keys = [key for key in transition_support if key[1][0] == branch]
        supported = [key for key in keys if len(transition_support[key]) >= 2]
        if len(supported) == 1 and all(
            not transition_support[key] for key in keys if key != supported[0]
        ):
            selected.add(supported[0])
    return {
        "branch_physical_nodes": len({node[0] for node in ambiguous}),
        "possible_transitions": len(transition_support),
        "transitions_with_1_pair": sum(bool(v) for v in transition_support.values()),
        "transitions_with_2_pairs": sum(len(v) >= 2 for v in transition_support.values()),
        "uniquely_resolved_transitions": len(selected),
        "conflicting_molecules": len(conflicting),
        "ambiguous_pair_placements": ambiguous_pair_placements,
        "projection_baseline_n50": _projected_n50(contigs, reduced, excluded, set()),
        "projection_only_n50": _projected_n50(contigs, reduced, excluded, selected),
    }


def run(left_path, right_path, output, maximum_pairs):
    if output.exists():
        raise ValueError("output already exists")
    left = load_reads(left_path, maximum_reads=maximum_pairs)
    right = load_reads(right_path, maximum_reads=maximum_pairs)
    if len(left) != len(right):
        raise ValueError("paired FASTQ files contain different record counts")
    if len(left) != maximum_pairs:
        raise ValueError(
            f"paired audit requires exactly {maximum_pairs} records per input"
        )
    mismatched = [
        index for index, pair in enumerate(zip(left, right))
        if _pair_name(pair[0].name) != _pair_name(pair[1].name)
    ]
    if mismatched:
        raise ValueError(f"paired FASTQ names disagree at record {mismatched[0] + 1}")
    reads = [read for pair in zip(left, right) for read in pair]
    molecules = [index for index in range(len(left)) for _ in range(2)]
    audit = {}

    def capture(contigs, reduced, excluded):
        audit.update(audit_paired_transitions(
            contigs, reduced, excluded, left, right,
        ))

    pool, graph = assemble(
        reads, molecule_ids=molecules, maximum_reads=len(reads),
        reduced_graph_audit=capture,
    )
    lengths = [len(record.current.sequence) for record in pool.active_derived]
    baseline_n50 = _n50(lengths)
    if audit["projection_baseline_n50"] != baseline_n50:
        raise ValueError("paired projection audit does not reproduce baseline N50")
    result = {
        "evidence": "reference_free_exact_mate_threading",
        "coordinate_system": "0-based; intervals half-open",
        "pairs": len(left),
        "reads": len(reads),
        "baseline_contigs": len(lengths),
        "baseline_n50": baseline_n50,
        "graph": graph,
        "audit": audit,
        "source_sha256": {
            "left": hashlib.sha256(Path(left_path).read_bytes()).hexdigest(),
            "right": hashlib.sha256(Path(right_path).read_bytes()).hexdigest(),
            **{
                str(path.relative_to(Path(__file__).resolve().parents[1])).replace("\\", "/"):
                hashlib.sha256(path.read_bytes()).hexdigest()
                for path in {
                    Path(__file__).resolve(),
                    Path(inspect.getfile(assemble)).resolve(),
                    Path(inspect.getfile(_projected_n50)).resolve(),
                }
            },
        },
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "paired-branch-opportunity.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8",
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left", required=True, type=Path)
    parser.add_argument("--right", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--maximum-pairs", type=int, default=25000)
    args = parser.parse_args()
    if args.maximum_pairs < 1:
        parser.error("maximum pairs must be positive")
    result = run(args.left, args.right, args.output, args.maximum_pairs)
    print(json.dumps({
        "pairs": result["pairs"],
        "baseline_n50": result["baseline_n50"],
        **result["audit"],
    }))


if __name__ == "__main__":
    main()
