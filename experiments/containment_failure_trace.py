"""Trace the consumed containment validation failure by immutable raw provenance.

All coordinates are 0-based and intervals are half-open. Reference truth is used
only after assembly to select and describe the target coordinate.
"""
import argparse
import hashlib
import inspect
import json
import sys
from dataclasses import asdict
from pathlib import Path

from anvaya.damage_string_graph import assemble
from anvaya.raw_consensus import DamageProfile, project_raw_consensus
from anvaya.sequences import reverse_complement

from controlled_assembly import score
from variable_length_validation import fixture


SEED = 96012
SCENARIO = "damage_error_strains_80_20"
REFERENCE = 1
POSITION_0BASED = 1539
DAMAGE = 0.25
ERROR_RATE = 0.002


def _resolved_transform(evidence):
    if "sequence_resolved_transform" in evidence:
        return evidence["sequence_resolved_transform"], "sequence_resolved"
    if len(evidence["transforms"]) == 1:
        return evidence["transforms"][0], "unique_origin"
    return None, None


def _target_outputs(pool, sequences, metrics, references, reads, origins):
    outputs = []
    for record, sequence, evidence in zip(pool.active_derived, sequences, metrics["contig_evidence"]):
        transform, resolution = _resolved_transform(evidence)
        if transform is None or transform["reference"] != REFERENCE:
            continue
        oriented_position = (
            len(references[REFERENCE]) - 1 - POSITION_0BASED
            if transform["reverse"] else POSITION_0BASED
        )
        contig_position = oriented_position - transform["start_0based"]
        if not 0 <= contig_position < len(sequence):
            continue
        base = sequence[contig_position]
        reference_base = reverse_complement(base) if transform["reverse"] else base
        observations = []
        for placement in record.raw_placements or ():
            position_in_oriented_read = contig_position - placement.offset
            raw_length = len(reads[placement.read_index].sequence)
            raw_position = (
                raw_length - 1 - position_in_oriented_read
                if placement.reverse else position_in_oriented_read
            )
            if not placement.read_start <= raw_position < placement.read_stop:
                continue
            read = reads[placement.read_index]
            raw_base = read.sequence[raw_position]
            contig_oriented_base = reverse_complement(raw_base) if placement.reverse else raw_base
            origin_reference, start, stop, origin_reverse = origins[placement.read_index]
            latent = references[origin_reference][start:stop]
            if origin_reverse:
                latent = reverse_complement(latent)
            latent_base = latent[raw_position]
            quality = read.qualities[raw_position] if read.qualities is not None else None
            observations.append({
                "raw_index": placement.read_index,
                "raw_name": read.name,
                "molecule_id": pool.records[placement.read_index].molecule_id,
                "raw_position_0based": raw_position,
                "distance_5p": raw_position,
                "distance_3p": raw_length - 1 - raw_position,
                "raw_base": raw_base,
                "latent_raw_oriented_base": latent_base,
                "contig_oriented_base": contig_oriented_base,
                "quality": quality,
                "difference_class": (
                    "none" if raw_base == latent_base else
                    "sequencing_error_or_post_damage_error" if quality == 10 else
                    "simulated_terminal_damage"
                ),
                "placement": asdict(placement),
                "origin": {
                    "reference": origin_reference, "start_0based": start,
                    "stop_exclusive": stop, "reverse": origin_reverse,
                },
            })
        outputs.append({
            "record_name": record.current.name,
            "sequence": sequence,
            "contig_position_0based": contig_position,
            "reference_oriented_base": reference_base,
            "resolution": resolution,
            "transform": transform,
            "raw_read_names": sorted({
                reads[placement.read_index].name for placement in record.raw_placements or ()
            }),
            "target_observations": observations,
        })
    return outputs


def _output_catalog(pool, sequences, metrics, references, reads):
    def target_candidates(sequence, evidence):
        candidates = []
        for transform in evidence["transforms"]:
            if transform["reference"] != REFERENCE:
                continue
            oriented_position = (
                len(references[REFERENCE]) - 1 - POSITION_0BASED
                if transform["reverse"] else POSITION_0BASED
            )
            contig_position = oriented_position - transform["start_0based"]
            if 0 <= contig_position < len(sequence):
                base = sequence[contig_position]
                candidates.append({
                    "contig_position_0based": contig_position,
                    "reference_oriented_base": (
                        reverse_complement(base) if transform["reverse"] else base
                    ),
                    "transform": transform,
                })
        return candidates

    return [
        {
            "record_name": record.current.name,
            "sequence": sequence,
            "raw_read_names": sorted({
                reads[placement.read_index].name
                for placement in record.raw_placements or ()
            }),
            "provenance_category": evidence["provenance_category"],
            "sequence_category": evidence["sequence_category"],
            "transforms": evidence["transforms"],
            "sequence_resolved_transform": evidence.get("sequence_resolved_transform"),
            "target_candidate_projections": target_candidates(sequence, evidence),
            "raw_placements": [asdict(placement) for placement in record.raw_placements or ()],
        }
        for record, sequence, evidence
        in zip(pool.active_derived, sequences, metrics["contig_evidence"])
    ]


def _capture_assembly_decisions(reads):
    captured = {}

    def trace(frame, event, _argument):
        if event != "return" or frame.f_code is not assemble.__code__:
            return trace
        local = frame.f_locals
        sequences = local.get("sequences", ())
        groups = local.get("groups", {})

        def raw_names(node):
            return sorted(reads[index].name for index, _ in groups.get(sequences[node], ()))

        def selected_paths(name):
            return [
                {
                    "path": path,
                    "edges": [
                        {
                            "source_node": edge.source[0],
                            "source_reverse": edge.source[1],
                            "source_raw_names": raw_names(edge.source[0]),
                            "target_node": edge.target[0],
                            "target_reverse": edge.target[1],
                            "target_raw_names": raw_names(edge.target[0]),
                            "shift": edge.shift,
                            "overlap": edge.overlap,
                            "corrections": edge.corrections,
                        }
                        for edge in edges
                    ],
                }
                for path, edges in sorted(local.get(name, {}).items())
            ]

        uniquely_owned = set(local.get("uniquely_owned", ()))
        captured["selected_paths"] = selected_paths("selected_edges")
        captured["unique_containment_baseline_selected_paths"] = selected_paths(
            "baseline_selected"
        )
        captured["containments"] = [
            {
                "child_node": child,
                "child_raw_names": raw_names(child),
                "ownership": "post_projection_unique" if child in uniquely_owned else "exact_unique",
                "parents": [
                    {
                        "parent_node": parent,
                        "parent_raw_names": raw_names(parent),
                        "child_reverse": reverse,
                        "offset_0based": offset,
                    }
                    for parent, reverse, offset in options
                ],
            }
            for child, options in sorted(local.get("containment_options", {}).items())
        ]
        captured["correction_protected_child_nodes"] = sorted(
            local.get("correction_protected", ())
        )
        return trace

    return captured, trace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--implementation-label", required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        parser.error("output already exists")

    references, reads, origins, events, errors = fixture(
        SEED, "strains", DAMAGE, ERROR_RATE, include_events=True,
    )
    rates = tuple(DAMAGE * 0.65 ** position for position in range(5))
    profile = DamageProfile(rates, rates)
    assembly_decisions, trace = _capture_assembly_decisions(reads)
    sys.settrace(trace)
    try:
        pool, graph = assemble(reads, damage_profile=profile)
    finally:
        sys.settrace(None)
    polished, consensus = project_raw_consensus(pool, profile)
    layout_sequences = [record.current.sequence for record in pool.active_derived]
    consensus_sequences = [read.sequence for read in polished]
    layout_metrics = score(
        pool, references, origins, damage_events=events,
        include_coordinate_states=True,
    )
    consensus_metrics = score(
        pool, references, origins, polished, damage_events=events,
        include_coordinate_states=True,
    )
    input_payload = [
        (read.name, read.sequence, read.qualities, origins[index])
        for index, read in enumerate(reads)
    ]
    source_path = Path(inspect.getfile(assemble)).resolve()
    output = {
        "status": "completed_consumed_validation_diagnosis",
        "seed_role": "consumed_validation_diagnosis",
        "implementation_label": arguments.implementation_label,
        "coordinate_system": "0-based reference positions; intervals half-open",
        "target": {
            "seed": SEED, "scenario": SCENARIO, "reference": REFERENCE,
            "position_0based": POSITION_0BASED,
            "reference_base": references[REFERENCE][POSITION_0BASED],
        },
        "input_sha256": hashlib.sha256(
            json.dumps(input_payload, sort_keys=True).encode()
        ).hexdigest(),
        "simulated_damage_events": len(events),
        "simulated_sequencing_errors": len(errors),
        "assembly_source": str(source_path),
        "assembly_source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "graph": graph,
        "assembly_decisions": assembly_decisions,
        "consensus_diagnostics": asdict(consensus),
        "layout": {
            "coordinate_state": next(
                (state for state in layout_metrics["nonduplicated_accuracy_coordinate_states"]
                 if state["reference"] == REFERENCE
                 and state["position_0based"] == POSITION_0BASED), None,
            ),
            "outputs": _target_outputs(
                pool, layout_sequences, layout_metrics, references, reads, origins,
            ),
            "output_catalog": _output_catalog(
                pool, layout_sequences, layout_metrics, references, reads,
            ),
        },
        "consensus": {
            "coordinate_state": next(
                (state for state in consensus_metrics["nonduplicated_accuracy_coordinate_states"]
                 if state["reference"] == REFERENCE
                 and state["position_0based"] == POSITION_0BASED), None,
            ),
            "outputs": _target_outputs(
                pool, consensus_sequences, consensus_metrics, references, reads, origins,
            ),
            "output_catalog": _output_catalog(
                pool, consensus_sequences, consensus_metrics, references, reads,
            ),
        },
        "limitations": [
            "The seed is consumed validation and this is diagnosis only",
            "Unitig names are deliberately excluded from cross-run joins",
            "Target selection uses private synthetic truth after assembly",
        ],
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
