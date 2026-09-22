"""Trace clean truth-known overlaps through the bounded string graph.

Truth coordinates are 0-based, half-open and are used only by this evaluator.
The assembler receives reads, qualities and molecule IDs only.
"""
import argparse
import hashlib
import importlib.util
import json
import platform
import resource
import shutil
import sys
import time
from collections import Counter, defaultdict, deque
from dataclasses import asdict
from pathlib import Path

from anvaya.damage_string_graph import assemble
from anvaya.overlap_assembly import (
    _Alignment, _MasterOverlapEdge, _candidate_alignments, _n50,
)
from anvaya.overlap_graph import _project_master_edges
from anvaya.overlap_index import _anchor_index
from anvaya.overlap_progressive_links import _add_bidirected_edge
from anvaya.reads import Read
from anvaya.sequences import reverse_complement



def _load_experiment(name):
    if name in sys.modules:
        return sys.modules[name]
    path = Path(__file__).with_name(f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


_load_experiment("controlled_assembly")
fixture = _load_experiment("variable_length_validation").fixture


def _reverse(node):
    return node[0], not node[1]


def _physical(source, target):
    return min((source, target), (_reverse(target), _reverse(source)))


def _edge_from_alignment(nodes, source_index, candidate):
    source_length = len(nodes[source_index].sequence)
    target_length = len(candidate.sequence)
    left = max(0, -candidate.offset)
    right = max(0, candidate.offset + target_length - source_length)
    target = (candidate.read_index, bool(candidate.reverse))
    if right and not left:
        return _MasterOverlapEdge(
            (source_index, False), target, candidate.offset,
            source_length - candidate.offset,
        )
    if left and not right:
        return _MasterOverlapEdge(
            target, (source_index, False), -candidate.offset,
            target_length + candidate.offset,
        )
    return None


def _candidate_edges(nodes, minimum_overlap=30):
    anchors = _anchor_index(nodes, 15, 0, 8, 100)
    indices = list(range(len(nodes)))
    position_bits = max(len(node.sequence) for node in nodes).bit_length()
    target_window = max(len(node.sequence) for node in nodes)
    edges, ambiguous, stages = {}, set(), {}
    stage_rank = {
        "anchor_support_rejected": 0, "short_overlap": 1,
        "dna_identity_rejected": 2, "ry_identity_rejected": 3,
        "molecule_selection_rejected": 4, "selected": 5,
    }
    for source_index, source in enumerate(nodes):
        trace = {}
        candidates = _candidate_alignments(
            source.sequence, nodes, indices, anchors, {source_index},
            anchor_k=15, anchors_per_read=8, maximum_anchor_occurrences=100,
            minimum_anchor_matches=1, minimum_overlap=minimum_overlap,
            minimum_identity=1.0, minimum_ry_identity=1.0,
            position_bits=position_bits, target_window=target_window,
            stage_trace=trace,
        )
        for (target_index, reverse, offset), stage in trace.items():
            candidate_sequence = nodes[target_index].sequence
            if reverse:
                candidate_sequence = reverse_complement(candidate_sequence)
            candidate = _Alignment(target_index, candidate_sequence, offset, 0, reverse)
            edge = _edge_from_alignment(nodes, source_index, candidate)
            if edge is not None:
                key = _physical(edge.source, edge.target)
                if stage_rank[stage] > stage_rank.get(stages.get(key, ""), -1):
                    stages[key] = stage
        for candidate in candidates:
            edge = _edge_from_alignment(nodes, source_index, candidate)
            if edge is not None:
                _add_bidirected_edge(edges, nodes, edge, ambiguous)
    return edges, ambiguous, stages


def _canonical_truth(reads, origins):
    groups = {}
    for raw_index, read in enumerate(reads):
        canonical = min(read.sequence, reverse_complement(read.sequence))
        groups.setdefault(canonical, []).append(raw_index)
    sequences = list(groups)
    placements = defaultdict(set)
    for node, sequence in enumerate(sequences):
        for raw_index in groups[sequence]:
            reference, start, stop, raw_reverse = origins[raw_index]
            canonical_reverse = reads[raw_index].sequence != sequence
            placements[node].add(
                (reference, start, stop, bool(raw_reverse) != canonical_reverse)
            )
    return sequences, groups, placements


def _truth_edges(placements, minimum_overlap=30):
    proper = defaultdict(list)
    boundary_containment = defaultdict(list)
    nodes = sorted(placements)
    for left in nodes:
        for right in nodes:
            if left == right:
                continue
            for first in placements[left]:
                for second in placements[right]:
                    ref1, start1, stop1, reverse1 = first
                    ref2, start2, stop2, reverse2 = second
                    overlap = stop1 - start2
                    if (ref1 == ref2 and start1 <= start2 and stop1 <= stop2
                            and (start1 < start2 or stop1 < stop2)
                            and overlap >= minimum_overlap):
                        source, target = (left, reverse1), (right, reverse2)
                        key = _physical(source, target)
                        destination = (proper if start1 < start2 and stop1 < stop2
                                       else boundary_containment)
                        destination[key].append({
                            "source": list(source), "target": list(target),
                            "reference": ref1, "source_interval": [start1, stop1],
                            "target_interval": [start2, stop2],
                            "shift": start2 - start1, "overlap": overlap,
                        })
    return proper, boundary_containment


def _containments(sequences):
    candidates = defaultdict(set)
    for child_index, child in enumerate(sequences):
        orientations = [(False, child)]
        reverse = reverse_complement(child)
        if reverse != child:
            orientations.append((True, reverse))
        for parent_index, parent in enumerate(sequences):
            if len(child) >= len(parent):
                continue
            for child_reverse, oriented in orientations:
                offset = parent.find(oriented)
                while offset >= 0:
                    candidates[child_index].add((parent_index, child_reverse, offset))
                    offset = parent.find(oriented, offset + 1)
    unique = {child: next(iter(options)) for child, options in candidates.items()
              if len(options) == 1}
    repeated_parent = {
        child for child, options in candidates.items()
        if any(count > 1 for count in Counter(parent for parent, _, _ in options).values())
    }
    return candidates, unique, repeated_parent


def _compatible_containment(child, candidate, placements, sequences):
    parent, child_reverse, offset = candidate
    length = len(sequences[child])
    for ref, parent_start, parent_stop, parent_reverse in placements[parent]:
        if parent_reverse:
            start = parent_stop - offset - length
            reverse = not child_reverse
        else:
            start = parent_start + offset
            reverse = child_reverse
        if (ref, start, start + length, reverse) in placements[child]:
            return ref, start, start + length, reverse
    return None


def _reduce(edges):
    outgoing = defaultdict(list)
    for edge in edges.values():
        outgoing[edge.source].append(edge)
    transitive = set()
    for edge in edges.values():
        for first in outgoing[edge.source]:
            if first.target == edge.target:
                continue
            for second in outgoing.get(first.target, ()):
                if second.target == edge.target and first.shift + second.shift == edge.shift:
                    transitive.add((edge.source, edge.target))
                    break
            if (edge.source, edge.target) in transitive:
                break
    reduced = {key: edge for key, edge in edges.items() if key not in transitive}
    outgoing, incoming = defaultdict(list), defaultdict(list)
    for edge in reduced.values():
        outgoing[edge.source].append(edge)
        incoming[edge.target].append(edge)
    reciprocal = {key: edge for key, edge in reduced.items()
                  if len(outgoing[edge.source]) == 1 and len(incoming[edge.target]) == 1}
    return transitive, reduced, reciprocal, outgoing, incoming


def _oriented_sequence(nodes, node):
    sequence = nodes[node[0]].sequence
    return reverse_complement(sequence) if node[1] else sequence


def _spell_edges(nodes, path):
    """Spell one directed path, including its accepted overlap corrections."""
    if not path:
        return None
    sequence = list(_oriented_sequence(nodes, path[0].source))
    current = path[0].source
    for edge in path:
        if edge.source != current:
            return None
        source = _oriented_sequence(nodes, edge.source)
        source_start = len(sequence) - len(source)
        if source_start < 0 or source_start + edge.shift != len(sequence) - edge.overlap:
            return None
        for position, base in edge.corrections:
            assembled_position = source_start + position
            if not 0 <= assembled_position < len(sequence):
                return None
            sequence[assembled_position] = base
        sequence.extend(_oriented_sequence(nodes, edge.target)[edge.overlap:])
        current = edge.target
    return "".join(sequence)


def _bounded_equivalent_witnesses(nodes, edges, direct, max_edges=8, max_states=10000):
    """Return at most two simple, positive-shift paths equivalent to ``direct``."""
    outgoing = defaultdict(list)
    for key, edge in edges.items():
        if key != (direct.source, direct.target) and edge.shift > 0:
            outgoing[edge.source].append(edge)
    expected = _spell_edges(nodes, (direct,))
    queue = deque([(direct.source, (), frozenset((direct.source,)), 0)])
    witnesses = []
    visited_states = 0
    while queue and visited_states < max_states and len(witnesses) < 2:
        node, path, visited, shift = queue.popleft()
        visited_states += 1
        if len(path) >= max_edges:
            continue
        for edge in outgoing.get(node, ()):
            total_shift = shift + edge.shift
            if total_shift > direct.shift or edge.target in visited:
                continue
            candidate = (*path, edge)
            if edge.target == direct.target:
                if total_shift == direct.shift and _spell_edges(nodes, candidate) == expected:
                    witnesses.append(candidate)
                continue
            queue.append((edge.target, candidate, visited | {edge.target}, total_shift))
    return witnesses, visited_states, bool(queue)


def _truth_joined(layouts, truth, sequences):
    source = tuple(truth["source"])
    target = tuple(truth["target"])
    overlap = truth["overlap"]
    mirror_source = (target[0], not target[1])
    mirror_target = (source[0], not source[1])
    mirror_shift = len(sequences[target[0]]) - overlap
    for layout in layouts.values():
        positions = {(node, reverse): offset for node, reverse, offset in layout}
        if source in positions and target in positions:
            if positions[target] - positions[source] == truth["shift"]:
                return True
        if mirror_source in positions and mirror_target in positions:
            if positions[mirror_target] - positions[mirror_source] == mirror_shift:
                return True
    return False


def _project_variant(nodes, edges, expected, sequences, name, excluded):
    layouts = {}
    projected, diagnostics = _project_master_edges(
        nodes, edges, name_prefix=name, path_layouts=layouts,
        verify_transitive_alleles=True, excluded_contigs=excluded,
    )
    active_expected = {
        key: truths for key, truths in expected.items()
        if key[0][0] not in excluded and key[1][0] not in excluded
    }
    joined = sum(
        any(_truth_joined(layouts, truth, sequences) for truth in truths)
        for truths in active_expected.values()
    )
    lengths = [len(read.sequence) for read in projected]
    return {
        "contigs": len(projected), "n50": _n50(lengths),
        "longest": max(lengths, default=0),
        "truth_dovetails_joined": joined,
        "truth_dovetails_total": len(active_expected),
        "diagnostics": asdict(diagnostics),
    }


def _maximal_nonbranching_variant(nodes, edges, expected, sequences, excluded):
    """Diagnostic maximal non-branching paths; does not alter the assembler."""
    _, reduced, _, outgoing, incoming = _reduce(edges)
    paths = []
    visited_edges = set()
    starts = sorted(
        node for node in set(outgoing) | set(incoming)
        if len(incoming[node]) != 1 or len(outgoing[node]) != 1
    )
    for start in starts:
        for first in sorted(outgoing.get(start, ()), key=lambda edge: edge.target):
            key = (first.source, first.target)
            if key in visited_edges:
                continue
            path = [first]
            visited_edges.add(key)
            seen = {first.source, first.target}
            node = first.target
            while len(incoming[node]) == 1 and len(outgoing[node]) == 1:
                edge = outgoing[node][0]
                key = (edge.source, edge.target)
                if key in visited_edges or edge.target in seen:
                    break
                path.append(edge)
                visited_edges.add(key)
                seen.add(edge.target)
                node = edge.target
            paths.append(tuple(path))

    canonical_paths = {}
    layouts = {}
    used_nodes = set()
    for path in paths:
        spelled = _spell_edges(nodes, path)
        if spelled is None:
            continue
        layout = [(path[0].source[0], path[0].source[1], 0)]
        offset = 0
        for edge in path:
            offset += edge.shift
            layout.append((edge.target[0], edge.target[1], offset))
        canonical = min(spelled, reverse_complement(spelled))
        used_nodes.update(node for node, _, _ in layout)
        if canonical not in canonical_paths:
            canonical_paths[canonical] = spelled
            layouts[f"unitig_{len(layouts) + 1}"] = tuple(layout)
    singletons = [
        node.sequence for index, node in enumerate(nodes)
        if index not in used_nodes and index not in excluded
    ]
    lengths = [len(sequence) for sequence in canonical_paths] + [len(s) for s in singletons]
    active_expected = {
        key: truths for key, truths in expected.items()
        if key[0][0] not in excluded and key[1][0] not in excluded
    }
    joined = sum(
        any(_truth_joined(layouts, truth, sequences) for truth in truths)
        for truths in active_expected.values()
    )
    return {
        "status": "diagnostic_only",
        "contigs": len(lengths), "n50": _n50(lengths),
        "longest": max(lengths, default=0),
        "path_sequences": len(canonical_paths),
        "nodes_reused_across_paths": sum(
            count > 1 for count in Counter(
                node for layout in layouts.values() for node, _, _ in layout
            ).values()
        ),
        "truth_dovetails_joined": joined,
        "truth_dovetails_total": len(active_expected),
        "unspelled_cycle_directed_edges": len(reduced) - len(visited_edges),
    }


def _phase_a_ablation(nodes, filtered, expected, sequences, contained_nodes):
    strict = {
        key: edge for key, edge in filtered.items()
        if 0 < edge.overlap < min(
            len(sequences[edge.source[0]]), len(sequences[edge.target[0]])
        )
    }
    boundary = set(filtered) - set(strict)
    _, reduced, _, outgoing, incoming = _reduce(filtered)
    residual_branch = {
        key: edge for key, edge in reduced.items()
        if len(outgoing[edge.source]) > 1 or len(incoming[edge.target]) > 1
    }

    removable = set()
    protected = set()
    witness_lengths = Counter()
    witness_statuses = Counter()
    search_limited = 0
    seen = set()
    candidates = sorted(
        residual_branch.values(), key=lambda edge: (-edge.shift, edge.source, edge.target)
    )
    for edge in candidates:
        physical = _physical(edge.source, edge.target)
        if physical in seen:
            continue
        seen.add(physical)
        if physical in protected or physical in removable:
            continue
        if (edge.source, edge.target) not in strict:
            witness_statuses["boundary_relation"] += 1
            continue
        witnesses, _, limited = _bounded_equivalent_witnesses(nodes, strict, edge)
        search_limited += limited
        if limited:
            witness_statuses["search_limit"] += 1
            continue
        if not witnesses:
            witness_statuses["no_equivalent_witness"] += 1
            continue
        if len(witnesses) > 1:
            witness_statuses["multiple_equivalent_witnesses"] += 1
            continue
        mirror = (_reverse(edge.target), _reverse(edge.source))
        mirror_witnesses, _, mirror_limited = _bounded_equivalent_witnesses(
            nodes, strict, strict[mirror]
        )
        search_limited += mirror_limited
        if mirror_limited:
            witness_statuses["mirror_search_limit"] += 1
            continue
        if len(mirror_witnesses) != 1:
            witness_statuses["mirror_not_unique"] += 1
            continue
        removable.add(physical)
        witness_statuses["unique_bidirectional_witness"] += 1
        for witness_edge in (*witnesses[0], *mirror_witnesses[0]):
            protected.add(_physical(witness_edge.source, witness_edge.target))
        witness_lengths[len(witnesses[0])] += 1

    without_boundary = {key: edge for key, edge in filtered.items() if key not in boundary}
    without_redundant = {
        key: edge for key, edge in filtered.items()
        if _physical(edge.source, edge.target) not in removable
    }
    return {
        "directed_edges": len(filtered),
        "physical_edges": len(filtered) // 2,
        "strict_directed_edges": len(strict),
        "boundary_directed_edges": len(boundary),
        "residual_branch_directed_edges": len(residual_branch),
        "verified_redundant_physical_edges": len(removable),
        "protected_witness_physical_edges": len(protected),
        "witness_path_edge_counts": dict(sorted(witness_lengths.items())),
        "residual_branch_physical_edge_statuses": dict(sorted(witness_statuses.items())),
        "bounded_search_limit_hits": search_limited,
        "baseline": _project_variant(
            nodes, filtered, expected, sequences, "baseline", contained_nodes
        ),
        "without_boundary_edges": _project_variant(
            nodes, without_boundary, expected, sequences, "no_boundary", contained_nodes
        ),
        "without_verified_redundant_edges": _project_variant(
            nodes, without_redundant, expected, sequences, "no_redundant", contained_nodes
        ),
        "maximal_nonbranching_paths": _maximal_nonbranching_variant(
            nodes, filtered, expected, sequences, contained_nodes
        ),
    }


def containment_branch_fixture(reverse=False):
    """Four clean reads where an ambiguously contained read creates a branch."""
    import random
    rng = random.Random(95302)
    reference = "".join(rng.choice("ACGT") for _ in range(160))
    reads = [Read(name, reference[start:stop]) for name, start, stop in (
        ("left", 0, 100), ("middle", 20, 120),
        ("contained", 30, 90), ("right", 60, 140),
    )]
    if reverse:
        reads = [Read(read.name, reverse_complement(read.sequence)) for read in reads]
    return reference, reads


def trace_containment_branch_fixture(reverse=False):
    reference, reads = containment_branch_fixture(reverse)
    pool, diagnostics = assemble(reads)
    return {
        "reverse_input": reverse,
        "expected_maximal_span": reference[:140],
        "expected_maximal_span_length": 140,
        "output_lengths": sorted(
            len(record.current.sequence) for record in pool.active_derived
        ),
        "ambiguous_containments": diagnostics["ambiguous_containments"],
        "ambiguous_ends": diagnostics["ambiguous_ends"],
    }


def trace_case(seed, coverage=10):
    references, reads, origins = fixture(seed, "unique", 0.0, 0.0, coverage=coverage)
    sequences, groups, placements = _canonical_truth(reads, origins)
    nodes = [Read(f"node_{index}", sequence) for index, sequence in enumerate(sequences)]
    expected, boundary_containment = _truth_edges(placements)
    accepted, ambiguous, stages = _candidate_edges(nodes)
    containment_candidates, containments, repeated_parent = _containments(sequences)

    filtered = {}
    filtered_reasons = {}
    for key, edge in accepted.items():
        if edge.source[0] in containments or edge.target[0] in containments:
            filtered_reasons[_physical(*key)] = "containment_retired"
        elif any(node in repeated_parent and edge.overlap >= len(sequences[node])
                 for node in (edge.source[0], edge.target[0])):
            filtered_reasons[_physical(*key)] = "ambiguous_full_containment"
        else:
            filtered[key] = edge
    transitive, reduced, reciprocal, outgoing, incoming = _reduce(filtered)
    transitive_physical = {_physical(*key) for key in transitive}
    reciprocal_physical = {_physical(*key) for key in reciprocal}
    reduced_physical = {_physical(*key) for key in reduced}
    accepted_physical = {_physical(*key) for key in accepted}

    pool, graph = assemble(reads)
    raw_outputs = defaultdict(set)
    for output_index, record in enumerate(pool.active_derived):
        for placement in record.raw_placements or ():
            raw_outputs[placement.read_index].add(output_index)
    node_outputs = {node: set().union(*(raw_outputs[index] for index in raw_indices))
                    for node, raw_indices in enumerate(groups.values())}

    dispositions = Counter()
    separated_counts = Counter()
    actionable_examples = []
    branching_with_ambiguous_containment = 0
    for key, truths in expected.items():
        if key not in accepted_physical:
            disposition = stages.get(key, "no_anchor_vote")
        elif key in filtered_reasons:
            disposition = filtered_reasons[key]
        elif key in transitive_physical:
            disposition = "transitive_removed"
        elif key in reciprocal_physical:
            disposition = "reciprocal"
        elif key in reduced_physical:
            disposition = "branching_after_reduction"
        else:
            disposition = "unclassified"
        dispositions[disposition] += 1
        source, target = key
        if disposition == "branching_after_reduction" and (
            source[0] in containment_candidates and source[0] not in containments
            or target[0] in containment_candidates and target[0] not in containments
        ):
            branching_with_ambiguous_containment += 1
        separated = node_outputs[source[0]].isdisjoint(node_outputs[target[0]])
        if separated:
            separated_counts[disposition] += 1
            if disposition in ("branching_after_reduction", "no_anchor_vote"):
                actionable_examples.append({
                    "disposition": disposition, "truth": truths[0],
                    "source_outputs": sorted(node_outputs[source[0]]),
                    "target_outputs": sorted(node_outputs[target[0]]),
                    "source_out_degree": len(outgoing[source]),
                    "target_in_degree": len(incoming[target]),
                    "competing_targets": [list(edge.target) for edge in outgoing[source]],
                    "source_sequence": sequences[source[0]],
                    "target_sequence": sequences[target[0]],
                })

    containment_classes = Counter()
    for child, options in containment_candidates.items():
        loci = {_compatible_containment(child, option, placements, sequences)
                for option in options}
        if None in loci:
            containment_classes["truth_incompatible_candidate"] += 1
        elif len(loci) == 1 and len(options) > 1:
            containment_classes["multiple_representations_one_truth_locus"] += 1
        elif len(loci) > 1:
            containment_classes["multiple_truth_loci"] += 1
        else:
            containment_classes["unique_truth_placement"] += 1

    selected_extras = accepted_physical - set(expected) - set(boundary_containment)
    all_contained = set(containment_candidates)
    counterfactual_edges = {
        key: edge for key, edge in accepted.items()
        if edge.source[0] not in all_contained and edge.target[0] not in all_contained
    }
    counterfactual, counterfactual_diagnostics = _project_master_edges(
        nodes, counterfactual_edges, name_prefix="truth_only",
        excluded_contigs=all_contained,
    )
    counterfactual_lengths = [len(read.sequence) for read in counterfactual]
    actionable_examples.sort(key=lambda item: (
        len(item["source_sequence"]) + len(item["target_sequence"]),
        item["truth"]["source_interval"], item["truth"]["target_interval"],
    ))
    result = {
        "seed": seed, "coordinate_system": "0-based half-open",
        "target_coverage": coverage,
        "reads": len(reads), "nodes": len(nodes),
        "expected_physical_dovetails": len(expected),
        "truth_boundary_containment_edges": len(boundary_containment),
        "accepted_physical_edges": len(accepted_physical),
        "selected_edges_outside_truth_dovetails": len(selected_extras),
        "ambiguous_physical_edges": len(ambiguous),
        "expected_edge_dispositions": dict(sorted(dispositions.items())),
        "separated_output_dovetails": dict(sorted(separated_counts.items())),
        "branching_edges_incident_to_ambiguous_containment":
            branching_with_ambiguous_containment,
        "containment_nodes": len(containment_candidates),
        "containment_classes": dict(sorted(containment_classes.items())),
        "unique_containments_retired": len(containments),
        "truth_only_all_contained_removed": {
            "contigs": len(counterfactual), "n50": _n50(counterfactual_lengths),
            "longest": max(counterfactual_lengths, default=0),
            "diagnostics": asdict(counterfactual_diagnostics),
        },
        "graph": graph, "broken_examples": actionable_examples[:10],
        "reference_length": len(references[0]),
    }
    if seed == 199001 and coverage == 20:
        result["phase_a_ablation"] = _phase_a_ablation(
            nodes, filtered, expected, sequences, set(containments)
        )
    return result


def main():
    started = time.perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seeds", nargs="+", type=int, default=[95031, 95032])
    parser.add_argument("--coverage", type=float, default=10)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists")
    cases = [trace_case(seed, args.coverage) for seed in args.seeds]
    args.output.mkdir(parents=True)
    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__), Path(__file__).with_name("variable_length_validation.py"),
               *sorted((root / "src/anvaya").glob("*.py"))]
    snapshot_root = args.output / "source"
    for path in sources:
        destination = snapshot_root / path.relative_to(root)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
    summary = {
        "status": "complete",
        "seed_role": "development", "settings": {"minimum_overlap": 30,
            "anchor_k": 15, "anchors_per_read": 8, "occurrence_cap": 100},
        "cases": cases,
        "smallest_reproducer": [
            trace_containment_branch_fixture(False),
            trace_containment_branch_fixture(True),
        ],
        "command": "PYTHONPATH=src:experiments python experiments/layout_loss_trace.py "
                   "--output <new-directory> --seeds 199001 --coverage 20",
        "environment": {"python": sys.version, "platform": platform.platform()},
        "wall_seconds": time.perf_counter() - started,
        "maximum_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "source_sha256": {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in sources},
        "limitations": ["Clean unique synthetic fixtures only",
            "Truth classifies losses but is unavailable to assembly",
            "Two development seeds are correlated diagnostics, not validation",
            "Containment and edge counts are associations until a bounded ablation"],
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"cases": [{k: v for k, v in case.items()
                                 if k not in ("graph", "broken_examples")}
                                for case in cases]}, indent=2))


if __name__ == "__main__":
    main()
