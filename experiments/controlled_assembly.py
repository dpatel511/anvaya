"""Small end-to-end progressive assembly mechanism benchmark with private truth.

Coordinates are 0-based half-open. No truth is passed to assembly or consensus.
"""
import argparse
import hashlib
import itertools
import json
import random
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from anvaya.overlap_assembly import _candidate_alignments, _n50
from anvaya.overlap_index import _anchor_index
from anvaya.overlap_progressive import (ProgressiveSequencePool,
    discover_progressive_raw_clusters, extend_progressive_raw_clusters,
    iterate_progressive_raw_extension)
from anvaya.raw_consensus import DamageProfile, project_raw_consensus
from anvaya.reads import Read
from anvaya.sequences import reverse_complement


def fixture(seed, panel, depth, damage, error_rate=0.0):
    if not 0.0 <= error_rate <= 1.0:
        raise ValueError("sequencing error rate must be between 0 and 1")
    rng = random.Random(seed)
    genome = "".join(rng.choice("ACGT") for _ in range(1200))
    if panel == "repeat":
        genome = genome[:700] + genome[200:350] + genome[850:]
    references = [genome]
    if panel == "strains":
        alternate = list(genome)
        for p in range(20, len(genome), 50):
            alternate[p] = next(b for b in "ACGT" if b != alternate[p])
        references.append("".join(alternate))
    # Separate streams ensure clean/damaged runs have identical molecules/origins.
    placement_rng = random.Random(seed + 1)
    damage_rng = random.Random(seed + 2)
    error_rng = random.Random(seed + 3)
    reads, origins, events, errors = [], [], [], []
    for ref, sequence in enumerate(references):
        for _ in range(depth * len(sequence) // 80):
            start = placement_rng.randrange(len(sequence) - 80 + 1)
            reverse = bool(placement_rng.randrange(2))
            raw = sequence[start:start + 80]
            if reverse:
                raw = reverse_complement(raw)
            observed = list(raw)
            for p, base in enumerate(raw):
                distance = p if base == "C" else len(raw) - 1 - p
                draw = damage_rng.random()
                if base in "CG" and distance < 5 and draw < damage * .65 ** distance:
                    observed[p] = "T" if base == "C" else "A"
                    events.append((len(reads), p, base, observed[p]))
            for p, base in enumerate(tuple(observed)):
                if error_rng.random() < error_rate:
                    observed[p] = error_rng.choice([candidate for candidate in "ACGT" if candidate != base])
                    errors.append((len(reads), p, base, observed[p]))
            reads.append(Read(f"molecule_{len(reads)}", "".join(observed), (35,) * 80))
            origins.append((ref, start, start + 80, reverse))
    return references, reads, origins, events, errors


def score(pool, references, origins, sequences=None, damage_events=(),
          reference_ambiguous_intervals=(), include_coordinate_states=False):
    """Score output with explicit 0-based origin and sequence compatibility evidence."""
    lengths, coverage, sequence_coverage, contig_evidence = [], set(), set(), []
    coordinate_observations = {}
    counts = Counter()
    ambiguous_intervals = {}
    for reference, start, stop in reference_ambiguous_intervals:
        if not 0 <= start <= stop <= len(references[reference]):
            raise ValueError("reference-ambiguous interval is outside its reference")
        ambiguous_intervals.setdefault(reference, []).append((start, stop))
    damage_by_read = {}
    for read_index, raw_position, latent, observed in damage_events:
        damage_by_read.setdefault(read_index, []).append(
            (raw_position, latent, observed)
        )
    if sequences is not None and len(sequences) != len(pool.active_derived):
        raise ValueError("consensus membership mismatch")
    for index, record in enumerate(pool.active_derived):
        sequence = record.current.sequence if sequences is None else sequences[index].sequence
        lengths.append(len(sequence))
        transforms = set()
        damage_support = {}
        for placement in record.raw_placements or ():
            ref, start, end, raw_reverse = origins[placement.read_index]
            reverse = raw_reverse != placement.reverse
            offset = (len(references[ref]) - end if reverse else start) - placement.offset
            transforms.add((ref, reverse, offset))
            read_length = len(pool.records[placement.read_index].raw.sequence)
            for raw_position, latent, observed in damage_by_read.get(placement.read_index, ()):
                if not placement.read_start <= raw_position < placement.read_stop:
                    continue
                contig_position = placement.offset + (
                    read_length - 1 - raw_position if placement.reverse else raw_position
                )
                if not 0 <= contig_position < len(sequence):
                    continue
                if placement.reverse:
                    latent = reverse_complement(latent)
                    observed = reverse_complement(observed)
                damage_support.setdefault(contig_position, set()).add((latent, observed))
        evidence = dict(
            contig_id=f"unitig_{index + 1}",
            record_name=record.current.name,
            length=len(sequence),
            transforms=[],
            diagnostic_positions_0based=[],
            diagnostic_alleles=[],
        )
        if not transforms:
            counts["missing_provenance_contigs"] += 1
            counts["missing_provenance_bases"] += len(sequence)
            evidence.update(
                provenance_category="missing_provenance",
                sequence_category="not_evaluable",
                compatible_transform_count=0,
            )
            counts["not_evaluable_contigs"] += 1
            counts["not_evaluable_bases"] += len(sequence)
            contig_evidence.append(evidence)
            continue

        evidence["provenance_category"] = (
            "unique_transform" if len(transforms) == 1 else
            "same_reference_multiple_transforms"
            if len({transform[0] for transform in transforms}) == 1 else
            "cross_reference_transforms"
        )
        counts[f"{evidence['provenance_category']}_contigs"] += 1
        if len(transforms) > 1:
            counts["conflicting_origin_contigs"] += 1

        candidates = []
        for ref, reverse, start in sorted(transforms):
            oriented = reverse_complement(references[ref]) if reverse else references[ref]
            in_bounds = 0 <= start and start + len(sequence) <= len(oriented)
            candidate = oriented[start:start + len(sequence)] if in_bounds else None
            mismatch_positions = (
                [position for position, (a, b) in enumerate(zip(sequence, candidate)) if a != b]
                if candidate is not None else None
            )
            mismatches = len(mismatch_positions) if mismatch_positions is not None else None
            damage_supported = (
                [position for position in mismatch_positions
                 if (candidate[position], sequence[position]) in damage_support.get(position, ())]
                if mismatch_positions is not None else None
            )
            unexplained = (
                [position for position in mismatch_positions if position not in damage_supported]
                if mismatch_positions is not None else None
            )
            evidence["transforms"].append(dict(
                reference=ref,
                reverse=reverse,
                start_0based=start,
                stop_exclusive=start + len(sequence),
                in_bounds=in_bounds,
                mismatches=mismatches,
                damage_supported_mismatches=(len(damage_supported) if damage_supported is not None else None),
                unexplained_mismatch_positions_0based=unexplained,
                damage_compatible=(not unexplained if unexplained is not None else False),
            ))
            if candidate is not None:
                candidates.append((ref, reverse, start, candidate, mismatches, unexplained))

        if not candidates:
            counts["out_of_bounds_contigs"] += 1
            counts["out_of_bounds_bases"] += len(sequence)
            evidence.update(sequence_category="not_evaluable", compatible_transform_count=0)
            counts["not_evaluable_contigs"] += 1
            counts["not_evaluable_bases"] += len(sequence)
            contig_evidence.append(evidence)
            continue

        expected_sequences = {candidate[3] for candidate in candidates}
        exact_candidates = [candidate for candidate in candidates if candidate[4] == 0]
        damage_candidates = [candidate for candidate in candidates if not candidate[5]]
        mismatch_counts = [candidate[4] for candidate in candidates]
        diagnostic_positions = [
            position for position in range(len(sequence))
            if len({candidate[3][position] for candidate in candidates}) > 1
        ]
        evidence["diagnostic_positions_0based"] = diagnostic_positions
        evidence["diagnostic_alleles"] = [
            dict(position_0based=position,
                 alleles=sorted({candidate[3][position] for candidate in candidates}))
            for position in diagnostic_positions
        ]
        evidence["compatible_transform_count"] = len(exact_candidates)
        evidence["damage_compatible_transform_count"] = len(damage_candidates)
        evidence["candidate_mismatch_lower_bound"] = min(mismatch_counts)
        evidence["candidate_mismatch_upper_bound"] = max(mismatch_counts)
        counts["candidate_evaluable_bases"] += len(sequence)
        counts["candidate_mismatch_lower_bound"] += min(mismatch_counts)
        counts["candidate_mismatch_upper_bound"] += max(mismatch_counts)

        evaluated_mismatches = None
        coordinate_candidate = candidates[0] if len(transforms) == 1 else None
        sequence_resolved_candidate = None
        if len(transforms) == 1:
            ref, reverse, start, _, mismatches, _ = candidates[0]
            evidence["sequence_category"] = "unique_origin_sequence"
            counts["resolved_contigs"] += 1
            counts["resolved_bases"] += len(sequence)
            counts["mismatches"] += mismatches
            evaluated_mismatches = mismatches
            oriented_length = len(references[ref])
            for position in range(start, start + len(sequence)):
                coordinate = oriented_length - 1 - position if reverse else position
                coverage.add((ref, coordinate))
                sequence_coverage.add((ref, coordinate))
        elif len(expected_sequences) == 1:
            evidence["sequence_category"] = "identical_candidate_sequences"
            evaluated_mismatches = candidates[0][4]
        elif exact_candidates:
            evidence["sequence_category"] = "matches_candidate_sequence"
            evaluated_mismatches = 0
            if len(exact_candidates) == 1:
                sequence_resolved_candidate = exact_candidates[0]
        elif damage_candidates:
            damage_mismatch_counts = {candidate[4] for candidate in damage_candidates}
            if len(damage_mismatch_counts) == 1:
                evidence["sequence_category"] = "matches_damage_event_candidate"
                evaluated_mismatches = next(iter(damage_mismatch_counts))
                if len(damage_candidates) == 1:
                    sequence_resolved_candidate = damage_candidates[0]
            else:
                evidence["sequence_category"] = "damage_compatible_candidates_disagree"
        else:
            matched_references = set()
            all_diagnostic_bases_explained = True
            for position in diagnostic_positions:
                matching = {
                    candidate[0] for candidate in candidates
                    if candidate[3][position] == sequence[position]
                }
                if not matching:
                    all_diagnostic_bases_explained = False
                    break
                matched_references.update(matching)
            evidence["sequence_category"] = (
                "diagnostic_mosaic"
                if diagnostic_positions
                and all_diagnostic_bases_explained
                and len(matched_references) > 1
                else "incompatible_candidate_sequences"
            )

        if sequence_resolved_candidate is not None:
            ref, reverse, start = sequence_resolved_candidate[:3]
            coordinate_candidate = sequence_resolved_candidate
            oriented_length = len(references[ref])
            for position in range(start, start + len(sequence)):
                coordinate = oriented_length - 1 - position if reverse else position
                sequence_coverage.add((ref, coordinate))
            evidence["sequence_resolved_transform"] = dict(
                reference=ref,
                reverse=reverse,
                start_0based=start,
                stop_exclusive=start + len(sequence),
            )
            counts["sequence_resolved_ambiguous_origin_contigs"] += 1

        if coordinate_candidate is not None:
            ref, reverse, start = coordinate_candidate[:3]
            reference_length = len(references[ref])
            for position, base in enumerate(sequence):
                oriented_position = start + position
                coordinate = (
                    reference_length - 1 - oriented_position
                    if reverse else oriented_position
                )
                reference_oriented_base = reverse_complement(base) if reverse else base
                coordinate_observations.setdefault((ref, coordinate), []).append(
                    (reference_oriented_base, evidence["contig_id"],
                     "unique_origin" if len(transforms) == 1 else "sequence_resolved")
                )

        counts[f"{evidence['sequence_category']}_contigs"] += 1
        counts[f"{evidence['sequence_category']}_bases"] += len(sequence)
        if evaluated_mismatches is None:
            counts["accuracy_unresolved_bases"] += len(sequence)
        else:
            counts["accuracy_evaluated_contigs"] += 1
            counts["accuracy_evaluated_bases"] += len(sequence)
            counts["accuracy_mismatches"] += evaluated_mismatches
        contig_evidence.append(evidence)

    total_bases = sum(lengths)
    counts["nonduplicated_accuracy_projected_bases"] = sum(
        len(observations) for observations in coordinate_observations.values()
    )
    counts["nonduplicated_accuracy_duplicate_observations"] = sum(
        len(observations) - 1 for observations in coordinate_observations.values()
    )
    discordant_sites = []
    coordinate_states = []
    for (reference, coordinate), observations in coordinate_observations.items():
        observed_bases = {base for base, _, _ in observations}
        state = dict(
            reference=reference,
            position_0based=coordinate,
            reference_base=references[reference][coordinate],
            observed_bases=sorted(observed_bases),
            observation_count=len(observations),
            origin_resolution=sorted({resolution for _, _, resolution in observations}),
        )
        if len(observed_bases) != 1:
            counts["nonduplicated_accuracy_discordant_reference_bases"] += 1
            state["status"] = "discordant"
            discordant_sites.append({
                "reference": reference,
                "position_0based": coordinate,
                "observations": [
                    {"base": base, "contig_id": contig_id}
                    for base, contig_id, _ in observations
                ],
            })
            coordinate_states.append(state)
            continue
        counts["nonduplicated_accuracy_evaluated_reference_bases"] += 1
        if next(iter(observed_bases)) != references[reference][coordinate]:
            counts["nonduplicated_accuracy_mismatches"] += 1
            state["status"] = "concordant_incorrect"
        else:
            state["status"] = "concordant_correct"
        coordinate_states.append(state)
    counts["unresolved_bases"] = total_bases - counts["resolved_bases"]
    counts["accuracy_unresolved_bases"] += (
        total_bases
        - counts["accuracy_evaluated_bases"]
        - counts["accuracy_unresolved_bases"]
    )
    for metric in (
        "resolved_contigs", "resolved_bases", "mismatches",
        "accuracy_evaluated_contigs", "accuracy_evaluated_bases",
        "accuracy_mismatches", "accuracy_unresolved_bases",
        "nonduplicated_accuracy_projected_bases",
        "nonduplicated_accuracy_duplicate_observations",
        "nonduplicated_accuracy_evaluated_reference_bases",
        "nonduplicated_accuracy_mismatches",
        "nonduplicated_accuracy_discordant_reference_bases",
        "diagnostic_mosaic_contigs", "diagnostic_mosaic_bases",
        "incompatible_candidate_sequences_contigs",
        "incompatible_candidate_sequences_bases",
        "matches_damage_event_candidate_contigs",
        "matches_damage_event_candidate_bases",
        "damage_compatible_candidates_disagree_contigs",
        "damage_compatible_candidates_disagree_bases",
        "not_evaluable_contigs", "not_evaluable_bases",
        "candidate_evaluable_bases", "candidate_mismatch_lower_bound",
        "candidate_mismatch_upper_bound",
    ):
        counts.setdefault(metric, 0)
    sequence_ambiguous_coverage = {
        (reference, position)
        for reference, position in sequence_coverage
        if any(start <= position < stop
               for start, stop in ambiguous_intervals.get(reference, ()))
    }
    result = dict(counts, contigs=len(lengths), bases=total_bases, n50=_n50(lengths),
                longest=max(lengths, default=0), resolved_unique_reference_bases=len(coverage),
                sequence_resolved_unique_reference_bases=len(sequence_coverage),
                sequence_resolved_nonambiguous_reference_bases=(
                    len(sequence_coverage) - len(sequence_ambiguous_coverage)
                ),
                sequence_resolved_ambiguous_reference_bases=len(sequence_ambiguous_coverage),
                nonduplicated_accuracy_discordant_sites=discordant_sites,
                contig_evidence=contig_evidence)
    if include_coordinate_states:
        result["nonduplicated_accuracy_coordinate_states"] = sorted(
            coordinate_states,
            key=lambda state: (state["reference"], state["position_0based"]),
        )
    return result


def run_case(seed, panel, depth, damage, output, reuse_raw_evidence=False, constrained=False,
             string_layout=False, damage_string_layout=False, error_rate=0.0):
    references, reads, origins, events, errors = fixture(
        seed, panel, depth, damage, error_rate
    )
    reference_ambiguous_intervals = (
        [(0, 200, 350), (0, 700, 850)] if panel == "repeat" else []
    )
    output.mkdir(parents=True, exist_ok=False)
    (output / "truth.json").write_text(json.dumps(dict(
        references=references, origins=origins, damage_events=events,
        sequencing_errors=errors,
        reference_ambiguous_intervals=reference_ambiguous_intervals,
    )))
    anchors = _anchor_index(reads, 15, 0, 8, 100)
    ids = list(range(len(reads)))
    stages = Counter()
    for index, (ref, start, end, reverse) in enumerate(origins):
        trace = {}
        _candidate_alignments(reads[index].sequence, reads, ids, anchors, {index},
            anchor_k=15, anchors_per_read=8, maximum_anchor_occurrences=100,
            minimum_anchor_matches=1, minimum_overlap=30, minimum_identity=.9,
            minimum_ry_identity=.99, position_bits=7, target_window=80, stage_trace=trace)
        for other, (r, s, e, rev) in enumerate(origins):
            if r == ref and 0 < abs(s - start) <= 50:
                offset = end - e if reverse else s - start
                stages[trace.get((other, reverse != rev, offset), "no_anchor_vote")] += 1
    pool = ProgressiveSequencePool.from_reads(reads)
    clusters, clustering = discover_progressive_raw_clusters(pool, minimum_anchor_matches=1)
    owners = {member: number for number, cluster in enumerate(clusters) for member in cluster.member_indices}
    ownership = Counter()
    for i, (ref, start, end, reverse) in enumerate(origins):
        for j in range(i + 1, len(origins)):
            r, s, e, rev = origins[j]
            if ref == r and 0 < abs(start - s) <= 50:
                label = "unassigned_endpoint" if i not in owners or j not in owners else (
                    "same_cluster" if owners[i] == owners[j] else "different_clusters")
                ownership[label] += 1
    components = []
    for ref in range(len(references)):
        begin = stop = None
        for start, end in sorted((s, e) for r, s, e, rev in origins if r == ref):
            if stop is None or stop - start < 30:
                if stop is not None:
                    components.append(stop - begin)
                begin, stop = start, end
            else:
                stop = max(stop, end)
        if stop is not None:
            components.append(stop - begin)
    pool, extension = extend_progressive_raw_clusters(pool, clusters, minimum_anchor_matches=1,
                                                     track_raw_placements=True)
    pool, iteration = iterate_progressive_raw_extension(pool, minimum_anchor_matches=1,
                                                        audit_rejected_extensions=True,
                                                        reuse_raw_evidence=reuse_raw_evidence or constrained,
                                                        minimum_identity=1.0 if constrained else .9)
    retired = 0
    graph = None
    rates = tuple(damage * .65 ** p for p in range(5))
    if string_layout or damage_string_layout:
        from string_layout import assemble
        profile = DamageProfile(rates, rates) if damage_string_layout else None
        pool, graph = assemble(reads, damage_profile=profile)
    unconsolidated = score(
        pool, references, origins, damage_events=events,
        reference_ambiguous_intervals=reference_ambiguous_intervals,
    )
    if constrained:
        from compatible_layout import consolidate
        pool, retired = consolidate(pool)
    before = score(
        pool, references, origins, damage_events=events,
        reference_ambiguous_intervals=reference_ambiguous_intervals,
    )
    polished, consensus = project_raw_consensus(pool, DamageProfile(rates, rates))
    result = dict(seed=seed, panel=panel, depth=depth, damage=damage,
                  error_rate=error_rate, reads=len(reads),
                  reuse_raw_evidence=reuse_raw_evidence,
                  constrained=constrained, consolidated_contigs=retired, unconsolidated=unconsolidated,
                  string_layout=string_layout, damage_string_layout=damage_string_layout, graph=graph,
                  truth_overlap_components=components, truth_pair_ownership=dict(ownership),
                  expected_directed_overlap_stages=dict(stages), clustering=asdict(clustering),
                  extension=asdict(extension), iteration=asdict(iteration),
                  before=before,
                  after=score(
                      pool, references, origins, polished, damage_events=events,
                      reference_ambiguous_intervals=reference_ambiguous_intervals,
                  ),
                  consensus=asdict(consensus))
    (output / "summary.json").write_text(json.dumps(result, indent=2))
    return result


def run(output, reuse_raw_evidence=False, constrained=False, string_layout=False,
        damage_string_layout=False, seeds=(93001, 93002), error_rate=0.0):
    output.mkdir(parents=True, exist_ok=False)
    results = []
    for seed, panel, depth, damage in itertools.product(seeds,
            ("unique", "repeat", "strains"), (3, 10, 20), (0., .4)):
        result = run_case(seed, panel, depth, damage,
                          output / f"{seed}-{panel}-{depth}-{damage}", reuse_raw_evidence,
                          constrained, string_layout, damage_string_layout, error_rate)
        results.append(result)
        progress = {k: result[k] for k in ("seed", "panel", "depth", "damage")}
        progress["before"] = {
            key: value for key, value in result["before"].items()
            if key != "contig_evidence"
        }
        progress["after"] = {
            key: value for key, value in result["after"].items()
            if key != "contig_evidence"
        }
        print(json.dumps(progress), flush=True)
    root = Path(__file__).resolve().parents[1]
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in [Path(__file__), Path(__file__).with_name("compatible_layout.py"), Path(__file__).with_name("string_layout.py"),
                        *sorted((root / "src/anvaya").glob("*.py"))]}
    (output / "summary.json").write_text(json.dumps(dict(cases=results, source_sha256=hashes,
        error_rate=error_rate,
        limitations=["Small synthetic mechanism benchmark, not a realistic metagenome",
                     "80bp fragments, equal strain abundance, two reference seeds, no indels/PCR duplicates",
                     "Independent substitution errors are optional; Q35 remains a nominal consensus input",
                     "Matched generating damage profile is optimistic",
                     "Progressive clustering/extension and one consensus pass; optional linking/rescues excluded",
                     "Conflicting origins may reflect repeats or strains, not necessarily sequence errors",
                     "Origin-consistent contigs only contribute accuracy and reference recovery metrics"]), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reuse-raw-evidence", action="store_true")
    parser.add_argument("--constrained", action="store_true", help="exact reuse and placement-compatible containment")
    parser.add_argument("--string-layout", action="store_true", help="exact read-overlap path control")
    parser.add_argument("--damage-string-layout", action="store_true",
                        help="terminal-damage-compatible read-overlap path experiment")
    parser.add_argument("--seeds", type=int, nargs="+", default=(93001, 93002),
                        help="simulation seeds (default: 93001 93002)")
    parser.add_argument("--error-rate", type=float, default=0.0,
                        help="independent substitution-error probability (default: 0)")
    args = parser.parse_args()
    if (args.string_layout or args.damage_string_layout) and (args.constrained or args.reuse_raw_evidence):
        parser.error("string layout is a separate experiment")
    if args.string_layout and args.damage_string_layout:
        parser.error("choose exact or damage-aware string layout")
    run(args.output, args.reuse_raw_evidence, args.constrained, args.string_layout,
        args.damage_string_layout, tuple(args.seeds), args.error_rate)
