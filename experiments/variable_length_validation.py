"""Development variable-length validation for the damage-aware string graph.

Coordinates and origins are 0-based, half-open. Truth is used only by scoring.
"""
import argparse
import hashlib
import inspect
import json
import random
from dataclasses import asdict
from pathlib import Path
from time import perf_counter

from anvaya.damage_string_graph import assemble
from anvaya.exact_string_graph import (
    assemble_damage_aware_unitigs,
    assemble_exact_unitigs,
    evidence_pool,
)
from anvaya.raw_consensus import DamageProfile, project_raw_consensus
from anvaya.reads import Read
from anvaya.sequences import reverse_complement

from controlled_assembly import score


SCENARIOS = (
    ("clean_unique", "unique", 0.0, 0.0, 1.0),
    ("damage_unique", "unique", 0.25, 0.0, 1.0),
    ("damage_error_unique", "unique", 0.25, 0.002, 1.0),
    ("profile_mismatch_unique", "unique", 0.25, 0.0, 0.5),
    ("damage_error_repeat", "repeat", 0.25, 0.002, 1.0),
    ("damage_error_strains_80_20", "strains", 0.25, 0.002, 1.0),
)
ERROR_ONLY_SCENARIOS = (
    ("error_only_unique", "unique", 0.0, 0.002, 1.0),
    ("error_only_repeat", "repeat", 0.0, 0.002, 1.0),
    ("error_only_strains_80_20", "strains", 0.0, 0.002, 1.0),
)
DEVELOPMENT_SEEDS = (199001, 199002)
HELD_OUT_SEEDS = (95201, 95202, 95203, 95204)
CONSUMED_VALIDATION_SEEDS = frozenset((
    94011, 94012, *HELD_OUT_SEEDS, 96011, 96012,
    197101, 197102, 197103, 197104,
    199001, 199002,
    20261401, 20261402,
))
REGISTERED_VALIDATION_SEEDS = frozenset()
SCHEMA_VERSION = 2


def _sha256_json(value):
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _case_identities(references, reads, origins, events, errors, profile):
    """Hash actual ordered evaluator inputs with an unambiguous JSON encoding."""
    return {
        "input_sha256": _sha256_json([
            {
                "read_id": read.name,
                "molecule_id": read.name,
                "sequence": read.sequence,
                "qualities": read.qualities,
            }
            for read in reads
        ]),
        "reference_sha256": _sha256_json(references),
        "truth_sha256": _sha256_json({
            "origins": origins,
            "damage_events": events,
            "sequencing_errors": errors,
        }),
        "damage_profile_sha256": _sha256_json({
            "five_prime_ct": profile.five_prime_ct,
            "three_prime_ga": profile.three_prime_ga,
        }),
    }


def fixture(seed, panel, damage, error_rate, include_events=False, coverage=10):
    if coverage <= 0:
        raise ValueError("coverage must be positive")
    reference_rng = random.Random(seed)
    genome = "".join(reference_rng.choice("ACGT") for _ in range(2400))
    if panel == "repeat":
        genome = genome[:1400] + genome[400:700] + genome[1700:]
    references = [genome]
    depths = [coverage]
    if panel == "strains":
        alternate = list(genome)
        for position in range(30, len(genome), 75):
            alternate[position] = next(base for base in "ACGT" if base != alternate[position])
        references.append("".join(alternate))
        depths = [coverage * 0.8, coverage * 0.2]

    placement_rng = random.Random(seed + 1)
    damage_rng = random.Random(seed + 2)
    error_rng = random.Random(seed + 3)
    reads, origins, events, errors = [], [], [], []
    expected_length = (31 + 200) / 2
    for reference_index, (reference, depth) in enumerate(zip(references, depths)):
        for _ in range(round(depth * len(reference) / expected_length)):
            length = placement_rng.randint(31, 200)
            start = placement_rng.randrange(len(reference) - length + 1)
            reverse = bool(placement_rng.randrange(2))
            source = reference[start:start + length]
            if reverse:
                source = reverse_complement(source)
            observed = list(source)
            qualities = [35] * length
            for position, base in enumerate(source):
                distance = position if base == "C" else length - 1 - position
                if (
                    base in "CG" and distance < 5
                    and damage_rng.random() < damage * 0.65 ** distance
                ):
                    observed[position] = "T" if base == "C" else "A"
                    events.append((len(reads), position, base, observed[position]))
            for position, base in enumerate(observed):
                if error_rng.random() < error_rate:
                    observed[position] = error_rng.choice(
                        [candidate for candidate in "ACGT" if candidate != base]
                    )
                    qualities[position] = 10
                    errors.append((len(reads), position, base, observed[position]))
                    events = [
                        event for event in events
                        if event[:2] != (len(reads), position)
                    ]
            reads.append(Read(f"molecule_{len(reads)}", "".join(observed), tuple(qualities)))
            origins.append((reference_index, start, start + length, reverse))
    result = references, reads, origins
    return (*result, events, errors) if include_events else result


def evaluate(
    seed, name, panel, damage, error_rate, profile_scale, coverage=10,
    implementation="baseline",
):
    references, reads, origins, events, errors = fixture(
        seed, panel, damage, error_rate, include_events=True, coverage=coverage,
    )
    true_rates = tuple(damage * 0.65 ** position for position in range(5))
    supplied = tuple(rate * profile_scale for rate in true_rates)
    profile = DamageProfile(supplied, supplied)
    identities = _case_identities(
        references, reads, origins, events, errors, profile,
    )
    results = {}
    reference_ambiguous_intervals = (
        [(0, 400, 700), (0, 1400, 1700)] if panel == "repeat" else []
    )
    for mode, graph_profile in (("exact", None), ("damage_aware", profile)):
        started = perf_counter()
        if implementation == "baseline":
            pool, graph = assemble(reads, damage_profile=graph_profile)
        elif implementation == "candidate":
            if graph_profile is None:
                unitig_graph = assemble_exact_unitigs(reads)
                graph = {
                    "topology_nodes": len(unitig_graph.topology_sequences),
                    "candidate_directed_edges": len(unitig_graph.candidate_edges),
                    "reduced_directed_edges": len(unitig_graph.reduced_edges),
                    "unitigs": len(unitig_graph.unitigs),
                    "contained_groups": unitig_graph.contained_groups,
                    "ambiguous_evidence_reads": unitig_graph.ambiguous_evidence_reads,
                }
            else:
                unitig_graph, graph = assemble_damage_aware_unitigs(
                    reads, graph_profile,
                )
                graph.update({
                    "topology_nodes": len(unitig_graph.topology_sequences),
                    "candidate_directed_edges": len(unitig_graph.candidate_edges),
                    "reduced_directed_edges": len(unitig_graph.reduced_edges),
                    "unitigs": len(unitig_graph.unitigs),
                    "contained_groups": unitig_graph.contained_groups,
                    "ambiguous_evidence_reads": unitig_graph.ambiguous_evidence_reads,
                })
            pool, placement_diagnostics = evidence_pool(reads, unitig_graph)
            graph["evidence"] = placement_diagnostics
        else:
            raise ValueError(f"unknown implementation: {implementation}")
        assembly_seconds = perf_counter() - started
        polished, consensus = project_raw_consensus(pool, profile)
        results[mode] = {
            "layout": score(
                pool, references, origins, damage_events=events,
                reference_ambiguous_intervals=reference_ambiguous_intervals,
                include_coordinate_states=True,
            ),
            "consensus": score(
                pool, references, origins, polished, damage_events=events,
                reference_ambiguous_intervals=reference_ambiguous_intervals,
                include_coordinate_states=True,
            ),
            "graph": graph,
            "consensus_diagnostics": asdict(consensus),
            "assembly_seconds": assembly_seconds,
        }
    return {
        "seed": seed,
        "scenario": name,
        "panel": panel,
        "damage": damage,
        "sequencing_error_rate": error_rate,
        "supplied_profile_scale": profile_scale,
        "target_coverage": coverage,
        "reads": len(reads),
        "simulated_damage_events": len(events),
        "simulated_sequencing_errors": len(errors),
        "minimum_read_length": min(len(read.sequence) for read in reads),
        "maximum_read_length": max(len(read.sequence) for read in reads),
        **identities,
        "results": results,
    }


def experiment_manifest(seeds, seed_role, scenarios, coverages):
    root = Path(__file__).resolve().parents[1]
    evaluator_path = Path(inspect.getfile(score)).resolve()
    return {
        "schema_version": SCHEMA_VERSION,
        "seed_role": seed_role,
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "evaluator_sha256": hashlib.sha256(evaluator_path.read_bytes()).hexdigest(),
        "modes": ["exact", "damage_aware"],
        "stages": ["layout", "consensus"],
        "scenario_definitions": [list(scenario) for scenario in scenarios],
        "target_coverages": list(coverages),
        "expected_cases": [
            {"seed": seed, "scenario": scenario[0], "target_coverage": coverage}
            for seed in seeds for coverage in coverages for scenario in scenarios
        ],
        "generator_path": Path(__file__).resolve().relative_to(root).as_posix(),
        "evaluator_path": evaluator_path.relative_to(root).as_posix(),
    }


def resolve_seed_plan(seeds=None, held_out=False, error_only=False, seed_role=None):
    if held_out and seeds:
        raise ValueError("--held-out and --seeds are mutually exclusive")
    if held_out and seed_role:
        raise ValueError("--held-out already denotes consumed validation diagnosis")
    if seed_role == "validation" and not seeds:
        raise ValueError("fresh validation requires explicit --seeds")
    selected = HELD_OUT_SEEDS if held_out else (
        tuple(seeds) if seeds else DEVELOPMENT_SEEDS
    )
    role = "consumed_validation_diagnosis" if held_out else (seed_role or "development")
    consumed = sorted(set(selected) & CONSUMED_VALIDATION_SEEDS)
    if role == "validation" and consumed:
        raise ValueError(f"consumed validation seeds cannot be fresh validation: {consumed}")
    unregistered = sorted(set(selected) - REGISTERED_VALIDATION_SEEDS)
    if role == "validation" and unregistered:
        raise ValueError(f"fresh validation seeds are not preregistered: {unregistered}")
    if error_only:
        role += "_error_only"
    return selected, role


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seeds", nargs="+", type=int)
    parser.add_argument("--held-out", action="store_true")
    parser.add_argument("--error-only", action="store_true")
    parser.add_argument("--depths", nargs="+", type=float, default=(2.0, 5.0, 10.0, 20.0))
    parser.add_argument(
        "--implementation", choices=("baseline", "candidate"), required=True,
    )
    parser.add_argument(
        "--seed-role",
        choices=("development", "validation", "consumed_validation_diagnosis"),
    )
    arguments = parser.parse_args()
    try:
        seeds, seed_role = resolve_seed_plan(
            arguments.seeds, arguments.held_out, arguments.error_only,
            arguments.seed_role,
        )
    except ValueError as error:
        parser.error(str(error))
    if arguments.output.exists():
        parser.error("output already exists")
    if any(depth <= 0 for depth in arguments.depths):
        parser.error("depths must be positive")
    scenarios = ERROR_ONLY_SCENARIOS if arguments.error_only else SCENARIOS
    manifest = experiment_manifest(seeds, seed_role, scenarios, arguments.depths)
    arguments.output.mkdir(parents=True)
    (arguments.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8",
    )
    cases = []
    for seed in seeds:
        for depth in arguments.depths:
            for scenario in scenarios:
                case = evaluate(
                    seed, *scenario, coverage=depth,
                    implementation=arguments.implementation,
                )
                cases.append(case)
                exact = case["results"]["exact"]["consensus"]
                aware = case["results"]["damage_aware"]["consensus"]
                omitted = {
                    "contig_evidence", "nonduplicated_accuracy_coordinate_states",
                    "nonduplicated_accuracy_discordant_sites",
                }
                print(json.dumps({
                    "seed": seed,
                    "target_coverage": depth,
                    "scenario": scenario[0],
                    "exact": {key: value for key, value in exact.items() if key not in omitted},
                    "damage_aware": {
                        key: value for key, value in aware.items() if key not in omitted
                    },
                }), flush=True)
    root = Path(__file__).resolve().parents[1]
    imported_objects = {
        "assemble": assemble,
        "assemble_damage_aware_unitigs": assemble_damage_aware_unitigs,
        "assemble_exact_unitigs": assemble_exact_unitigs,
        "evidence_pool": evidence_pool,
        "project_raw_consensus": project_raw_consensus,
        "DamageProfile": DamageProfile,
        "Read": Read,
        "reverse_complement": reverse_complement,
        "score": score,
    }
    imported_module_paths = {
        name: str(Path(inspect.getfile(value)).resolve().relative_to(root))
        for name, value in imported_objects.items()
    }
    package_root = Path(inspect.getfile(assemble)).resolve().parent
    sources = sorted({
        Path(__file__).resolve(),
        Path(inspect.getfile(score)).resolve(),
        *package_root.glob("*.py"),
    })
    source_sha256 = {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sources
    }
    assembly_path = Path(inspect.getfile(
        assemble if arguments.implementation == "baseline"
        else assemble_damage_aware_unitigs
    )).resolve()
    assembly_source_sha256 = hashlib.sha256(assembly_path.read_bytes()).hexdigest()
    (arguments.output / "summary.json").write_text(
        json.dumps({
            "seeds": seeds,
            "seed_role": seed_role,
            "implementation": arguments.implementation,
            "consumed_validation_seeds": sorted(CONSUMED_VALIDATION_SEEDS),
            "registered_validation_seeds": sorted(REGISTERED_VALIDATION_SEEDS),
            "scenarios": [scenario[0] for scenario in scenarios],
            "target_coverages": arguments.depths,
            "experiment_manifest": manifest,
            "assembly_source_sha256": assembly_source_sha256,
            "cases": cases,
            "imported_module_paths": imported_module_paths,
            "source_sha256": source_sha256,
            "limitations": [
                "Small synthetic references with 31-200 bp fragments",
                "Substitution errors are independent and assigned Q10",
                "No indels, PCR duplicates, environmental contamination or abundance spectrum",
                "Repeat origin conflicts are unresolved evidence, not automatically false joins",
            ],
        }, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
