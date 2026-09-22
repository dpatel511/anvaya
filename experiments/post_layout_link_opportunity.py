"""Bounded exact-link opportunity check on damage-string-graph output.

Synthetic origins are generated for later evaluation but never passed to assembly
or link discovery. Coordinates are 0-based and intervals are half-open.
"""
import argparse
import hashlib
import inspect
import json
import platform
import sys
from dataclasses import asdict
from pathlib import Path
from time import perf_counter

from anvaya.damage_string_graph import assemble
from anvaya.overlap_assembly import _n50
from anvaya.overlap_progressive_links import audit_raw_supported_progressive_links
from anvaya.raw_consensus import DamageProfile

from variable_length_validation import SCHEMA_VERSION, _case_identities, fixture


CASES = (
    ("clean_unique", 0.0, 20.0),
    ("damage_unique", 0.25, 10.0),
)
SEED = 199001


def _length_summary(lengths):
    ordered = sorted(lengths)
    return {
        "count": len(ordered),
        "bases": sum(ordered),
        "minimum": min(ordered, default=0),
        "median": ordered[len(ordered) // 2] if ordered else 0,
        "maximum": max(ordered, default=0),
        "n50": _n50(ordered),
    }


def evaluate_case(name, damage, coverage):
    references, reads, origins, events, errors = fixture(
        SEED, "unique", damage, 0.0, include_events=True, coverage=coverage,
    )
    rates = tuple(damage * 0.65 ** position for position in range(5))
    profile = DamageProfile(rates, rates)
    started = perf_counter()
    pool, graph = assemble(
        reads, damage_profile=profile if damage else None,
    )
    assembly_seconds = perf_counter() - started
    input_lengths = [len(record.current.sequence) for record in pool.active_derived]
    started = perf_counter()
    projection, links = audit_raw_supported_progressive_links(
        pool,
        allow_near_exact=False,
        minimum_read_support=2,
    )
    link_seconds = perf_counter() - started
    output_lengths = [len(read.sequence) for read in projection]
    return {
        "seed": SEED,
        "scenario": name,
        "target_coverage": coverage,
        "damage": damage,
        "reads": len(reads),
        **_case_identities(references, reads, origins, events, errors, profile),
        "input_contigs": _length_summary(input_lengths),
        "projected_contigs": _length_summary(output_lengths),
        "n50_ratio": (
            _n50(output_lengths) / _n50(input_lengths) if _n50(input_lengths) else None
        ),
        "changed_contig_count": len(input_lengths) - len(output_lengths),
        "graph": graph,
        "links": asdict(links),
        "coordinate_safety": "not_evaluable_without_composed_raw_placements",
        "seconds": {"assembly": assembly_seconds, "link_audit": link_seconds},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    if arguments.output.exists():
        parser.error("output already exists")
    arguments.output.mkdir(parents=True)
    root = Path(__file__).resolve().parents[1]
    source_paths = sorted({
        Path(__file__).resolve(),
        Path(inspect.getfile(fixture)).resolve(),
        *Path(inspect.getfile(assemble)).resolve().parent.glob("*.py"),
    })
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "hypothesis": (
            "baseline contigs retain exact dovetails bridged by at least two "
            "unambiguous raw molecules"
        ),
        "rejection_condition": (
            "no accepted useful join, clean projected N50 below 1.5x baseline, "
            "or any later coordinate-safety failure"
        ),
        "coordinate_system": "0-based; intervals half-open",
        "seed": SEED,
        "cases": [
            {"scenario": name, "damage": damage, "target_coverage": coverage}
            for name, damage, coverage in CASES
        ],
        "link_settings": {
            "allow_near_exact": False,
            "minimum_read_support": 2,
            "anchor_k": 15,
            "anchors_per_read": 8,
            "maximum_anchor_occurrences": 100,
            "minimum_anchor_matches": 2,
            "minimum_overlap": 30,
            "minimum_identity": 0.90,
            "minimum_ry_identity": 0.99,
        },
        "source_sha256": {
            path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in source_paths
        },
        "python": sys.version,
        "platform": platform.platform(),
    }
    (arguments.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8",
    )
    cases = [evaluate_case(*case) for case in CASES]
    result = {"manifest": manifest, "cases": cases}
    (arguments.output / "opportunity.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8",
    )
    for case in cases:
        print(json.dumps({
            "scenario": case["scenario"],
            "reads": case["reads"],
            "input_n50": case["input_contigs"]["n50"],
            "projected_n50": case["projected_contigs"]["n50"],
            "exact_dovetails": case["links"]["exact_dovetails"],
            "geometrically_spannable": case["links"]["geometrically_spannable_dovetails"],
            "supported_dovetails": case["links"]["supported_dovetails"],
            "accepted_paths": case["links"]["linear_paths"],
        }), flush=True)


if __name__ == "__main__":
    main()
