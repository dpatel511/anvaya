"""Diagnose phase evidence in the consumed 95202 mosaic; never a validation run."""

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from controlled_assembly import fixture, score
from anvaya.damage_string_graph import assemble
from anvaya.phase_blocks import phase_component
from anvaya.raw_consensus import DamageProfile, project_raw_consensus


SEED = 95202
PANEL = "strains"
DEPTH = 3
DAMAGE = 0.4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    if arguments.output.exists():
        parser.error("output already exists")

    references, reads, origins, damage_events, _ = fixture(
        SEED, PANEL, DEPTH, DAMAGE,
    )
    rates = tuple(DAMAGE * 0.65 ** position for position in range(5))
    profile = DamageProfile(rates, rates)
    pool, graph = assemble(reads, damage_profile=profile)
    polished, consensus = project_raw_consensus(pool, profile)
    metrics = score(
        pool, references, origins, polished, damage_events=damage_events,
    )

    targets = []
    active = pool.active_derived
    for index, contig in enumerate(metrics["contig_evidence"]):
        if contig["sequence_category"] not in {
            "diagnostic_mosaic", "incompatible_candidate_sequences",
        }:
            continue
        record = active[index]
        evidence = phase_component(pool, record, profile)
        targets.append({
            "contig_index_1based": index + 1,
            "record_name": record.current.name,
            "sequence": polished[index].sequence,
            "evaluation": contig,
            "raw_placements": [asdict(item) for item in evidence.placements],
            "phase": asdict(evidence.phase),
            "phase_diagnostics": asdict(evidence.diagnostics),
        })

    root = Path(__file__).resolve().parents[1]
    sources = [
        Path(__file__),
        Path(__file__).with_name("controlled_assembly.py"),
        *sorted((root / "src/anvaya").glob("*.py")),
    ]
    arguments.output.mkdir(parents=True)
    (arguments.output / "summary.json").write_text(json.dumps({
        "status": "completed_diagnosis_not_validation",
        "seed_role": "consumed_validation_diagnosis",
        "settings": {
            "seed": SEED, "panel": PANEL, "depth": DEPTH, "damage": DAMAGE,
            "profile": {"five_prime_ct": rates, "three_prime_ga": rates},
        },
        "reads": len(reads),
        "contigs": len(active),
        "graph": graph,
        "consensus_diagnostics": asdict(consensus),
        "diagnostic_mosaic_contigs": metrics["diagnostic_mosaic_contigs"],
        "diagnostic_mosaic_bases": metrics["diagnostic_mosaic_bases"],
        "targets": targets,
        "source_sha256": {
            str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sources
        },
        "limitations": [
            "The seed was already consumed and this run is diagnosis, not validation",
            "Truth selects target contigs only after truth-blind assembly and consensus",
            "Final path placements omit alternative graph nodes rejected during layout",
            "Absence of variable final-path columns does not prove absence in whole input",
        ],
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
