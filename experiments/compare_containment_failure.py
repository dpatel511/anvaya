"""Join two containment failure traces by immutable raw-read membership."""
import argparse
import json
from pathlib import Path


def compare(before, after):
    if before["target"] != after["target"]:
        raise ValueError("trace targets differ")
    if before["input_sha256"] != after["input_sha256"]:
        raise ValueError("trace inputs differ")
    stages = {}
    for stage in ("layout", "consensus"):
        old_outputs = before[stage]["outputs"]
        old_catalog = before[stage]["output_catalog"]
        old_target_by_record = {output["record_name"]: output for output in old_outputs}
        rows = []
        for new in after[stage]["outputs"]:
            new_names = set(new["raw_read_names"])
            candidates = []
            for old in old_catalog:
                old_names = set(old["raw_read_names"])
                union = old_names | new_names
                overlap = old_names & new_names
                candidates.append((len(overlap) / len(union) if union else 1.0, old))
            similarity, old = max(candidates, key=lambda item: item[0], default=(0.0, None))
            old_target = old_target_by_record.get(old["record_name"]) if old else None
            if old is not None and old_target is None:
                old_names = set(old["raw_read_names"])
                old_target = next(
                    (output for output in old_outputs
                     if set(output["raw_read_names"]) == old_names),
                    None,
                )
            old_candidate_bases = (
                {item["reference_oriented_base"]
                 for item in old["target_candidate_projections"]}
                if old else set()
            )
            if old is None or similarity == 0:
                classification = "new_provenance"
            elif old_target is None and new["reference_oriented_base"] in old_candidate_bases:
                classification = (
                    "newly_resolved_existing_spelling_same_raw_membership"
                    if old["raw_read_names"] == new["raw_read_names"] else
                    "newly_resolved_existing_spelling_membership_changed"
                )
            elif old_target is None and old["raw_read_names"] == new["raw_read_names"]:
                classification = "newly_projected_same_raw_membership"
            elif old_target is None:
                classification = "newly_projected_membership_changed"
            elif old["raw_read_names"] == new["raw_read_names"]:
                classification = (
                    "spelling_changed_same_raw_membership"
                    if old_target["reference_oriented_base"] != new["reference_oriented_base"]
                    else "same_spelling_same_raw_membership"
                )
            else:
                classification = (
                    "spelling_and_membership_changed"
                    if old_target["reference_oriented_base"] != new["reference_oriented_base"]
                    else "membership_changed_same_spelling"
                )
            rows.append({
                "after_record_name": new["record_name"],
                "after_reference_oriented_base": new["reference_oriented_base"],
                "best_before_record_name": old["record_name"] if old else None,
                "before_reference_oriented_base": (
                    old_target["reference_oriented_base"] if old_target else None
                ),
                "before_provenance_category": old["provenance_category"] if old else None,
                "before_sequence_category": old["sequence_category"] if old else None,
                "before_target_candidate_bases": sorted(old_candidate_bases),
                "raw_membership_jaccard": similarity,
                "classification": classification,
            })
        stages[stage] = {
            "before_coordinate_state": before[stage]["coordinate_state"],
            "after_coordinate_state": after[stage]["coordinate_state"],
            "paired_outputs": rows,
        }
    return {
        "coordinate_system": before["coordinate_system"],
        "target": before["target"],
        "input_sha256": before["input_sha256"],
        "before_source": before["assembly_source"],
        "after_source": after["assembly_source"],
        "stages": stages,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    if arguments.output.exists():
        parser.error("output already exists")
    result = compare(
        json.loads(arguments.before.read_text()),
        json.loads(arguments.after.read_text()),
    )
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
