"""Compare synthetic assembly results at 0-based reference coordinates.

The gate uses private synthetic truth. It must not be used inside assembly.
"""
import argparse
import json
from collections import Counter
from pathlib import Path


MODES = ("exact", "damage_aware")
STAGES = ("layout", "consensus")
SCHEMA_VERSION = 2
CASE_IDENTITY_FIELDS = (
    "input_sha256",
    "reference_sha256",
    "truth_sha256",
    "damage_profile_sha256",
)
INCREASE_IS_FAILURE = (
    "diagnostic_mosaic_bases",
    "incompatible_candidate_sequences_bases",
    "accuracy_unresolved_bases",
    "nonduplicated_accuracy_discordant_reference_bases",
)
DECREASE_IS_FAILURE = (
    "nonduplicated_accuracy_evaluated_reference_bases",
)


def _case_key(case):
    return case["seed"], case["scenario"], case.get("target_coverage")


def _expected_case_keys(manifest):
    expected = set()
    for case in manifest["expected_cases"]:
        key = _case_key(case)
        if key in expected:
            raise ValueError(f"duplicate expected case: {key}")
        expected.add(key)
    return expected


def _validate_manifest(summary):
    manifest = summary.get("experiment_manifest")
    if not isinstance(manifest, dict):
        raise ValueError("missing experiment manifest")
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"experiment manifest schema must be {SCHEMA_VERSION}")
    for field in ("generator_sha256", "evaluator_sha256", "expected_cases"):
        if not manifest.get(field):
            raise ValueError(f"experiment manifest is missing {field}")
    if tuple(manifest.get("modes", ())) != MODES:
        raise ValueError(f"experiment manifest modes must be {MODES}")
    if tuple(manifest.get("stages", ())) != STAGES:
        raise ValueError(f"experiment manifest stages must be {STAGES}")
    sources = summary.get("source_sha256")
    if not isinstance(sources, dict) or not sources:
        raise ValueError("missing source hash manifest")
    for kind in ("generator", "evaluator"):
        path = manifest.get(f"{kind}_path")
        if not path or sources.get(path) != manifest[f"{kind}_sha256"]:
            raise ValueError(f"{kind} hash does not match source manifest")
    assembly_hash = summary.get("assembly_source_sha256")
    if not isinstance(assembly_hash, str) or len(assembly_hash) != 64:
        raise ValueError("missing valid assembly_source_sha256")
    return manifest, _expected_case_keys(manifest)


def _case_index(summary, expected):
    cases = {}
    for case in summary["cases"]:
        key = _case_key(case)
        if key in cases:
            raise ValueError(f"duplicate case: {key}")
        for field in CASE_IDENTITY_FIELDS:
            value = case.get(field)
            if not isinstance(value, str) or len(value) != 64:
                raise ValueError(f"case {key} is missing valid {field}")
        cases[key] = case
    if set(cases) != expected:
        missing = sorted(expected - set(cases))
        extra = sorted(set(cases) - expected)
        raise ValueError(f"cases differ from manifest; missing={missing}, extra={extra}")
    return cases


def _state_index(metrics):
    states = {}
    for state in metrics["nonduplicated_accuracy_coordinate_states"]:
        key = state["reference"], state["position_0based"]
        if key in states:
            raise ValueError(f"duplicate coordinate state: {key}")
        states[key] = state
    return states


def compare_stage(before, after):
    """Return coverage-aware transitions and prospective safety failures."""
    before_states = _state_index(before)
    after_states = _state_index(after)
    missing = {"status": "not_uniquely_evaluated", "reason": "not_projected"}
    transitions = Counter()
    changed_sites = []
    failures = []
    correct_gain = 0
    for reference, position in sorted(before_states.keys() | after_states.keys()):
        old = before_states.get((reference, position), missing)
        new = after_states.get((reference, position), missing)
        if (
            "reference_base" in old and "reference_base" in new
            and old["reference_base"] != new["reference_base"]
        ):
            raise ValueError(f"reference base differs at {(reference, position)}")
        transition = f"{old['status']}->{new['status']}"
        transitions[transition] += 1
        if old["status"] != new["status"] or old.get("observed_bases") != new.get("observed_bases"):
            changed_sites.append({
                "reference": reference,
                "position_0based": position,
                "before": old,
                "after": new,
            })
        if old["status"] == "concordant_correct" and new["status"] != "concordant_correct":
            failures.append(f"previously_correct_coordinate_regressed:{reference}:{position}")
        if old["status"] != "discordant" and new["status"] == "discordant":
            failures.append(f"new_discordant_coordinate:{reference}:{position}")
        if old["status"] == "not_uniquely_evaluated" and new["status"] == "concordant_incorrect":
            failures.append(f"newly_evaluated_incorrect_coordinate:{reference}:{position}")
        if old["status"] != "concordant_correct" and new["status"] == "concordant_correct":
            correct_gain += 1

    for metric in INCREASE_IS_FAILURE:
        if after[metric] > before[metric]:
            failures.append(f"{metric}_increased")
    for metric in DECREASE_IS_FAILURE:
        if after[metric] < before[metric]:
            failures.append(f"{metric}_decreased")
    before_denominator = before["nonduplicated_accuracy_evaluated_reference_bases"]
    after_denominator = after["nonduplicated_accuracy_evaluated_reference_bases"]
    if before_denominator and not after_denominator:
        failures.append("nonduplicated_mismatch_rate_not_evaluable")
    elif before_denominator and after_denominator and (
        after["nonduplicated_accuracy_mismatches"] * before_denominator
        > before["nonduplicated_accuracy_mismatches"] * after_denominator
    ):
        failures.append("nonduplicated_mismatch_rate_increased")

    return {
        "coordinate_system": "0-based reference positions; intervals half-open",
        "transitions": dict(sorted(transitions.items())),
        "changed_sites": changed_sites,
        "correct_coordinate_gain": correct_gain,
        "n50": [before["n50"], after["n50"]],
        "longest": [before["longest"], after["longest"]],
        "failures": sorted(set(failures)),
    }


def compare_summaries(before, after):
    if before.get("seed_role") != after.get("seed_role"):
        raise ValueError("seed roles differ")
    before_manifest, before_expected = _validate_manifest(before)
    after_manifest, after_expected = _validate_manifest(after)
    if before_manifest != after_manifest:
        raise ValueError("experiment manifests differ")
    if before["assembly_source_sha256"] == after["assembly_source_sha256"]:
        raise ValueError("baseline and candidate assembly source hashes are identical")
    before_cases = _case_index(before, before_expected)
    after_cases = _case_index(after, after_expected)

    rows = []
    for key in sorted(before_cases):
        old_case, new_case = before_cases[key], after_cases[key]
        for field in (
            "panel", "damage", "sequencing_error_rate", "supplied_profile_scale",
            "reads", "minimum_read_length", "maximum_read_length",
            "simulated_damage_events", "simulated_sequencing_errors",
            *CASE_IDENTITY_FIELDS,
        ):
            if old_case.get(field) != new_case.get(field):
                raise ValueError(f"case {key} differs in {field}")
        if set(old_case["results"]) != set(MODES) or set(new_case["results"]) != set(MODES):
            raise ValueError(f"case {key} does not contain exactly {MODES}")
        for mode in MODES:
            for stage in STAGES:
                comparison = compare_stage(
                    old_case["results"][mode][stage],
                    new_case["results"][mode][stage],
                )
                baseline_seconds = old_case["results"][mode].get("assembly_seconds")
                candidate_seconds = new_case["results"][mode].get("assembly_seconds")
                comparison["assembly_seconds"] = [
                    baseline_seconds, candidate_seconds,
                ]
                if (
                    baseline_seconds is not None
                    and candidate_seconds is not None
                    and candidate_seconds > 2 * baseline_seconds
                ):
                    comparison["failures"].append("assembly_time_exceeded_2x")
                rows.append({
                    "seed": key[0], "scenario": key[1],
                    "target_coverage": key[2], "mode": mode,
                    "stage": stage, "damage": old_case["damage"], **comparison,
                })

    safety_failures = [row for row in rows if row["failures"]]
    changed_rows = [row for row in rows if row["changed_sites"] or row["n50"][0] != row["n50"][1]]
    global_failures = []
    control_panel = str(before.get("seed_role", "")).endswith("_error_only")
    if not control_panel:
        if not changed_rows:
            global_failures.append("candidate_is_no_op")
        seeds = sorted({key[0] for key in before_cases})
        for seed in seeds:
            clean_layouts = [
                row for row in rows
                if row["seed"] == seed and row["scenario"] == "clean_unique"
                and row["target_coverage"] == 20 and row["stage"] == "layout"
            ]
            if {row["mode"] for row in clean_layouts} != set(MODES) or any(
                row["n50"][1] * 2 < row["n50"][0] * 3 for row in clean_layouts
            ):
                global_failures.append(f"insufficient_clean_20x_n50_gain:{seed}")
            damaged_cases = {}
            for row in rows:
                if row["seed"] == seed and row["damage"] > 0 and row["mode"] == "damage_aware":
                    damaged_cases.setdefault(
                        (row["scenario"], row["target_coverage"]), {}
                    )[row["stage"]] = row
            damaged_benefit = any(
                set(stages) == set(STAGES)
                and not stages["layout"]["failures"]
                and not stages["consensus"]["failures"]
                and stages["layout"]["n50"][1] > stages["layout"]["n50"][0]
                and stages["consensus"]["n50"][1] >= stages["consensus"]["n50"][0]
                for stages in damaged_cases.values()
            )
            if not damaged_benefit:
                global_failures.append(f"no_damage_aware_contiguity_benefit:{seed}")
    return {
        "coordinate_system": "0-based reference positions; intervals half-open",
        "schema_version": SCHEMA_VERSION,
        "case_combinations": len(before_cases),
        "distinct_reference_seeds": sorted({key[0] for key in before_cases}),
        "correlated_mode_stage_rows": len(rows),
        "rows": rows,
        "failed_rows": safety_failures,
        "global_failures": global_failures,
        "control_panel": control_panel,
        "passed": not safety_failures and not global_failures,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    if arguments.output.exists():
        parser.error("output already exists")
    result = compare_summaries(
        json.loads(arguments.before.read_text()),
        json.loads(arguments.after.read_text()),
    )
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
