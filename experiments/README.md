# Experiments

The active assembly implementation is overlap-only. Generated data, reports
and assemblies belong in ignored `results/` or `generated/` directories.

| Entry point | Purpose |
|---|---|
| `overlap_regression.py` | Deterministic five-seed/six-scenario CLI regression with before/after FASTA hash comparison |
| `overlap_regression.json` | Versioned fixture data and recovery configuration |
| `19_overlap_correction_truth_audit.py` | Reference-assisted, one-edit-at-a-time audit of accepted polishing corrections |
| `20_restore_unpolished_overlap.py` | Restore an exact topology-matched pre-polish FASTA from its correction report |

Run from the repository root:

```bash
PYTHONPATH=src python3 experiments/overlap_regression.py \
  --output-dir results/overlap-baseline
PYTHONPATH=src python3 experiments/overlap_regression.py \
  --output-dir results/overlap-recheck \
  --compare results/overlap-baseline/manifest.json
```

The runner reuses the deterministic simulator in `tests/helpers/simulation.py`.
Its checksums validate structural changes, not biological accuracy. It reports
contig counts and total bases for context and deliberately does not label
non-exact contigs as false. Output directories must be new.

Reference truth is used only by the external correction audit, never by the
assembler. The correction report uses 0-based base offsets; PAF intervals are
0-based and half-open. Resolved edit precision excludes unresolved mappings and
must be reported with its resolved fraction.

Contig merging, support-two rescue and paired scaffolding remain experimental.
See the [benchmark plan](../docs/benchmark_plan.md) for measured results,
limitations and pending gates. The [research history](archive.md) records why
earlier assembly directions were retired.
