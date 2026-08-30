# Experiments

This directory contains reproducible research programs. Generated reads,
reports, assemblies, and comparison outputs must not be committed.

## Active acceptance ladder

`18_validation_ladder.py` is the frozen truth-labelled acceptance baseline for
new assembly backends. Its versioned configuration runs five seeds across
clean, terminal-damage, sequencing-error, combined damage/error, rare-strain,
and contamination scenarios at fixed `k=31` and `min_count=2`.

```bash
PYTHONPATH=src:tests python3 experiments/18_validation_ladder.py \
  --output-dir experiments/validation/generated/ladder
```

Use `--dry-run` to inspect the matrix without generating data. The runner
writes a real tab-separated `summary.tsv`, a JSON summary, and a manifest with
configuration, Git, runtime, and checksum provenance.

The acceptance contract is:

- no clean-data regression;
- no increase in false contigs or misassemblies;
- no primary-reference loss;
- less than 1% rare-strain recovery loss;
- at least 25% N50 improvement in damage and combined damage/error scenarios;
- resource use compatible with the public-data environment.

## Current decisions

| Experiment | Decision |
|---|---|
| Support-only and damage-aware tip cleaning | Retain as conservative opt-in comparators; continuity gain is small. |
| Paired extension and direct read threading | Retain as validated research paths; public-data N50 gain is only 38 to 40 bp. |
| P1 dead-end attribution and P2 correction audit | Retain as diagnostics; neither changes topology. |
| Singleton filtered-edge rescue | Rejected; truth-valid edge rescue did not produce safe contig gains. |
| Damage-aware read consensus | Retain as opt-in infrastructure; synthetic gains are strong, but rare-strain loss is 8.98% and public-data N50 remains 38 bp. |
| Next architecture | Iterative whole-fragment, R/Y-aware cluster-and-extend assembly with safe consensus. |

Public datasets remain realism and resource checks unless suitable reference
truth is available. Public-data N50 alone is not an acceptance signal.

## Historical record

The complete chronological experiment notebook, including commands, negative
results, rejected thresholds, and public-data measurements, is preserved in
[archive.md](archive.md). Individual scripts `01` through `17` remain available
to reproduce those results; they are historical evidence, not the active
development roadmap.
