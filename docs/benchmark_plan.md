# Overlap benchmark plan

## Structural regression

`experiments/overlap_regression.py` replaces the retired graph-only runner as
the executable regression entry point. `experiments/overlap_regression.json`
contains five seeds and six scenarios: clean, damage, sequencing error, combined
damage/error, rare strain and contamination. References are 5,000 bp, reads are
100 bp in both orientations, and major coverage is 20x.

The runner records input identities and exact FASTA checksums for primary,
progressive, selective and support-two outputs. It uses uniform Q30 qualities
rather than assigning quality from simulated error truth. Configurations and
Python/Git provenance accompany each run. Output directories cannot be reused.

```bash
PYTHONPATH=src python3 experiments/overlap_regression.py \
  --output-dir results/overlap-baseline
PYTHONPATH=src python3 experiments/overlap_regression.py \
  --output-dir results/overlap-recheck \
  --compare results/overlap-baseline/manifest.json
```

This is a behavior-preservation matrix, not an accuracy or performance claim.
It does not execute external comparators or measure false joins. Historical
graph-recovery measurements are not acceptance evidence for the overlap path.

The overlap-only structural cleanup matched all 120 FASTA hashes across the
30 pre/post fixtures. An EMN001 100k rerun also matched all four saved FASTAs
from `results/carpedeam_p0_100k_support_two_rescue/`. Local verification artifacts
are in `results/overlap_only_refactor_before/`, `results/overlap_only_refactor_after/`
and `results/overlap_only_refactor_emn100k/`. The unchanged primary checksum is
`e123731c5fee06d2c6303f9c3b33783a3d7f2ef846eeefab08cbecaa883e284f`;
the support-two checksum is
`f599ac0b6a34cdc651e15c4368b8045f3a9d35eba9171ec1fc73962b46d38429`.

## Recorded EMN001 100k position

These are previously recorded reference-backed results, not outputs of the
small regression fixtures:

| Metric | Progressive | Support-two | CarpeDeam safe |
|---|---:|---:|---:|
| Contigs | 5,388 | 16,137 | 12,491 |
| Total bases | 717,381 | 1,567,417 | 1,534,868 |
| Aligned length | 702,018 | 1,320,332 | 1,431,468 |
| Genome fraction (%) | 0.272 | 0.478 | 0.512 |
| N50 / NA50 | 139 / 137 | 106 / 105 | 134 / 133 |
| Largest contig | 325 | 325 | 660 |
| Mismatches / 100 kbp | 1,272.76 | 988.31 | 873.79 |
| Indels / 100 kbp | 3.99 | 2.80 | 3.21 |
| Unaligned length | 10,647 | 238,103 | 92,737 |
| Reported misassemblies | 0 | 0 | 0 |

The saved report is under
`results/carpedeam_p0_100k_support_two_rescue_metaquast/combined_reference/`.
The input is
`data/carpedeam_p0/simulated/synth-EMN001/upload/synth-EMN001/p0-subsets/EMN001-simulated-100k-min31.fq.gz`;
the reference is the same dataset's `refs/all.fa`. Large inputs and result
artifacts remain unversioned. The recovery command is in the project README.

Support-two reaches 92.2% of CarpeDeam's aligned length but remains more
fragmented. Its 238,103 unaligned bases include 174,252 bases in contigs shorter
than the evaluator's 65 bp minimum alignment. Report these separately from
longer unaligned output; unaligned is not synonymous with false sequence.
Zero reported misassemblies on short contigs does not establish zero false joins.

For reproduction, preserve `--min-contig 31`, `--min-identity 90`, reference
identity, preprocessing and subset identity. Record the actual alignment-length
and ambiguity settings from the evaluator log, plus tool versions. Supplement
historical MetaQUAST settings with short-sequence and origin-truth validation
instead of silently changing the comparator's evaluation.

## Other checkpoints

At 500k, the previously recorded reciprocal/no-contig-merge output reached
N50/NA50 133/132, 4,896,280 aligned bases and zero reported misassemblies.
Damage ranking reached 4,913,114 aligned bases with N50/NA50 133/131. A
raw-confirmed graph projection reached N50/NA50 141/139 and 5,305,903 aligned
bases. These are different configurations; support-two has not been validated
at 500k.

TAF016 paired-library mate ablations validate execution only. R1-only and R2-only
are not genuine single-end biological benchmarks. Support-two N50 remained
76 bp. Reducing the R1 anchor length to 11 reduced validated alignments and
clustered reads; that direction was rejected.

## Next accuracy work

1. Add regression tests and fixes for later-recruit quality and seed-boundary
   policy, without conflating cluster support with per-base support.
2. Add an independent truth evaluator for overlap outputs: false joins,
   ordinary versus damage-compatible errors, attainable reference recovery,
   and strain-specific allele retention. Use 0-based, half-open intervals
   internally, normalize reverse orientations and keep reference IDs separate.
3. Expand held-out simulations to variable fragment lengths, repeats, abundance
   mixtures and realistic stochastic quality/error relationships. Separate
   merged fragments from genuinely truncated single-end input.
4. Compare unchanged overlap configurations, damage-off controls and CarpeDeam
   safe with identical inputs and evaluation settings.
5. Attribute true-overlap loss to candidate discovery, read ownership,
   recruitment and consensus before changing thresholds or adding graph logic.
6. Run the fixed support-two configuration at 500k and profile runtime/RSS.

Accuracy gates require no clean-data regression, no increase in false joins,
no primary-reference loss and less than 1% rare-strain loss using an
allele-sensitive metric. Continuity improvements must include NA50 or aligned
continuity, not N50 alone. Historical 25% damage-scenario N50 targets are
hypotheses for the new accuracy ladder, not satisfied by matching fixture hashes.
Report independent replicate results and do not select thresholds on EMN001
alone. New calibration must state sampling, dependence and multiple-testing
assumptions explicitly.

## Recording a run

Use a new result directory with the full command/configuration, input hashes,
Git revision and dirty state, Python and external-tool versions, stdout/stderr,
timing and peak RSS, output checksums and evaluation reports. Label decisions
accepted, experimental or rejected. Keep generated data outside Git and retain
enough metadata to distinguish projections and reproduce each comparison.
