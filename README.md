# Anvaya

Anvaya is a reference-free, damage-aware overlap assembly research prototype
for ancient bacterial and archaeal metagenomes. It retains raw fragment evidence
while exploring conservative overlap clustering, extension, and consensus.

The supported research scope starts with untreated double-stranded libraries
and merged fragments. Paired merging and scaffolding are experimental. General
single-end biological validation and overall superiority to CarpeDeam have not
been demonstrated.

## Installation

Python 3.13 or later is required. The runtime uses the Python standard library.

```bash
python3 -m pip install -e .
anvaya --help
```

## Assembly

The repository includes a reproducible CarpeDeam safe-mode comparator runner.
It validates and fingerprints the external binary, input and damage profiles and
preserves command logs/diagnostics:

```bash
anvaya carpedeam-assemble \
  -i merged.fastq.gz \
  --profile-prefix /path/to/sample_damage_ \
  -o results/carpedeam-safe/contigs.fasta \
  --temporary-directory results/carpedeam-safe/tmp \
  --diagnostics results/carpedeam-safe/diagnostics.json \
  --executable /absolute/path/to/carpedeam --threads 2
```

Safe mode is mandatory. CarpeDeam remains a separate GPL-3.0 executable and is
not installed by Anvaya. This command is benchmark infrastructure, not Anvaya's
assembly engine. Its FASTA does not contain raw-read placements, so the comparison
does not claim Anvaya consensus/provenance. See
[the audited comparator contract](docs/carpedeam_backend.md).

The internal Python assembler remains a bounded research baseline. Use its
explicit conservative configuration:

```bash
anvaya overlap-assemble -i merged.fastq.gz \
  --max-rounds 3 --max-contig-iterations 0 \
  --min-cluster-size 5 --min-output-length 31 \
  --ranked-extension --extension-consensus \
  --reciprocal-best-extension --damage-aware-ranking \
  -o contigs.fasta
```

FASTA and FASTQ inputs may be gzipped. Progress goes to stderr and diagnostics
to stdout. Existing overlap options and defaults are unchanged by the structural
cleanup: a bare invocation is not the conservative recipe above. In particular,
contig merging remains experimental and must be disabled explicitly with
`--max-contig-iterations 0` for this configuration.

Damage-aware ranking currently uses fixed terminal mismatch penalties; it is
not a calibrated probability of assembly correctness. Frozen-layout terminal
polishing is opt-in via `--damage-end-window`; zero disables it. See
`anvaya overlap-assemble --help` for projection and audit options.

For experimental paired input, use `-1 reads_1.fastq.gz -2 reads_2.fastq.gz`
instead of `-i`, optionally with `--merge-overlapping-pairs`. Inputs must already
be synchronized. The current loader checks equal record counts but not matching
read IDs. Unmerged mates do not necessarily expose both physical molecule ends.

## Progressive recovery

The recovery checkpoint uses separate projections and leaves the primary
`contigs.fasta` independent of those outputs. The versioned fixture configuration
in `experiments/overlap_regression.json` records its options. An example matching
the recorded 100k recovery configuration is:

```bash
anvaya overlap-assemble -i merged.fastq.gz \
  --min-cluster-size 5 --max-rounds 3 --min-anchor-matches 1 \
  --min-output-length 31 --max-contig-iterations 0 \
  --ranked-extension --damage-aware-ranking --min-overlap-confidence-margin 0 \
  --reciprocal-best-extension \
  --progressive-raw-phase-audit --max-progressive-raw-iterations 3 \
  --progressive-raw-phase-projection progressive-contigs.fasta \
  --selective-support-rescue-audit --adaptive-rescue-min-support 3 \
  --selective-support-rescue-projection selective-rescue-contigs.fasta \
  --high-confidence-support-two-rescue-audit \
  --high-confidence-support-two-rescue-projection support-two-contigs.fasta \
  --raw-confirmed-master-min-base-quality 20 \
  --raw-confirmed-master-damage-end-window 5 -o contigs.fasta
```

Support-two remains experimental. Initial partner overhangs and later recruit
overhangs require the configured minimum base quality; missing quality fails
this check. Whole-seed rejection is disabled by default after the 100k ablation
showed substantial aligned recovery loss. Overlapping bases still use the
existing count-based consensus, and later recruits do not yet repeat the seed
admission mismatch checks.

For controlled comparisons, enable whole-seed rejection with
`--support-two-seed-quality-filter` or disable later-recruit checking with
`--no-support-two-recruit-quality-filter`. Disabling both restores the earlier
quality policy, including its original checks on initial partner overhangs.
The provisional default is recruit filtering only.

The literature/repository review and prioritized implementation plan are in
[the prototype roadmap](docs/prototype_roadmap.md).

Optional `--support-two-trim-seed-tails` trims exposed unsupported low-quality
runs from rescue contig ends **after extension**, before redundancy projection.
It stops at an independently covered base or a base meeting the quality threshold;
it does not remove internal low-quality bases or trim partner/recruit overhangs.
Missing qualities are treated as insufficient quality. This option cannot be
combined with whole-seed rejection. Results below `--min-output-length` are
excluded from the trimming projection and counted as `trimmed_short_contigs`.

Trimming does not change overlap scoring or read recruitment, so no artificial
trimmed boundary enters the damage-ranking model. The raw read remains intact.
Derived records track a retained original `seed_interval` (0-based, half-open)
and `seed_offset`: contig position = original seed position + offset. These are
in-memory mappings to original **read** ends, not proof of physical molecule ends;
FASTA carries sequences only. Later untracked transformations clear the mapping.

Trimming diagnostics count affected rescue contigs and removed left/right bases
before redundancy projection and the minimum-length filter. Extension diagnostics
describe the assembly before trimming. Run the baseline versus trimming comparison:

```bash
bash experiments/support_two_quality_ablation.sh merged.fastq.gz results/seed-trimming trimming
```

To isolate trimming from changes in redundancy/extension classification, use
`fixed` instead of `trimming`. This runs assembly once, selects output membership
on untrimmed contigs, and writes matched `fixed/before-contigs.fasta` and
`fixed/support-two-contigs.fasta` in the same order. The CLI option is
`--support-two-fixed-membership-before BEFORE_FASTA`, together with trimming and
a support-two projection path. Contigs falling below the minimum length after
trimming are excluded from **both** matched outputs and counted separately as
`fixed_membership_short_exclusions`. This is a paired evaluation subset, not an
unconditional recovery benchmark. MetaQUAST may still change its alignments.

`primary_extension_exclusions` counts rescue candidates excluded because they
extend a primary contig; `ambiguous_primary_extension_exclusions` is a subset
of that count, not another disjoint category. In fixed mode these are baseline
selection decisions; novel/projected base counts describe the trimmed output.
The accounting identity is rescue candidates = primary-contained + rescue-redundant
+ primary-extension-excluded + novel + fixed-membership-short-excluded.

Run the four-way comparison in WSL Bash (the output directory must be new):

```bash
bash experiments/support_two_quality_ablation.sh merged.fastq.gz results/quality-ablation
```

Set `ANVAYA_PYTHON` to the desired Python executable if it is not `python3`.
The script runs identical configurations with neither new filter, seed only,
recruit only, and both, saving FASTAs, diagnostics and timing per configuration.
The user runs this dataset experiment; unit tests use small synthetic fixtures.

`seed_quality_rejections` counts clusters rejected at the seed check, after
the existing partner checks. `seed_quality_rejected_bases` sums the **whole seed
lengths** of those clusters, not just failing bases or lost output bases.
`recruit_quality_rejections` counts rejected candidate evaluations;
`recruit_quality_rejected_overhang_bases` sums their full overhang lengths.
The same read can be evaluated again in another round or cluster, so these are
not unique-read or unique-base counts. `boundary_quality_rejections` retains
its aggregate seed-admission meaning (initial partner failures plus seed failures).

## Raw-fragment consensus experiment

An opt-in consensus projection now uses tracked raw-read placements, original
qualities and a supplied terminal damage profile. It leaves the progressive
layout fixed and emits a separate FASTA. A zero-damage quality control isolates
the profile's contribution. See [model assumptions and outputs](docs/raw_consensus.md).

```bash
bash experiments/raw_consensus_comparison.sh merged.fastq.gz /path/to/profile-prefix results/raw-consensus
```

Use merged untreated double-stranded fragments with retained molecule ends.
The example profile bundled with CarpeDeam is a comparator setting, not a fitted
sample profile. This model remains experimental and has not been validated on the
actual FASTQ benchmark yet.

## Testing and regression

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v -b
PYTHONPATH=src python3 experiments/overlap_regression.py \
  --output-dir results/overlap-baseline
PYTHONPATH=src python3 experiments/overlap_regression.py \
  --output-dir results/overlap-recheck \
  --compare results/overlap-baseline/manifest.json
```

Each regression output directory must be new. Five seeds and six scenarios
exercise the primary, progressive, selective, and support-two FASTAs. Input and
output checksums must match between structural revisions. These small synthetic
fixtures verify reproducibility, not biological accuracy or strain preservation.

The retired `assemble` and `calibrate-events` commands and graph-specific Python
APIs are removed. Existing FASTA `unitig_N` identifiers remain for byte-level
compatibility with saved overlap outputs and correction reports.

## Documentation

- [Research question](docs/research_question.md)
- [Overlap architecture](docs/design.md)
- [Benchmark plan and results](docs/benchmark_plan.md)
- [Audited CarpeDeam backbone](docs/carpedeam_backend.md)
- [Experiment entry points](experiments/README.md)
- [Research history](experiments/archive.md)
- [Literature notes](notes/literature.md)
