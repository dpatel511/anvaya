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

Use the explicit conservative research configuration:

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

Support-two remains experimental: seed admission checks do not yet enforce
quality on later recruits or all retained seed flanks. The structural cleanup
preserves that behavior so correctness changes can be evaluated separately.

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
- [Experiment entry points](experiments/README.md)
- [Research history](experiments/archive.md)
- [Literature notes](notes/literature.md)
