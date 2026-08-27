# Anvaya

Anvaya is a research project exploring damage-aware de Bruijn graph assembly for ancient bacterial and archaeal metagenomes.

The central idea is to retain evidence from the original DNA molecules—such as position within the fragment, orientation, and base quality—and use it to distinguish postmortem-damage-induced graph paths from genuine biological variation.

## Current status

Early Python prototype. Anvaya supports sequence-file input, directed and experimental orientation-aware de Bruijn graph assembly, compacted unitig-graph construction, topology analysis, conservative tip cleaning, simple-bubble detection, weak-tip-to-backbone matching, unitig-level bubble scoring, and non-destructive evidence-based classification of graph alternatives.

Correctness and scientific validation remain the priorities. The orientation-aware graph now uses a compact pure-Python representation that preserves the validated assembly while reducing runtime and memory.

## Initial scope

- ancient bacterial and archaeal metagenomes;
- Illumina short reads;
- untreated double-stranded libraries initially;
- reference-free assembly;
- no damage profile required as input.

## Implemented

- plain and gzipped FASTA/FASTQ input with Phred scores;
- DNA validation, reverse complements, and canonical sequences;
- overlapping k-mer extraction with ambiguous k-mers excluded;
- directed de Bruijn graph construction and edge multiplicity;
- experimental orientation-aware graph construction using canonical nodes and oriented handles;
- rolling 2-bit canonical k-mer encoding and compact integer graph arrays;
- combined forward and reverse-complement k-mer support;
- optional minimum k-mer support filtering;
- experimental conservative tip cleaning for short, weak dead-end paths;
- non-destructive detection of bounded, simple graph bubbles;
- compact read-end evidence collected during orientation-aware graph construction;
- side-specific base qualities retained for molecule-linked terminal observations;
- report-only Phred-derived sequencing-error log likelihoods with explicit missing-quality counts;
- deterministic five-fold held-out damage profiles for matched weak-tip loci;
- report-only damage/error/variation likelihood rankings with explicit model scope and margins;
- sequence-level TSV reports for weak tips and bubble paths;
- damage-compatible C→T/G→A annotation for simple bubble alternatives;
- non-destructive weak-tip matching to equal-length local backbone paths with DNA identity, RY identity, relative coverage, substitutions, and damage compatibility;
- explainable weak-tip evidence scores and conservative damage-like, error-like, variation-like, or ambiguous labels;
- non-destructive detection, backbone matching, and classification of bounded weak branches not represented as standard tips or bubbles;
- linear-time node degree calculation;
- source, sink, and branch identification;
- maximal non-branching unitig extraction;
- oriented unitig-level graph links with aggregate coverage and terminal evidence;
- bounded unitig-level bubble detection with local coverage and sequence-similarity scores;
- labelled unitig-bubble validation across clean, error, damage, and rare-strain simulations;
- unitig-path terminal/internal evidence, terminal enrichment, and strand-balance scores;
- explainable `error-like`, `damage-like`, `variation-like`, and `ambiguous` unitig-path labels;
- graph topology and support summaries;
- `anvaya assemble` command-line workflow and FASTA output;
- deterministic test-only simulation helpers;
- unit tests for assembler primitives.

## Command-line usage

```bash
anvaya assemble -i reads.fastq.gz --k 21 --min-count 2 -o contigs.fasta
```

`--min-count` defaults to `1`, which preserves all observed k-mers. Progress is written to stderr and the final assembly summary to stdout.

The orientation-aware implementation can be tested explicitly:

```bash
anvaya assemble -i reads.fastq.gz --k 31 --min-count 2 \
  --orientation-aware -o contigs.fasta
```

This mode prevents forward and reverse-complement paths from being assembled separately. It remains experimental while graph-cleaning and damage-aware decisions are developed and validated.

Tip cleaning can be enabled for controlled comparisons:

```bash
anvaya assemble -i reads.fastq.gz --k 31 --min-count 2 \
  --orientation-aware --clean-tips -o contigs.fasta
```

Tip cleaning is not yet a default because it must be validated on damaged reads and low-abundance strains.

Simple bubbles can be reported without changing the graph or output unitigs:

```bash
anvaya assemble -i reads.fastq.gz --k 21 --min-count 2 \
  --orientation-aware --clean-tips --detect-bubbles \
  -o contigs.fasta
```

Initial controlled tests found that rare-strain SNPs and sequencing errors can both form simple bubbles, while terminal damage primarily formed tips and incomplete branches. Bubble detection is therefore an analysis primitive, not yet a graph-cleaning or damage-classification rule.

Read-end evidence and event sequences can be reported before cleaning:

```bash
anvaya assemble -i reads.fastq.gz --k 21 --min-count 2 \
  --orientation-aware --end-window 5 \
  --clean-tips --detect-bubbles --event-report events.tsv \
  --damage-profile-report damage-profile.json \
  -o contigs.fasta
```

Evidence reporting is non-destructive. Weak tips and bounded incomplete branches are matched to an equal-length locally competing linear backbone when one is available. The report includes DNA and RY identity, relative coverage, substitutions, oriented terminal evidence, source-read linkage across substitutions, retained base-quality summaries, a Phred-derived sequencing-error log likelihood, three heuristic scores, and an auditable classification. The error likelihood assumes equiprobable wrong nucleotides, `P(observed alternative | sequencing error) = 10^(-Q/10) / 3`, and is omitted rather than imputed when qualities are unavailable. When a damage-profile report is requested, matched damage-compatible events also receive report-only damage, sequencing-error, and biological-variation likelihoods. Weak-tip loci use deterministic five-fold damage profiles fitted without the fold containing the event; incomplete branches use the tip-derived profile and are therefore independent of its training loci. The variation explanation fits one constant local alternative frequency up to 0.5 and receives a one-parameter BIC penalty. Rankings and margins are diagnostic—not posterior probabilities, classifier inputs, or removal decisions.

The optional JSON damage profile compares damage-compatible alternative and backbone observations by the exact changed-base cycle at unique matched weak-tip loci; incomplete branches remain excluded from fitting pending separate profile validation. Schema version 6 retains the schema-5 candidate-conditioned beta-binomial geometric fit and records cross-fit status, fold count, assignment, and scoring scope. The fitted amplitude is not a calibrated whole-library damage rate. Linkage strengthens only multi-substitution evidence; neither linkage nor the fitted profile alone establishes biochemical causality. Reporting does not remove or retain paths automatically and does not require a supplied damage profile.

Compacted-graph bubble paths can be scored for threshold validation:

```bash
anvaya assemble -i reads.fastq.gz --k 21 --min-count 2 \
  --orientation-aware --clean-tips \
  --unitig-bubble-report unitig_bubbles.tsv \
  -o contigs.fasta
```

This report records path support, local coverage, sequence similarity, substitutions, and an evidence-based classification with auditable reasons. Classification is conservative and does not simplify the graph or change the output assembly.

## Testing

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

The clean tip-validation matrix can be run or safely resumed with:

```bash
python3 experiments/03_tip_cleaning_validation.py \
  --reference reference.fna --output-dir results --resume
```

Completed simulation, assembly, and QUAST stages are verified and skipped.

The controlled terminal-evidence validation can be reproduced with:

```bash
PYTHONPATH=src:tests python3 \
  experiments/04_terminal_evidence_ground_truth.py
```

The extended 20-seed unitig-bubble threshold matrix can be reproduced with:

```bash
PYTHONPATH=src:tests python3 \
  experiments/05_unitig_bubble_validation.py
```

The 20-seed weak-tip matching matrix covering separate 5′/3′ damage,
ordinary and terminally concentrated errors, and rare strains can be run with:

```bash
PYTHONPATH=src:tests python3 \
  experiments/06_tip_matching_validation.py
```

The 20-seed incomplete-branch matrix can be reproduced with:

```bash
PYTHONPATH=src:tests python3 \
  experiments/07_incomplete_branch_validation.py
```

Generated candidate and threshold tables are written under `experiments/validation/generated/` and remain ignored by Git.

## Documentation

- [Research question](docs/research_question.md)
- [Initial design](docs/design.md)
- [Benchmark plan](docs/benchmark_plan.md)
- [Literature notes](notes/literature.md)
- [Simulation experiments](experiments/README.md)
