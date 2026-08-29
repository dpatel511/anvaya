# Anvaya

Anvaya is a research project exploring damage-aware de Bruijn graph assembly for ancient bacterial and archaeal metagenomes.

The central idea is to retain evidence from the original DNA molecules—such as position within the fragment, orientation, and base quality—and use it to distinguish postmortem-damage-induced graph paths from genuine biological variation.

## Current status

Early Python prototype. Anvaya supports sequence-file input, directed and experimental orientation-aware de Bruijn graph assembly, compacted unitig-graph construction, topology analysis, support-only and damage-aware conservative tip cleaning, simple-bubble detection, weak-tip-to-backbone matching, unitig-level bubble scoring, and evidence-based classification of graph alternatives.

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
- experimental damage-aware tip cleaning that removes only one-sided
  error-like tips and protects damage-like, variation-like, ambiguous,
  unmatched, or bidirectionally supported alternatives;
- non-destructive detection of bounded, simple graph bubbles;
- compact read-end evidence collected during orientation-aware graph construction;
- side-specific base qualities retained for molecule-linked terminal observations;
- report-only Phred-derived sequencing-error log likelihoods with explicit missing-quality counts;
- deterministic five-fold held-out damage profiles for matched weak-tip loci;
- report-only damage/error/variation likelihood rankings with explicit model scope and margins;
- fragment-level event-ascertainment conditioning with auditable per-model conditioning probabilities;
- conservative ambiguity for high-quality singleton alternatives that cannot be separated from genuine rare variation;
- report-only graph-context features and held-out separation summaries for topology, path multiplicity, local coverage, fragment support, and terminal evidence;
- direct dominant-path likelihood comparison for equal-length bubble alternatives;
- ordinary-substitution likelihood scoring without inventing a damage channel;
- targeted whole-read molecule and base-quality recovery for ordinary
  substitutions, without retaining observations for every graph edge;
- physical-fragment identities shared by paired reads, with conflicting allele evidence excluded;
- optional report-only conformal confidence calibration with multiple-testing correction and conservative protection decisions;
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

The support-only cleaner is useful as a comparison policy but can erase genuine
low-abundance alternatives. The damage-aware policy adds read-end classification
and bidirectional-support gates:

```bash
anvaya assemble -i reads.fastq.gz --k 31 --min-count 2 \
  --orientation-aware --end-window 5 --damage-aware-clean-tips \
  -o contigs.fasta
```

The two cleaning modes are mutually exclusive and remain opt-in. In a 20-seed
controlled selection matrix, the damage-aware gate selected 82 simulated error
tips and no simulated damage or rare-strain tips. Four reference-backed clean
assemblies at 5× and 20× introduced no misassemblies, mismatches, indels, or
genome-fraction loss. This establishes conservative behavior in those tests,
not universal biological accuracy.

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

Evidence reporting is non-destructive. Weak tips and bounded incomplete branches are matched to an equal-length locally competing linear backbone when one is available; equal-length bubble alternatives are compared directly with the locally dominant bubble path. The report includes DNA and RY identity, relative coverage, substitutions, oriented terminal evidence, physical-fragment linkage across substitutions, retained base-quality summaries, a Phred-derived sequencing-error log likelihood, three heuristic scores, and an auditable classification. Paired mates retain separate read indices but share one fragment identifier, and a fragment with contradictory allele observations is excluded from both support counts. The error likelihood assumes equiprobable wrong nucleotides, `P(observed alternative | sequencing error) = 10^(-Q/10) / 3`, and is omitted rather than imputed when qualities are unavailable. Damage, error, and variation likelihoods are conditioned on the graph event containing at least one consistent alternative fragment and one consistent reference fragment; the report records the conditioning scope and probability for every model. When a damage-profile report is requested, exact C→T/5′ and G→A/3′ events receive all three likelihoods; ordinary substitutions receive error-versus-variation scoring without a damage explanation. Weak-tip loci use deterministic five-fold damage profiles fitted without the fold containing the event; incomplete branches and bubble paths use the tip-derived profile and are therefore independent of its training loci. The variation explanation fits one constant local alternative frequency up to 0.5 and receives a one-parameter BIC penalty. A high-quality alternative seen in only one physical fragment is reported as ambiguous when error otherwise ranks first, because the observations cannot distinguish it from genuine rare variation. Rankings and margins are diagnostic—not posterior probabilities, classifier inputs, or removal decisions.

Ordinary substitutions are no longer limited to the configured damage window.
After event matching, reporting rescans the loaded reads for only the changed
alternative and reference edges, retains one best-quality observation per
physical molecule, and applies the error-versus-variation models to that
whole-read evidence. This targeted second pass avoids an all-edge occurrence
index, but it is intentionally accuracy-first: on the 1,409,072-read public
dataset it increased report time from 171.65 to 547.07 seconds and the complete
run from 504.61 to 853.12 seconds. Assembly without `--event-report` does not
pay this cost.

Classification accuracy is assessed at three distinct levels. Controlled simulations provide exact generating-process labels for damage, independent errors, systematic terminal errors, and rare-strain variants; calibration and validation references use disjoint seeds. Public-data reruns test deterministic non-regression by requiring identical contigs, damage profiles, event counts, and evidence funnels. They do not provide biological ground truth. Claims of real-library damage accuracy therefore remain deferred until independently characterized untreated, partial-UDG, and full-UDG libraries reproduce the expected protocol-specific damage behavior.

The optional JSON damage profile compares damage-compatible alternative and backbone observations by the exact changed-base cycle at unique matched weak-tip loci; incomplete branches remain excluded from fitting pending separate profile validation. Schema version 6 retains the schema-5 candidate-conditioned beta-binomial geometric fit and records cross-fit status, fold count, assignment, and scoring scope. The fitted amplitude is not a calibrated whole-library damage rate. Linkage strengthens only multi-substitution evidence; neither linkage nor the fitted profile alone establishes biochemical causality. Reporting does not remove or retain paths automatically and does not require a supplied damage profile.

Likelihood reports can be calibrated separately using a model trained on independently simulated or validated samples:

```bash
anvaya calibrate-events \
  --input events.tsv \
  --model event-calibration.json \
  --alpha 0.01 \
  --output events-calibrated.tsv
```

The calibrator uses class-conditional conformal p-values, propagates the spread across cross-fit damage curves, and applies Benjamini–Hochberg correction before rejecting the damage and variation explanations. `eligible_error` requires both protective explanations to be rejected, at least two alternative molecules, at least five reference molecules, and stable damage evidence. Damage-like, variation-like, systematic-error-like, and unresolved events are protected or left insufficient. A calibration model is intentionally not bundled: its samples must match the intended library protocol, damage range, coverage, and error processes. The command only appends report fields and never edits the graph or contigs.

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
