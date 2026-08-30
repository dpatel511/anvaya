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
- report-only damage-projection candidates that require concordant classification,
  held-out likelihood, variation contrast, and complete quality evidence;
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
- optional reciprocal paired-read unitig extension with strong-winner and
  bounded-work ambiguity guards;
- optional direct read-thread extension with molecule-level strong-winner,
  reciprocal-orientation, and local repeat-coverage guards;
- bounded unitig-level bubble detection with local coverage and sequence-similarity scores;
- labelled unitig-bubble validation across clean, error, damage, and rare-strain simulations;
- unitig-path terminal/internal evidence, terminal enrichment, and strand-balance scores;
- explainable `error-like`, `damage-like`, `variation-like`, and `ambiguous` unitig-path labels;
- graph topology and support summaries;
- `anvaya assemble` command-line workflow and FASTA output;
- deterministic test-only simulation helpers;
- unit tests for assembler primitives.
- report-only raw-read attribution of compacted-graph dead ends;
- report-only solid-k-mer correction-candidate auditing with explicit
  protection of terminal damage-compatible reversals.

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

Compacted-graph boundaries can be attributed without changing the FASTA:

```bash
anvaya assemble -i reads.fastq.gz --k 31 --min-count 2 \
  --orientation-aware --fragmentation-report fragmentation.tsv \
  -o contigs.fasta
```

The TSV contains one row per physical unitig end and labels its retained graph
topology as `dead_end`, `unique_continuation`, or `ambiguous_branch`, together
with candidate handles and unitig support. It also labels each physical unitig
as `isolated`, `one_sided`, or `connected`, making disconnected sequence easy
to distinguish from branch-limited sequence. These labels describe the
retained graph only; they do not claim that a missing edge was caused by
filtering, damage, or insufficient coverage.

Raw reads can be rescanned to attribute those dead ends without changing the
graph or FASTA:

```bash
anvaya assemble -i reads.fastq.gz --k 31 --min-count 2 \
  --orientation-aware --dead-end-report dead-ends.tsv \
  -o contigs.fasta
```

The P1 report distinguishes read boundaries, apparent coverage gaps, unique
continuations removed by support filtering, conflicting filtered
continuations, and retained contexts observed elsewhere. It also records
terminal, low-quality, and missing-quality evidence for each raw extension.

P2 audits conservative solid-k-mer substitution candidates:

```bash
anvaya assemble -i reads.fastq.gz --k 31 --min-count 2 \
  --orientation-aware --correction-report corrections.tsv \
  --correction-max-quality 20 --correction-max-reads 20000 \
  -o contigs.fasta
```

Only low-quality positions are considered. A candidate is marked
`would_correct` only when one substitution uniquely restores every overlapping
k-mer to solid support. Terminal T-to-C and A-to-G reversals compatible with
ancient-DNA damage are protected, and ambiguous or incomplete rescues are
reported separately. Both reports are diagnostic: they do not edit reads,
graph topology, unitigs, or output contigs.

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

Paired-end reads can optionally resolve a conservative subset of compacted-graph
junctions after damage-aware cleaning:

```bash
anvaya assemble -1 left.fastq.gz -2 right.fastq.gz --k 31 --min-count 2 \
  --orientation-aware --end-window 5 --damage-aware-clean-tips \
  --paired-unitig-extension -o contigs.fasta
```

The default rule requires five independent pairs, 3:1 winner-to-runner-up
dominance, reciprocal orientation agreement, and at most 1,000 searched graph
states per junction. Work-limited or ambiguous junctions are left unresolved.
These defaults are deliberately stricter than the initial 3-pair/2:1 prototype:
on a clean 20× reference-backed validation the initial rule introduced small
base-level errors, while the strict rule increased N50 from 4,665 to 4,733 bp
and genome fraction from 92.952% to 93.014% with zero QUAST misassemblies,
mismatches, or indels. On the unlabeled `SRR32866683` library, paired extension
made only 17 joins at `k=21` and did not change its 32 bp N50, so this remains
an experimental repeat resolver rather than a general continuity solution.

Source reads can instead be threaded directly across compacted-graph junctions:

```bash
anvaya assemble -1 left.fastq.gz -2 right.fastq.gz --k 31 --min-count 2 \
  --orientation-aware --end-window 5 --damage-aware-clean-tips \
  --read-thread-extension --read-thread-report read-threads.tsv \
  -o contigs.fasta
```

The default policy requires five independent molecules, 3:1 dominance, and
reciprocal orientation agreement. It also breaks the weaker of two adjacent
joins when their internal unitig has at least twice the mean edge support of
both flanks, a local copy-number signal controlled by
`--thread-repeat-coverage-ratio`. The report-only form can be used without
`--read-thread-extension` and does not change the FASTA.

On the clean 20× `GCF_000007145.1` reference-backed dataset, direct threading
reduced 6,342 unitigs to 3,942 contigs, increased N50 from 4,856 to 18,021 bp,
increased the largest contig from 20,336 to 81,176 bp, and recovered 97.260% of
the reference. QUAST reported zero misassemblies, mismatches, and indels. The
coverage guard rejected one repeat-like join that otherwise introduced a 7 bp
indel, without reducing N50. A 20-run short-fragment, damage, sequencing-error,
related-strain, and contamination matrix accepted 6,989 links with 100% R/Y
damage-tolerant topology precision and no increase in false contigs.

On public `SRR32866683`, threading at `k=31` joined 32,663 links, reduced
639,369 unitigs to 606,706 contigs, increased N50 from 38 to 40 bp, and raised
the largest contig from 366 to 435 bp. The coverage guard rejected only six of
32,669 reciprocal candidates and preserved the v20 continuity metrics. This
library lacks assembly truth, so those numbers establish execution and
reference-free topology changes—not biological correctness. Read threading
therefore remains experimental and opt-in.

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

Evidence reporting is non-destructive. Weak tips and bounded incomplete branches are matched to an equal-length locally competing linear backbone when one is available; equal-length bubble alternatives are compared directly with the locally dominant bubble path. The report includes DNA and RY identity, relative coverage, substitutions, oriented terminal evidence, physical-fragment linkage across substitutions, retained base-quality summaries, a Phred-derived sequencing-error log likelihood, three heuristic scores, and an auditable classification. Paired mates retain separate read indices but share one fragment identifier, and a fragment with contradictory allele observations is excluded from both support counts. The error likelihood assumes equiprobable wrong nucleotides, `P(observed alternative | sequencing error) = 10^(-Q/10) / 3`, and is omitted rather than imputed when qualities are unavailable. Damage, error, and variation likelihoods are conditioned on the graph event containing at least one consistent alternative fragment and one consistent reference fragment; the report records the conditioning scope and probability for every model. When a damage-profile report is requested, exact C→T/5′ and G→A/3′ events receive all three likelihoods; ordinary substitutions receive error-versus-variation scoring without a damage explanation. Weak-tip loci use deterministic five-fold damage profiles fitted without the fold containing the event; incomplete branches and bubble paths use the tip-derived profile and are therefore independent of its training loci. The variation explanation fits one constant local alternative frequency up to 0.5 and receives a one-parameter BIC penalty. A high-quality alternative seen in only one physical fragment is reported as ambiguous when error otherwise ranks first, because the observations cannot distinguish it from genuine rare variation. A `would_project` report row additionally requires all of the damage-like classification, held-out likelihood, variation contrast, and base-quality gates; it is an audit signal only and does not alter graph topology or FASTA output. Rankings and margins are diagnostic—not posterior probabilities, classifier inputs, or removal decisions.

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

The frozen P0 validation ladder runs deterministic truth-labelled controls for
clean reads, terminal damage, sequencing error, combined damage and error,
rare-strain mixtures, and contamination. It records contig continuity,
reference recovery, false contigs, damage-event counts, runtime provenance, and
the exact benchmark configuration:

```bash
PYTHONPATH=src:tests python3 experiments/18_validation_ladder.py \
  --output-dir experiments/validation/generated/ladder
```

Use `--dry-run` to inspect the frozen matrix without generating data. The
versioned configuration is `experiments/benchmark_ladder.json`; generated
summaries and manifests remain ignored. This ladder is the acceptance baseline
for future assembly changes. Public data remains a realism/performance check,
not exact biological truth.

For collaborative validation, Codex supplies the exact test command and the
user runs it in the target environment. Codex then analyzes the saved console
output and generated result files; test execution by Codex is not assumed.

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
