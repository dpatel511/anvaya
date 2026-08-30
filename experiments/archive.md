# Experiments

Generated FASTQ files and result directories must not be committed. Publication-level experiments will use established simulators with recorded versions, commands, configurations, seeds, and reference checksums.

## Frozen P0 validation ladder

`18_validation_ladder.py` is the reproducible acceptance baseline for future
assembly changes. Its versioned configuration runs five seeds across clean,
terminal-damage, sequencing-error, combined-damage/error, rare-strain, and
contamination scenarios at fixed `k=31` and `min_count=2`. It writes
truth-labelled continuity and recovery metrics plus a manifest with the
configuration checksum, Git state, runtime, and output checksums.

```bash
PYTHONPATH=src:tests python3 experiments/18_validation_ladder.py \
  --output-dir experiments/validation/generated/ladder
```

Use `--dry-run` to inspect the six-scenario, five-seed matrix without running
assembly. The generated ladder directory remains ignored by Git. Established
assemblers are comparator slots in the frozen plan; they must be run on the
same generated reads and reported separately from Anvaya's exact-truth metrics.

The initial 30-run baseline (current orientation-aware graph, `k=31`,
`min_count=2`) was:

| Scenario | Mean N50 (bp) | Mean largest contig (bp) | Mean major-reference fraction | False contigs (total) |
|---|---:|---:|---:|---:|
| Clean | 4,980 | 4,980 | 1.0000 | 0 |
| Sequencing error | 1,315 | 1,890 | 0.7700 | 33 |
| Damage | 472 | 991 | 0.7400 | 87 |
| Damage + error | 292 | 966 | 0.7200 | 120 |
| Rare strain | 69 | 335 | 0.8900 | 120 |
| Contamination | 377 | 991 | 0.7400 | 87 |

These are baseline measurements, not claims of biological accuracy for the
public library. The clean-versus-damaged contrast confirms that terminal
damage is a high-value continuity target for the next phase, while the
rare-strain row supplies the protection gate against simply deleting shorter
paths.

## Graph-disturbance development test

`01_graph_disturbance.py` is a deterministic implementation check, not a scientific benchmark. It uses a random 500 bp reference, 60 bp reads placed every 10 bp, `k=21`, and terminal damage rates of 0.50, 0.30, 0.15, 0.08, and 0.04 from each end.

```bash
PYTHONPATH=src:tests python3 experiments/01_graph_disturbance.py
```

Results for seed 42:

| Condition | Nodes | Edges | Branches | Unitigs | Longest unitig |
|---|---:|---:|---:|---:|---:|
| Clean | 481 | 480 | 0 | 1 | 500 bp |
| Terminal damage | 531 | 530 | 24 | 49 | 85 bp |

The clean reads exactly reconstruct the reference. Twenty-seven terminal damage events fragment the graph and reduce the longest unitig from 500 bp to 85 bp. This validates the experimental plumbing only.

## Graph-traversal optimization test

The CLI was run with `k=21` on a local gzipped FASTQ containing 19,643 reads. This measured implementation performance, not assembly accuracy.

| Measurement | Before | After |
|---|---:|---:|
| Total runtime | 1324.06 s | 8.74 s |
| Unitig extraction | 383.00 s | 0.24 s |
| Peak memory after optimization | — | 132,724 KiB |

The optimized run was approximately 151 times faster. Both runs produced identical graph statistics and byte-identical FASTA output. The improvement came from calculating all node degrees in one pass, using constant-time degree lookups, and avoiding duplicate unitig extraction.

## Baseline accuracy experiment

### Question

How well does the current non-damage-aware Anvaya prototype reconstruct a known genome from ordinary Illumina-like reads compared with mature metagenomic assemblers?

This is a baseline test of the current assembler, not an ancient-DNA or metagenomic benchmark.

### Setup

- Reference: `GCF_000007145.1.fna`, 5,076,188 bp, 65.07% GC.
- Simulation: ART paired 150 bp HiSeq 2500 reads, 20× coverage, 300 ± 30 bp fragments, seed 42.
- Input: 676,696 FASTQ reads across two mate files.
- Anvaya: `k=21`; mates were pooled without using pairing links.
- Comparators: MEGAHIT and metaSPAdes with four threads; metaSPAdes used `--only-assembler`.
- Evaluation: QUAST 5.3.0 against the known reference, minimum contig length 200 bp.

Run the reproducible workflow with:

```bash
python3 experiments/02_baseline_accuracy.py \
  --reference ../annostat/data/GCF_000007145.1.fna \
  --output-dir experiments/baseline/results/run_01
```

### Accuracy results

| Metric | Anvaya | MEGAHIT | metaSPAdes |
|---|---:|---:|---:|
| Contigs ≥200 bp | 4,620 | 106 | 103 |
| Largest contig | 657 bp | 448,867 bp | 449,487 bp |
| Total length ≥200 bp | 1,182,280 bp | 5,011,570 bp | 4,993,374 bp |
| N50 | 249 bp | 151,608 bp | 130,616 bp |
| Genome fraction | 21.444% | 98.601% | 98.347% |
| Misassemblies | 0 | 2 | 0 |
| Mismatches per 100 kbp | 0.68 | 2.02 | 3.74 |
| Duplication ratio | 1.086 | 1.001 | 1.000 |

Anvaya produced 547,833 unitigs in total, but only 4,620 were at least 200 bp and none reached 1,000 bp. Its low mismatch and misassembly counts do not indicate a superior assembly: the contigs are too short and recover only 21.444% of the reference.

### Performance results

| Tool | Wall time | Peak memory |
|---|---:|---:|
| Anvaya | 4 min 39 s | 7.1 GiB |
| MEGAHIT | 38 s | 299 MiB |
| metaSPAdes | 1 min 20 s | 1.5 GiB |

Anvaya constructed 13,424,232 nodes and 13,584,362 edges from 87,970,480 k-mer observations. The large graph and fragmented output are consistent with unsupported sequencing-error paths being retained.

### What this establishes

- The simulation, assembly, timing, and reference-based evaluation workflow works end to end.
- The current graph preserves accurate short local sequences but fails to recover a complete or contiguous genome.
- The current prototype is slower and substantially more memory intensive than both comparison tools.
- This is not a feature-matched comparison: Anvaya lacks graph cleaning, orientation-aware construction, multi-k assembly, and paired-link traversal.
- Optional low-support filtering was selected as the next controlled change and is evaluated below.

### Reports

Generated reports remain ignored under `experiments/baseline/results/run_01/`:

- `quast/report.txt`: complete comparison table;
- `quast/report.html`: interactive QUAST report;
- `timing/*.txt`: wall time and peak-memory measurements;
- `commands.txt`: exact commands used;
- `manifest.json`: parameters, git state, runtime, resolved tool paths, command
  arguments, and SHA-256 checksums for available inputs and primary outputs;
- `*.log`: tool output and diagnostics.

QUAST reported a non-fatal Minimap2 version warning. Tool versions should be pinned before publication-level experiments.

## Minimum-support filtering experiment

### Question

Does removing k-mers with weak abundance support improve the baseline without introducing assembly errors?

The same reads, reference, `k=21`, and QUAST settings were used for `min_count` values 1, 2, and 3. A value of 1 is the unchanged baseline.

```bash
anvaya assemble -1 simulated1.fq -2 simulated2.fq \
  --k 21 --min-count 2 -o contigs.fasta
```

### Results

| Metric | min=1 | min=2 | min=3 |
|---|---:|---:|---:|
| Graph nodes | 13,424,232 | 9,902,440 | 9,811,260 |
| Branching nodes | 348,205 | 13,665 | 12,289 |
| Total unitigs | 547,833 | 24,336 | 26,989 |
| Contigs ≥200 bp | 4,620 | 8,332 | 10,108 |
| Largest contig | 657 bp | 10,197 bp | 7,929 bp |
| N50 | 249 bp | 1,698 bp | 1,299 bp |
| Genome fraction | 21.444% | 94.940% | 94.870% |
| Duplication ratio | 1.086 | 2.005 | 1.978 |
| Misassemblies | 0 | 0 | 0 |
| Mismatches per 100 kbp | 0.68 | 0.01 | 0.00 |
| Wall time | 4 min 39 s | 3 min 15 s | 3 min 19 s |
| Peak memory | 7.1 GiB | 5.83 GiB | 5.79 GiB |

### Interpretation

`min_count=2` produced the best accuracy-contiguity trade-off. It increased genome recovery from 21.444% to 94.940%, increased N50 nearly sevenfold, reduced unitig fragmentation by about 96%, reduced runtime by about 30%, and reduced peak memory by about 18%. Threshold 3 removed more connecting k-mers and reduced contiguity without improving genome recovery.

The approximately 2.0 duplication ratio after filtering was consistent with forward and reverse-complement paths being assembled separately. This result motivated the orientation-aware and compact graph experiments below.

The full comparison remains ignored under `experiments/baseline/results/min_count_comparison/quast/`.

## Orientation-aware graph experiment

### Question

Does representing each k-mer and its reverse complement as one physical graph edge remove strand duplication without reducing assembly quality?

The baseline simulated reads were assembled with `k=31` and `min_count=2`. The directed and orientation-aware runs used identical reads and evaluation settings.

```bash
anvaya assemble -1 simulated1.fq -2 simulated2.fq \
  --k 31 --min-count 2 --orientation-aware -o contigs.fasta
```

### Results

| Metric | Directed | Orientation-aware |
|---|---:|---:|
| Graph nodes | 9,923,379 | 5,006,437 |
| Branching nodes | 4,669 | 3,756 |
| Total unitigs | 10,708 | 6,355 |
| Contigs >=200 bp | 4,135 | 1,664 |
| Largest contig | 28,204 bp | 20,336 bp |
| N50 | 3,804 bp | 4,813 bp |
| NGA50 | 5,874 bp | 4,715 bp |
| Genome fraction | 97.309% | 96.904% |
| Duplication ratio | 1.999 | 1.005 |
| Misassemblies | 0 | 0 |
| Mismatches per 100 kbp | 0.01 | 0.00 |
| Wall time | 202.07 s | 193.90 s |
| Peak memory | 6.13 GiB | 6.92 GiB |

### Interpretation

Canonical physical edges reduced the node count by about half and corrected the duplication ratio from 1.999 to 1.005. N50 improved by about 27%, and the assembly remained free of detected misassemblies. This validates the orientation-aware graph model.

The result is not yet an overall performance improvement. Genome fraction decreased by 0.405 percentage points, the largest contig and NGA50 decreased, and peak memory increased by about 13%. The current implementation stores two oriented adjacency maps per node, Python tuple edge keys, counters, and one support object per physical edge. Graph construction also repeatedly creates and reverse-complements string k-mers.

The initial orientation-aware mode therefore remained experimental pending a compact representation. That optimization is evaluated below.

The complete generated QUAST report remains ignored under `experiments/baseline/results/orientation_aware/quast/`.

## Compact orientation-aware graph experiment

### Question

Can compact graph storage reduce runtime and memory without changing the validated orientation-aware assembly?

The same 676,696 reads were assembled with `k=31`, `min_count=2`, and `--orientation-aware`. The compact implementation uses rolling 2-bit canonical k-mers, packed strand support, typed integer arrays, and byte-array traversal state.

### Results

| Metric | Initial orientation-aware | Compact orientation-aware | Change |
|---|---:|---:|---:|
| Graph construction | 161.26 s | 92.73 s | 43% faster |
| Unitig extraction | 25.01 s | 16.97 s | 32% faster |
| Wall time | 193.90 s | 113.99 s | 41% faster |
| Peak memory | 6.92 GiB | 2.78 GiB | 60% lower |
| Graph nodes | 5,006,437 | 5,006,437 | unchanged |
| Total unitigs | 6,355 | 6,355 | unchanged |
| N50 | 4,813 bp | 4,813 bp | unchanged |
| Genome fraction | 96.904% | 96.904% | unchanged |
| Duplication ratio | 1.005 | 1.005 | unchanged |
| Misassemblies | 0 | 0 | unchanged |

### Interpretation

The optimization met its goal. Peak memory fell below the directed `k=31` baseline of 6.13 GiB, and wall time improved by about 1.7-fold. The canonical unitig multisets were identical even though FASTA ordering and strand presentation differed. QUAST reported identical values for every assembly metric, confirming that the storage change did not alter accuracy.

The implementation remains CPU-bound and slower than mature assemblers. Further optimization should be guided by profiling; stable hot loops can later be moved to Cython without rewriting the complete Python workflow.

Generated timing data and the QUAST comparison remain ignored under `experiments/baseline/results/compact_orientation/`.

## Conservative tip-cleaning experiment

### Question

Can short, weak dead-end paths be removed to improve contiguity without sacrificing reference recovery, accuracy, or memory efficiency?

The compact orientation-aware baseline was rerun on the same 676,696 reads with `k=31`, `min_count=2`, and `--clean-tips`. Cleaning removes paths no longer than `2 × k` only when their strongest edge has at most 20% of competing support.

```bash
anvaya assemble -1 simulated1.fq -2 simulated2.fq \
  --k 31 --min-count 2 --orientation-aware --clean-tips \
  -o contigs.fasta
```

### Results

| Metric | Compact baseline | Tip cleaning |
|---|---:|---:|
| Active graph nodes | 5,006,437 | 4,993,920 |
| Active graph edges | 5,007,930 | 4,995,413 |
| Branching nodes | 3,756 | 2,967 |
| Total unitigs | 6,355 | 4,760 |
| Contigs ≥200 bp | 1,664 | 1,081 |
| Largest contig | 20,336 bp | 46,465 bp |
| N50 | 4,813 bp | 8,056 bp |
| NGA50 | 4,715 bp | 7,876 bp |
| Genome fraction | 96.904% | 97.022% |
| Duplication ratio | 1.005 | 1.001 |
| Misassemblies | 0 | 0 |
| Wall time | 113.99 s | 117.52 s |
| Peak memory | 2.78 GiB | 2.77 GiB |

### Interpretation

Cleaning removed 806 tips containing 12,517 edges in 0.92 seconds. N50 increased by 67%, the largest contig increased by 128%, unitig count fell by 25%, and genome fraction increased by 0.118 percentage points. QUAST reported no misassemblies, mismatches, or indels. Peak memory remained effectively unchanged, while total runtime increased by about 3%.

This was a practically important result but not evidence of general or statistical significance by itself. The clean multi-dataset validation below tests whether the behaviour reproduces across genomes, seeds, and coverage levels.

Generated timing data and the QUAST comparison remain ignored under `experiments/baseline/results/tip_cleaning/`.

## Clean tip-validation matrix

### Question

Does conservative tip cleaning behave consistently across bacterial genomes, simulation seeds, and sequencing coverage?

The current assembler was tested unchanged on two bacterial references, seeds 42 and 43, and 5× and 20× ART coverage. Each of the eight datasets was assembled with tip cleaning disabled and enabled, then evaluated with QUAST.

```bash
python3 experiments/03_tip_cleaning_validation.py \
  --reference reference_1.fna --reference reference_2.fna \
  --output-dir experiments/validation/results/clean_pilot \
  --seeds 42 43 --coverages 5 20
```

Use `--resume` to verify outputs, skip completed stages, and regenerate the combined summary without repeating completed work.

### Results

| Coverage | Conditions | N50 effect | Largest-contig effect | Genome fraction | Misassemblies |
|---|---:|---:|---:|---:|---:|
| 5× | 4 | 0–0.14% increase | unchanged in 3/4; 36% increase in 1/4 | increased by 0.011–0.019 points | 0 in all runs |
| 20× | 4 | 56–67% increase | 52–128% increase | increased by 0.100–0.138 points | 0 in all runs |

At 20×, duplication ratio also improved from 1.005 to 1.001–1.002. At 5×, only 34–116 tips were removed per dataset, compared with 727–1,329 at 20×. Peak memory was effectively unchanged, and runtime showed no consistent cleaning overhead beyond normal run-to-run variation.

### Interpretation

The rule behaves conservatively when coverage is weak and reproducibly improves contiguity when support separation is stronger. No clean validation condition lost genome fraction or gained a reported misassembly. This strengthens the clean-isolate evidence, but the sample remains too small for a broad statistical claim.

The remaining safety tests are low-abundance related strains and realistic ancient-DNA damage. Tip cleaning remains opt-in until those tests show that weak genuine paths and damage-supported molecules are preserved appropriately.

Generated data and reports remain ignored under `experiments/validation/results/clean_pilot/`. The combined generated table is `summary.tsv` in that directory.

## Simple-bubble validation

### Question

Does the non-destructive detector find simple alternatives reproducibly, and how do clean reads, terminal damage, sequencing errors, and low-abundance strain variation affect the detected topology?

### Setup

- Random 20,000 bp reference with fixed seed 901.
- Fixed 60 bp fragments at 20× total coverage.
- `k=21`, `min_count=2`, orientation-aware graph, then conservative tip cleaning.
- Three simulation seeds for clean, damage-only, sequencing-error, 2× rare-strain, and mixed conditions.
- Terminal C→T/G→A rates of 0.45, 0.30, 0.20, 0.12, and 0.06 from each end.
- Twelve spaced SNPs for the related-strain sensitivity test, repeated across five seeds and 1×, 2×, 5×, and 10× minor-strain coverage.

### Condition results

| Condition | Mean simple bubbles | Range |
|---|---:|---:|
| Clean | 0.00 | 0–0 |
| Terminal damage | 0.00 | 0–0 |
| Sequencing errors | 3.33 | 1–6 |
| Rare strain at 2× | 2.33 | 1–5 |
| Damage and rare strain | 2.00 | 0–5 |

Damage-only runs contained 3,704–3,862 recorded terminal events and 312–374 removed tips. Under these settings, terminal damage therefore appeared mainly as tips or incomplete branches rather than simple bubbles.

### Rare-strain sensitivity

| Minor-strain coverage | Mean bubbles detected from 12 introduced SNPs | Range |
|---|---:|---:|
| 1× | 1.0 | 0–3 |
| 2× | 2.6 | 1–5 |
| 5× | 8.4 | 6–11 |
| 10× | 11.6 | 11–12 |

Tip cleaning did not change the detected bubble count in any of the 20 strain-sensitivity runs. Very low-coverage alternatives were nevertheless often absent after `min_count=2`, showing that a fixed support threshold limits rare-strain recovery before topology analysis begins.

### Real-data smoke test

On the local 19,643-read FASTQ used for prior performance checks, the tip-cleaned graph contained 1,778 branching nodes and 10 detected simple bubbles. Detection took 0.03 seconds. Five control and five detection runs produced byte-identical FASTA outputs, confirming that reporting does not alter assembly.

### Interpretation

- The detector is deterministic, non-destructive, and inexpensive.
- Clean controls remained free of detected simple bubbles in this matrix.
- Sequencing errors and genuine strain variants can create the same simple topology, so coverage alone is not a sufficient classifier.
- Terminal damage must be evaluated in tips and incomplete branches as well as bubbles.
- The next scientific step is compact read-end and orientation evidence collected before tip cleaning, followed by conservative damage/error/variation scoring.

These are controlled development results, not publication-level sensitivity or specificity estimates. Exact event-level precision and recall require bubble path sequences and provenance-linked ground truth.

## Terminal-evidence ground-truth validation

### Question

Can read-end evidence recover graph edges created by terminal aDNA damage while distinguishing them from clean data and ordinary sequencing errors?

### Setup

- Ten random 20,000 bp references using seeds 901–910.
- Identical 60 bp, 20× source molecules for clean, damage, and error conditions within each seed.
- `k=21`, `min_count=2`, orientation-aware graph, and a five-base terminal window.
- C→T/G→A terminal rates of 0.45, 0.30, 0.20, 0.12, and 0.06.
- Independent sequencing-error rate of 0.005.
- Ground truth defined as changed k-mer edges that remained after `min_count=2`.

Run with:

```bash
PYTHONPATH=src:tests python3 \
  experiments/04_terminal_evidence_ground_truth.py
```

The ignored output directory contains `parameters.tsv`, `per_run.tsv`, and `summary.tsv`.

### Results

| Metric | Clean | Terminal damage | Sequencing error |
|---|---:|---:|---:|
| Introduced base events | 0 | 37,908 | 20,149 |
| Retained changed edges | 0 | 7,257 | 2,376 |
| Candidate edges | 0 | 5,460 | 1,753 |
| Ground-truth candidate edges | 0 | 5,460 | 1,501 |
| Terminal-only candidate edges | 0 | 5,460 | 82 |
| Tips | 0 | 3,583 | 146 |
| Simple bubbles | 0 | 0 | 12 |

Across retained damage-derived edges, terminal-only detection had 75.24% sensitivity. Against the clean and sequencing-error controls, pooled terminal-only precision was 98.52%. Damage reduced the mean longest unitig from 19,849 bp in clean graphs to 305 bp. All damage candidates were tips; none formed a simple bubble.

### Interpretation

- Topology alone is insufficient because sequencing errors also form candidate paths.
- Terminal localization strongly separates the simulated damage candidates from error candidates.
- Evaluation must use retained graph truth: most introduced base events do not survive support filtering as distinct graph edges.
- The current detector misses about one quarter of retained damage edges and does not yet cover incomplete branches.
- This validates the evidence representation under controlled conditions, not a complete damage-aware cleaning decision.

### Weak-tip-to-backbone matching extension

The same deterministic ten-seed matrix was extended to match each weak tip to an equal-length locally competing linear backbone. Matches report DNA identity, purine/pyrimidine (RY) identity, relative coverage, substitutions, and whether C→T/G→A changes have terminal evidence in the expected orientation.

The complete 139-test suite passed. A fresh run wrote `parameters.tsv`, `per_run.tsv`, and `summary.tsv` under `experiments/validation/generated/backbone_matching_review/` and reproduced the preceding user-run reports byte-for-byte.

| Metric | Clean | Terminal damage | Sequencing error |
|---|---:|---:|---:|
| Weak tips | 0 | 3,583 | 146 |
| Tips matched to a backbone | 0 | 3,538 | 146 |
| RY-exact matches | 0 | 3,538 | 47 |
| Damage-compatible matches | 0 | 3,538 | 9 |

Damage-tip match recall was 98.74%. Every matched damage tip was RY-exact and damage-compatible, but 9 of 146 random-error tips (6.16%) also passed the compatibility rule. Candidate recovery is therefore strong in this simulation, while the compatibility label is not specific enough to justify correction or graph removal. The next matrix must add reverse-orientation and 3′ G→A controls, multiple substitutions, incidental terminal evidence, alternative backbone topologies, low-abundance strain tips, and incomplete branches.

## Expanded weak-tip matching matrix

### Question

Does weak-tip matching remain sensitive across separate 5′ and 3′ damage profiles, and does the current compatibility rule protect ordinary errors and genuine rare-strain tips?

The matrix follows established ideas rather than treating one simulated damage profile as universal: position-specific terminal substitution profiles follow [mapDamage](https://academic.oup.com/bioinformatics/article/27/15/2153/404129), RY-space comparison follows the damage-tolerant strategy used by [CarpeDeam](https://doi.org/10.1186/s13059-025-03839-5), and rare-strain controls reflect the local-coverage ambiguity described by [metaSPAdes](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411777/).

### Setup

- Twenty random 20,000 bp references using seeds 1501–1520.
- Bidirectional 60 bp reads at 20× total coverage.
- `k=21`, `min_count=2`, and a five-base terminal window.
- Low, standard, and high damage profiles, each applied at the 5′ end only, 3′ end only, and both ends.
- Independent sequencing-error rates of 0.001, 0.005, and 0.010.
- A terminal-error control with ordinary substitutions concentrated in the same five-base end windows used as damage evidence.
- Related strains at 1×, 2×, 5×, and 10× coverage with 12 spaced SNPs.

Run with:

```bash
PYTHONPATH=src:tests python3 \
  experiments/06_tip_matching_validation.py
```

The ignored output directory `experiments/validation/generated/tip_matching_matrix/` contains `parameters.tsv`, `per_run.tsv`, `candidates.tsv`, and `summary.tsv`.

### Results

| Condition | Weak tips | Matched | Damage-compatible |
|---|---:|---:|---:|
| Clean | 0 | 0 | 0 |
| Terminal damage | 35,202 | 34,762 (98.75%) | 34,762 |
| Ordinary sequencing error | 1,081 | 1,058 | 31 |
| Terminally concentrated error | 583 | 583 | 46 (7.89%) |
| Genuine rare-strain tips | 110 | 110 | 5 (4.55%) |

Separate 5′ and 3′ profiles had comparable recovery, and explicit regression fixtures verified oriented 3′ G→A matching and a tip containing both C→T and G→A substitutions. Clean controls again produced no weak tips.

### Interpretation

- Equal-length local backbone recovery is stable across damage intensity, end, and read orientation.
- RY-exact, terminally compatible matching is highly sensitive for the simulated damage tips.
- Terminal localization alone is not specific: ordinary errors placed near molecule ends passed more often than uniformly distributed errors.
- Five genuine rare-strain tips also passed the current Boolean damage-compatibility rule, so automatic removal would cause measurable biological loss.
- The next production change must associate substitutions more directly with their supporting molecule/end evidence and compare damage, error, and variation scores. Graph modification remains deferred.
- Incomplete branches, indels, and non-linear competing backbones still require a later topology matrix.

### Evidence-scoring follow-up

The same 360-run matrix now evaluates the non-destructive matched-tip classifier. Each substitution is associated with its exact oriented tip edge. Candidate tables include terminal/internal observations, expected-end distance, terminal enrichment, strand balance, three heuristic scores, the selected label, and auditable reasons.

| Truth/control | Tips | Damage-like | Error-like | Variation-like | Ambiguous | Unmatched |
|---|---:|---:|---:|---:|---:|---:|
| Terminal damage | 35,202 | 30,960 | 0 | 0 | 3,802 | 440 |
| Ordinary sequencing error | 1,081 | 5 | 97 | 0 | 956 | 23 |
| Terminally concentrated error | 583 | 27 | 0 | 0 | 556 | 0 |
| Genuine rare-strain tips | 110 | 0 | 1 | 0 | 109 | 0 |

Among matched damage tips, damage-like recall was 89.06%. Against matched error and genuine rare-strain controls, 30,960 of 30,992 damage-like calls were true simulated damage (99.90%). The 32 false calls show that biochemical compatibility plus aggregate terminal support cannot prove that the same molecule carried the substitution. No genuine rare tip was called damage-like.

The error-like gate requires terminal depletion. It therefore identifies only a small, high-specificity subset and keeps most low-coverage alternatives ambiguous; this avoids relabelling rare biological paths as errors on coverage alone. Variation-like tip calls were withheld by the minimum-coverage gate in this matrix. Scores are heuristics rather than calibrated probabilities, and graph/FASTA output remains unchanged.

## Incomplete-branch validation

### Question

Can bounded weak branches outside the existing short-tip and simple-bubble representations recover additional damage-derived graph edges without labelling clean, error, or rare-strain controls as damage-like?

### Setup

- Twenty random 20,000 bp references using seeds 1501–1520.
- Bidirectional 60 bp reads at 20× coverage, `k=21`, `min_count=2`, and a five-base terminal window.
- Clean, standard two-ended damage, 1% sequencing-error, and 2× related-strain conditions.
- Branches bounded at 64 edges and at no more than 50% of matched-backbone support.
- Exact retained k-mer truth for damage and error; exact reference membership for rare-strain candidates.

Run with:

```bash
PYTHONPATH=src:tests python3 \
  experiments/07_incomplete_branch_validation.py
```

### Results

| Condition | Candidates | Truth candidates | Damage-like | Ambiguous |
|---|---:|---:|---:|---:|
| Clean | 0 | 0 | 0 | 0 |
| Terminal damage | 1,482 | 1,482 | 557 | 925 |
| Sequencing error | 72 | 72 | 0 | 67 |
| Genuine rare strain | 2 | 2 | 0 | 2 |

The five remaining error candidates were classified as two error-like and three variation-like. Existing tip and bubble candidates covered 10,584 of 14,194 retained damage-derived edges (74.57%). Incomplete branches contributed 2,394 additional truth edges, increasing combined coverage to 12,978 of 14,194 (91.43%). This recovered 66.32% of the damage truth previously missed by those two topology detectors.

The detector is intentionally broader than the tip cleaner: detection accepts support ratios through 0.50 because it does not modify the graph, while cleaning retains its 0.20 rule. The current validation supports non-destructive reporting only. Unequal-length alternatives, complex tangles, molecule-linked evidence, and additional damage/error/strain profiles remain future work.

## Empirical non-UDG/UDG pilot

Five non-overlapping 100,000-read R1 subsets were compared for untreated and full-UDG libraries from the same PES001.B dental-calculus sample. Reads shorter than `k=19` were excluded identically. Both treatments used `min_count=2`, orientation-aware construction, tip cleaning, bubble detection, and a five-base terminal window.

| Metric per million retained observations | Non-UDG | Full-UDG | Exact permutation p |
|---|---:|---:|---:|
| Bubbles | 11.73 | 5.66 | 0.0159 |
| Branching nodes | 391.04 | 269.15 | 0.0079 |
| Unitigs | 16,663.80 | 16,014.96 | 0.0079 |
| Tips | 8.65 | 7.51 | 0.4444 |
| Event terminal fraction | 9.5% | 7.6% | 0.3730 |

Non-UDG data had more normalized bubbles, branching nodes, and unitigs in all five replicates, supporting greater graph disturbance without UDG treatment. Tip and terminal-enrichment differences were not significant, and the libraries differ in layout and read-length composition. This pilot therefore supports a graph-level effect but does not independently prove that the observed difference is caused only by terminal damage.

## Compacted unitig-graph validation

### Question

Can the orientation-aware k-mer graph be compacted into a connected unitig-level graph without changing assembly sequences or losing coverage and terminal evidence?

### Validation

- All 113 unit tests passed.
- A 676,696-read clean bacterial benchmark was assembled with `k=31`, `min_count=2`, orientation-aware construction, and tip cleaning.
- The historical, initial compacted, and optimized compacted FASTA files had the same SHA-256 hash.
- Ten clean, ten terminal-damage, ten sequencing-error, and five rare-strain simulations were checked independently.
- All 35 controlled graphs preserved exact unitig order and sequences, total observations, terminal observations, and reverse-symmetric links.

### Large clean benchmark

| Metric | Historical extraction | Initial compaction | Optimized compaction |
|---|---:|---:|---:|
| Unitigs | 4,760 | 4,760 | 4,760 |
| Oriented unitig links | — | 12,986 | 12,986 |
| Largest contig | 46,465 bp | 46,465 bp | 46,465 bp |
| Raw FASTA N50 | 7,721 bp | 7,721 bp | 7,721 bp |
| Wall time | 117.52 s | 134.37 s | 130.64 s |
| Peak memory | 2,903,640 KiB | 2,937,744 KiB | 2,915,928 KiB |
| Compaction stage | — | 21.11 s | 6.82 s |

Constant-time oriented base extraction and removal of the temporary edge-owner array made the compaction stage 68% faster and reduced its measured peak-memory overhead by about 21 MiB relative to the initial implementation. The optimized complete run remained about 11% slower and 0.4% higher in peak memory than the historical run; read loading and k-mer graph construction also varied between runs, so these totals are directional rather than a controlled microbenchmark.

### Interpretation

- Scientific output and evidence accounting are unchanged.
- Compaction itself does not improve N50; it creates the connected assembly graph required for simplification.
- The current implementation briefly holds both graph representations, so construction-time peak memory does not decrease yet.
- Subsequent simplification can operate on 4,760 unitigs instead of roughly five million active k-mer edges.
- Generated validation outputs remain ignored under `experiments/validation/generated/unitig_graph/`.

## Unitig-level bubble-scoring smoke test

### Question

Can bounded alternatives be scored on the compacted graph without changing assembly output or adding meaningful runtime overhead?

The known two-path control contained one 10× path and one 3× path. The detector found the bubble and measured the weaker path at a local coverage ratio of 0.30. The complete test suite passed 119 tests, including unchanged graph state and byte-identical CLI output checks.

The local 19,643-read FASTQ was then assembled with `k=21`, `min_count=2`, orientation-aware construction, and tip cleaning:

```bash
anvaya assemble -i example.fastq.gz --k 21 --min-count 2 \
  --orientation-aware --clean-tips \
  --unitig-bubble-report unitig_bubbles.tsv \
  -o contigs.fasta
```

### Results

| Metric | Result |
|---|---:|
| Compacted unitigs | 3,637 |
| Unitig-level bubbles | 8 |
| Alternative paths | 8 |
| Alternatives with coverage ratio ≤0.30 and similarity ≥0.90 | 6 |
| Detection and report time | 0.03 s |
| Complete run time | 4.19 s |
| Complete-run peak memory | 111,344 KiB |
| FASTA compared with reporting disabled | Byte-identical |

### Interpretation

The detector is fast and non-destructive, and its raw scores separate several weak, highly similar alternatives from their locally strongest paths. The six threshold-matching alternatives are candidates only: this smoke test does not identify whether they are sequencing errors, damage-derived paths, or genuine low-abundance variation. N50 and contig lengths are unchanged by design.

The next step is a labelled validation matrix for sequencing errors, terminal damage, and rare strains. Bubble removal will be implemented only after that matrix establishes thresholds that improve contiguity without sacrificing accuracy or biological variation.

## Unitig-level bubble validation matrix

### Question

Can local coverage and sequence similarity identify removable sequencing-error paths without removing genuine low-abundance strain paths?

Five random 20,000 bp references were tested with matched clean, sequencing-error, terminal-damage, and related-strain conditions. Reads were 60 bp at 20× total coverage with `k=21`, `min_count=2`, orientation-aware construction, and conservative tip cleaning. Related-strain mixtures contained 15× major-strain and 5× minor-strain reads with 12 spaced SNPs.

```bash
PYTHONPATH=src:tests python3 \
  experiments/05_unitig_bubble_validation.py
```

### Detected alternatives

| Condition | Bubbles | Labelled weaker paths |
|---|---:|---:|
| Clean | 0 | 0 |
| Sequencing error | 17 | 17 error paths |
| Terminal damage | 0 | 0 |
| Rare strain | 40 | 39 rare paths and 1 major path |

Error-path relative coverage ranged from 0.123 to 0.378, while rare paths ranged from 0.175 to 0.827. Both classes had approximately 0.976 sequence similarity because the detected alternatives represented single-base differences.

### Threshold trade-off

The table reports selection among detected weaker paths. Selecting a rare path represents biological loss.

| Maximum relative coverage | Errors selected | Rare paths selected |
|---|---:|---:|
| 0.10 | 0/17 | 0/39 |
| 0.20 | 13/17 (76%) | 2/39 (5%) |
| 0.30 | 16/17 (94%) | 11/39 (28%) |
| 0.50 | 17/17 (100%) | 28/39 (72%) |

Similarity floors of 0.90 and 0.95 produced identical selections. A floor of 0.98 rejected every tested error and rare path, so similarity did not improve discrimination in this SNP-only matrix.

### Interpretation

- Coverage is useful but overlaps substantially between errors and low-abundance strains.
- A 0.20 cutoff is the least harmful tested trade-off, not a validated cleaning rule.
- Increasing the cutoff would probably merge more unitigs and improve N50, but rapidly increases genuine variant loss.
- Terminal damage again produced no simple bubbles, so damage-aware assembly also requires tip and incomplete-branch decisions.
- These percentages describe detected bubble paths, not all introduced sequencing errors or SNPs.

Automatic bubble removal remains deferred. The next scoring stage must incorporate terminal/internal and orientation evidence, then repeat this matrix before contiguity changes are accepted. Generated tables remain ignored under `experiments/validation/generated/unitig_bubbles/`.

## Extended combined-evidence validation

### Question

Does terminal/internal or strand evidence produce a bubble rule that removes detected sequencing errors without losing genuine strain paths across broader conditions?

The matrix was expanded to 20 random 20,000 bp references and bidirectional 60 bp reads. It included three independent sequencing-error rates (`0.001`, `0.005`, and `0.010`), low/standard/high terminal-damage profiles, and related strains at 1×, 2×, 5×, and 10× within 20× total coverage. Each related strain contained 12 spaced SNPs. In total, 220 graph conditions and 39,600 threshold combinations were evaluated.

The scorer retained orientation-correct forward/reverse support, left/right terminal observations, internal support, terminal fraction, terminal enrichment relative to the strongest local path, and strand balance. The run completed in 69.38 seconds with a 67,824 KiB peak for the experiment process.

### Detected paths

| Scenario | Introduced events/SNPs | Detected weaker paths |
|---|---:|---:|
| Clean | 0 | 0 |
| Error rate 0.001 | 7,987 | 2 error paths |
| Error rate 0.005 | 40,362 | 47 error paths |
| Error rate 0.010 | 80,448 | 133 error paths |
| Damage, all intensities | 222,262 | 0 bubble paths |
| Rare strain 1× | 240 SNPs | 11 biological paths |
| Rare strain 2× | 240 SNPs | 46 biological paths |
| Rare strain 5× | 240 SNPs | 166 biological paths |
| Rare strain 10× | 240 SNPs | 228 biological paths |

Detection is topology- and support-dependent, so these counts are not recall against every introduced base event. Rare-path recovery increased strongly with coverage. At 10×/10×, 117 weaker paths matched the nominal major strain and 111 matched the nominal rare strain, demonstrating that “weaker” is not a biological label when abundances are similar.

### Rule comparison

All rules required relative coverage no greater than 0.30 and sequence similarity at least 0.95. The terminal rules additionally required the weaker path's terminal fraction to be lower than the strongest path by the stated amount.

| Rule | Errors selected | Biological paths selected |
|---|---:|---:|
| Coverage only | 174/182 (95.6%) | 84/451 (18.6%) |
| Terminal enrichment ≤−0.05 | 159/182 (87.4%) | 39/451 (8.6%) |
| Terminal enrichment ≤−0.15 | 111/182 (61.0%) | 20/451 (4.4%) |

For the conservative `−0.15` rule, biological selection was 5/11 at 1×, 10/46 at 2×, 4/166 at 5×, and 1/228 at 10×. It therefore protects moderate/high-coverage strains much better than coverage alone but does not safely distinguish extremely low-coverage strains from errors. Strand-balance restrictions reduced error sensitivity without improving the best overall trade-off.

No tested rule selected any error path while maintaining zero biological-path loss across the full matrix. The earlier five-seed zero-loss result did not reproduce and must be treated as a pilot observation rather than a validated rule.

### Real-data performance check

On the existing 19,643-read FASTQ, baseline and evidence-enabled runs produced byte-identical FASTA files. The evidence report contained eight bubbles and 16 paths.

| Metric | Baseline median | Evidence median |
|---|---:|---:|
| Wall time, three paired runs | 4.49 s | 4.77 s |
| Peak memory | 111,436 KiB | 114,224 KiB |

Combined evidence added approximately 6% median runtime and 2.5% peak memory in this small smoke benchmark. The three-run timing sample is sufficient to detect a modest implementation cost, not to establish general performance.

### Interpretation

- Terminal enrichment materially improves error-versus-variation separation.
- Low-coverage genuine variants remain confounded with sequencing errors.
- Strand balance was not independently useful in this matrix.
- Damage again appeared as tips or incomplete branches rather than simple bubbles.
- Hard bubble removal remains unjustified; ambiguous paths should remain until stronger molecule, quality, linkage, or damage evidence is available.

Generated tables remain ignored under `experiments/validation/generated/unitig_bubbles_extended/`.

## Non-destructive unitig-path classification

### Question

Can the validated evidence be converted into explainable path labels without changing the assembly, and can damage-like labels avoid confusing genuine strain variation with terminal damage?

The 220-condition combined-evidence matrix was rerun after adding `error-like`, `damage-like`, `variation-like`, and `ambiguous` labels. The locally strongest path is reported separately as `dominant`. Each TSV decision includes machine-readable reasons, substitutions relative to the strongest path, and damage compatibility.

### Results

| Truth class | Error-like | Damage-like | Variation-like | Ambiguous |
|---|---:|---:|---:|---:|
| Sequencing error (182 paths) | 111 | 0 | 2 | 69 |
| Biological variation (451 paths) | 20 | 0 | 272 | 159 |

The error-like label had 84.7% precision among selected paths and retained 95.6% of detected biological paths from error classification. An initial substitution-plus-enrichment damage rule incorrectly labelled 37 biological paths as damage-like. Requiring low local coverage and strand skew removed those false calls in this matrix.

On the local 19,643-read FASTQ, classification reported eight bubbles and eight weaker alternatives, all ambiguous. Reporting took 0.03 seconds and produced a byte-identical FASTA compared with reporting disabled.

### Interpretation

- Classification is explainable, fast, and non-destructive.
- The error-like label is useful for prioritization but is not safe enough for automatic deletion.
- The absence of false damage-like calls supports conservative specificity.
- No damage-derived bubble was detected, so damage-label sensitivity is unmeasured.
- The next damage-aware implementation must classify weak tips and incomplete branches before tip cleaning.

Generated classification tables remain ignored under `experiments/validation/generated/unitig_path_classification/`.

## Molecule-linked substitution ablation

### Question

Does requiring multiple damage-compatible substitutions to co-occur on the same source read improve damage-tip classification without increasing error or rare-strain calls?

`08_molecule_linkage_validation.py` reuses the expanded weak-tip generators and classifies every matched candidate twice on the same graph: once with molecule linkage disabled and once enabled. The 20-seed matrix contains 460 paired conditions: clean reads; low, standard, and high 5′, 3′, and bidirectional damage; three uniform error rates; terminal-only errors; deliberately correlated C→T/G→A terminal errors; four ordinary rare-strain coverages; and four molecule-spanning C→T/G→A rare-strain coverages. Graph nodes, edges, and observations are asserted identical between aggregate and linkage-enabled builds.

### Results

| Truth class | Candidates | Multi-substitution | Aggregate damage-like | Linked damage-like | Label changes |
|---|---:|---:|---:|---:|---:|
| Damage | 34,762 | 1,155 | 30,960 | 31,018 | 58 |
| Sequencing error | 1,810 | 0 | 166 | 166 | 0 |
| Rare strain | 124 | 0 | 1 | 1 | 0 |

All 58 changes were true damage tips moving from ambiguous to damage-like. They occurred in standard and high damage conditions across 18 of 20 seeds; no label moved in the opposite direction. Candidate-level exact McNemar testing gives `p=6.94e-18`. At the independent-seed level, 18 seeds improved and none regressed, giving a two-sided sign-test value of approximately `7.63e-6` among discordant seeds.

For multi-substitution damage tips, recall increased from 1,093/1,155 (94.63%) to 1,151/1,155 (99.65%), a 5.02 percentage-point gain. Overall damage recall increased from 89.06% to 89.23%, so the material effect is concentrated where linkage supplies genuinely new information. Opt-in linkage increased graph-build time by 2.04× and used about 672 KB mean compact-array storage per 20 kb, 20× graph.

### Interpretation and limitation

- Read-level linkage provides a statistically supported sensitivity improvement for multi-substitution simulated damage tips.
- It caused no observed label regression in this matrix and does not change graph topology or assembly output.
- It does not identify biochemical cause. Correlated ordinary terminal errors can have the same observed C→T/G→A pattern: 134/169 paired-terminal-error tips were damage-like under both scoring modes.
- The error and rare-strain generators produced their paired changes as separate single-substitution tips. Consequently, they did not exercise the multi-substitution linkage gate. Linked-versus-unlinked multi-substitution non-damage paths remain a required adversarial fixture before any evidence-based removal rule.
- The next implementation step is a sample-level positional damage profile, followed by calibrated likelihoods that combine linkage with profile fit, base quality, coverage, and variation evidence.

Generated tables remain ignored under `experiments/validation/generated/molecule_linkage_validation/`.

## Candidate-locus damage-profile validation

### Question

Can molecule-linked matched alternatives recover the positional shape of a sample's terminal damage without a reference genome or supplied profile?

`09_damage_profile_validation.py` simulates bidirectional low, standard, and high damage on 20 independent 20,000 bp references. For every graph it matches weak tips, infers the profile using only graph and source-read evidence, and then compares the inferred strand-symmetric candidate-locus fractions with the known simulator profile. Absolute fractions are not expected to equal the generating probabilities because only graph-discordant matched loci enter the denominator.

### Results

| Damage profile | Mean observed loci | Raw correlation | Equal-locus correlation | Capped correlation | Endpoint enriched | Exactly monotonic |
|---|---:|---:|---:|---:|---:|---:|
| Low | 46.8 | 0.908 | 0.910 | 0.908 | 20/20 | 12/20 |
| Standard | 351.4 | 0.966 | 0.963 | 0.965 | 20/20 | 20/20 |
| High | 686.0 | 0.974 | 0.971 | 0.974 | 20/20 | 20/20 |

Every inferred profile had a higher terminal than interior endpoint. That older endpoint check was incorrectly described as proof of a decreasing curve. Schema version 4 checks every adjacent pair: low-damage profiles averaged 0.4 upward steps and only 12/20 were exactly monotonic; standard and high profiles were monotonic in all runs. Correlations ranged from 0.769–0.984 for low damage, 0.919–0.993 for standard damage, and 0.960–0.991 for high damage. Canonical bidirected traversal represented the simulated events predominantly in one reverse-complement-equivalent channel, so validation uses the combined C→T/5′ plus G→A/3′ curve while retaining both raw channels in the report.

Schema version 4 reran the same 60 graphs using exact changed-base cycles and separate endpoint/monotonicity metrics. Correlations were unchanged from schema version 3. Equal-locus weighting slightly improved the mean low-damage correlation but slightly reduced the standard and high correlations. Capped weighting remained within 0.001 of the raw mean in every condition. These alternatives are therefore retained as coverage-sensitivity diagnostics rather than replacing the raw estimator.

### Interpretation

- Matched source-read evidence reliably recovers the positional decay shape when enough eligible loci exist.
- Low-damage samples are noisier because they contribute far fewer graph-discordant loci.
- The output is an empirical candidate-locus profile, not a calibrated whole-library probability: unmatched reference opportunities and uncertainty intervals remain future work.
- The profile is non-destructive and is not yet used to alter classification scores or graph topology.

Generated tables remain ignored under `experiments/validation/generated/damage_profile_validation/`.

### Real-data diagnostic

The first paired-end run on 1,409,072 reads from `SRR32866683` produced
14,615 matched weak tips and 4,517 eligible C→T/G→A loci. The combined
candidate-locus fractions at distances 0–4 were 0.190, 0.114, 0.119, 0.030,
and 0.104. This shows terminal enrichment but not a smooth decay curve.

Damage-profile schema version 2 was used to investigate the unusually large
distance-3 reference denominator. Distances 1–4 had 243, 246, 242, and 245
contributing reference edge/channel pairs respectively, while their reference
observation counts were 856, 757, 2,894, and 757. The distance-3 excess is
therefore not repeated counting of many more graph loci. It is coverage
concentration among a comparable number of reference edges: the largest edge
contributed 782 observations (27.0% of the bin), compared with maximum shares
of 24.1%, 14.0%, and 11.1% at distances 1, 2, and 4.

This single sample supports treating the output as an observation-weighted
candidate-locus diagnostic, not yet as a calibrated library-wide damage
estimate. Robust inference will require per-locus normalization or hierarchical
weighting, uncertainty intervals, and validation on additional real libraries
with independently estimated damage.

Schema version 3 retained the 4,517 per-locus count vectors and compared the
three summaries on the same sample:

| Distance | Raw | Equal locus | Coverage capped at 20 |
|---:|---:|---:|---:|
| 0 | 0.190 | 0.291 | 0.201 |
| 1 | 0.114 | 0.244 | 0.153 |
| 2 | 0.119 | 0.231 | 0.146 |
| 3 | 0.030 | 0.200 | 0.104 |
| 4 | 0.104 | 0.217 | 0.115 |

The capped summary reduced the distance-3 versus distance-4 discontinuity from
0.075 to 0.011 while retaining terminal enrichment. Equal-locus weighting also
reduced the discontinuity but flattened the terminal contrast more strongly.
Neither alternative produced a strictly decreasing empirical curve, so the
next stage remains a fitted overdispersed decay model with uncertainty rather
than choosing a raw summary as the final biological estimate.

An independent public-data rerun completed in 319.72 seconds and reproduced
the same 9,456,876 graph nodes, 9,423,529 edges, 14,615 matched tips, and 4,517
eligible loci. Its event TSV and contig FASTA were byte-identical to the earlier
outputs. The new summaries therefore change only profile interpretation, as
intended; they do not change graph topology or assembly output.

Schema version 4 found that the earlier real-data profile used k-mer-boundary
distance rather than the exact changed-base cycle. With `k=21` and a five-cycle
window, all 4,517 structural loci still matched, but only 513 (11.36%) had any
changed-base observation inside cycles 0–4; 4,004 (88.64%) were boundary-only.
The corrected C→T/5′ counts were zero in this window. The remaining combined
profile is therefore the G→A/3′ channel: 0.219, 0.197, 0.179, 0.033, and 0.105.
Equal-locus values were 0.364, 0.321, 0.267, 0.235, and 0.205; capped values
were 0.225, 0.204, 0.184, 0.113, and 0.107. This strengthens the terminal
signal but still does not produce a monotonic empirical curve.

The corrected run preserved all topology and report counts: 9,456,876 nodes,
9,423,529 edges, 14,615 matched tips, and 740,785 unitigs. It completed in
385.95 seconds versus 319.72 seconds for the earlier run, while input loading
and graph construction were also slower on the host; one unpaired wall-time
comparison is insufficient to infer a performance regression.

## Candidate damage-likelihood validation

Damage-profile schema version 5 fits the per-locus counts with a beta-binomial
geometric mean curve, `background + amplitude * decay**distance`, and compares
it with a constant-background beta-binomial null. A fit requires at least 10
observed loci, three observed distances, three alternative observations, and
20 total observations. Parameters receive bounded conditional
likelihood-support intervals. The model remains candidate-conditioned and
report-only.

The expanded `09_damage_profile_validation.py` matrix produced:

| Damage profile | Fitted runs | Mean fitted-shape correlation | Mean LR statistic | LR range |
|---|---:|---:|---:|---:|
| Low | 20/20 | 0.9994 | 26.40 | 12.78–44.52 |
| Standard | 20/20 | 0.9988 | 178.73 | 147.54–203.14 |
| High | 20/20 | 0.9994 | 431.39 | 371.89–512.16 |

`10_candidate_damage_likelihood_validation.py` adds non-damage controls. All
20 clean runs had zero observed loci, and all 20 ordinary terminal-error runs
remained insufficient with a mean of 2.3 loci. Correlated C→T/G→A terminal
errors were deliberately adversarial: all 20 fitted, with a mean LR statistic
of 5.62 and range 0.00–12.97. The maximum slightly exceeded the minimum
low-damage statistic, so these results do not justify a fixed classification
threshold.

Fitting the previously validated `SRR32866683` schema-version-4 locus counts
produced an LR statistic of 47.99, decay 0.807, and fitted probabilities 0.336,
0.271, 0.219, 0.176, and 0.142. The background optimum was at the lower bound,
with conditional support extending to 0.020. The fit took about 14 seconds for
4,517 matched and 513 observed loci. Its smooth curve is evidence for
candidate-locus positional structure, not proof of damage or a whole-library
damage estimate.

Generated control tables remain ignored under
`experiments/validation/generated/candidate_damage_likelihood_validation/`.

## Quality-evidence validation

`11_quality_evidence_validation.py` checks the exact observation path from
FASTQ Phred scores through compact molecule links to matched-tip evidence. The
20-condition matrix combines one to five supporting molecules with Phred 8,
20, and 35 alternatives plus missing-quality input. For every quality-bearing
condition, the reported per-observation log likelihood exactly matched
`log(10^(-Q/10) / 3)`: -2.9407 at Q8, -5.7038 at Q20, and -9.1577 at Q35.
Across one to five observations the joint log likelihood scaled linearly, as
required by the current conditional-independence model. All missing-quality
conditions reported zero quality observations, one to five missing
observations, and no invented likelihood.

These results validate propagation and arithmetic, not biological
classification. Low-quality alternative calls are more compatible with the
sequencing-error explanation, but high-quality damage and genuine variation
remain indistinguishable from quality alone. The likelihood is therefore
report-only and must next be compared with held-out damage and variation
models. Generated tables remain ignored under
`experiments/validation/generated/quality_evidence_validation/`.

The quality-aware build was also run on the paired `SRR32866683` dataset. It
preserved the previous public-run summary exactly: 9,456,876 nodes, 9,423,529
edges, 21,279 tips, 82,711 incomplete branches, 3,703 bubbles, 14,615 matched
tips, and 740,785 unitigs. The 111,399 event rows contained 7,561 quality
observations across 2,504 matched candidate rows, with zero missing-quality
observations and a weighted mean Phred score of 38.01. The candidate damage
fit also remained unchanged at 513 observed loci, decay 0.807, and LR statistic
47.99. The contig file's SHA-256 was
`8dff90623e1f5159f3052c85d50a67ada4dc129db38a9cf964bfe485e0a01da1`;
the earlier contig file was no longer available for a direct byte comparison.

The run completed in 379.94 seconds. This is 60.22 seconds slower than the
fastest v3 run but 6.01 seconds faster than the earlier corrected public run;
graph construction and reporting varied together across those unpaired runs.
The available timings therefore do not establish a quality-storage regression.
A paired repeated timing comparison is required if performance becomes the
next decision gate.

## Held-out event-likelihood validation

`12_event_likelihood_validation.py` compares three explanations using the same
terminal molecule observations and qualities: a fixed positional damage curve,
the Phred sequencing-error model, and a fitted constant-frequency variation
model with a one-parameter BIC penalty. Across 20 seeds per condition, all 20
damage simulations ranked damage, all 20 independent sequencing-error
simulations ranked sequencing error, and all 20 constant-frequency simulations
ranked variation.

The adversarial condition generated high-quality errors with exactly the same
positional distribution as the damage curve. All 20 ranked damage. This is an
expected identifiability failure: once substitution, position, and quality have
the same distribution, these observations contain no variable that identifies
the generating mechanism. The output is consequently a report-only ranking
with an explicit uncalibrated reason, not a replacement classifier.

`13_cross_fit_profile_validation.py` applied the production fold assignment to
the existing `SRR32866683` schema-version-5 profile. All five held-out models
fitted successfully. Each held out 903–904 of 4,517 loci and trained on
3,613–3,614. Their fitted cycle-0 probabilities ranged from 0.3292 to 0.3396;
cycle-4 probabilities ranged from 0.1315 to 0.1531, and LR statistics ranged
from 33.34 to 44.42. The similar curves show that the public candidate profile
is not controlled by any one fold.

A full schema-version-6 rerun on the same 1,409,072 public reads preserved the
assembly exactly: `contigs-v4.fasta` and `contigs-v5.fasta` had the same
43,447,947-byte size and SHA-256 digest. All 45 pre-existing TSV fields matched
across all 111,399 rows, and the full candidate damage fit was unchanged. The
new likelihood layer scored 1,478 events (1,209 incomplete branches and 269
tips): 968 ranked damage, 24 sequencing error, and 486 variation. It explicitly
left 95,848 matched events insufficient; 92,988 had no eligible terminal
observations, while 1,834 lacked alternative support and 1,026 lacked reference
support.

Most rankings remain weak. Of the 1,478 scored events, 1,054 (71.31%) had a
margin below 1, 424 had a margin of at least 1, 154 at least 2, and 57 at least
5. None of the damage-ranked events had a margin of at least 5. Substantial
disagreement with the existing heuristic labels therefore reinforces that the
new ranking must remain report-only until calibration and systematic-error and
overdispersion controls are implemented.

The schema-version-6 run completed in 423.07 seconds versus 379.94 seconds for
the preceding run, an increase of 43.13 seconds (11.35%). Event reporting took
131.18 seconds versus 66.24 seconds, although graph construction was 23.65
seconds faster in the later run. The event TSV grew by 50.30% because of the 15
new likelihood and provenance fields. These unpaired measurements identify
report generation and output size as optimization targets but are not a stable
microbenchmark.

Generated tables remain ignored under
`experiments/validation/generated/event_likelihood_validation/` and
`experiments/validation/generated/cross_fit_profile_validation/`.

## Event-confidence calibration validation

`14_event_confidence_validation.py` fits a class-conditional conformal model on
5,000 independently seeded samples: 1,000 each for the fitted damage curve,
damage-profile misspecification, independent Phred errors, systematic errors
with the same positional curve as damage, and constant-frequency variation.
The calibration corpus is separate from 500 test samples. Nonconformity uses
the competing likelihood gap divided by the square root of molecule count;
damage-score envelopes propagate cross-fit profile instability.

At `alpha=0.01`, all 100 fitted-damage tests were protected as damage. Of 100
misspecified-damage tests, 99 were protected as damage and one as variation.
All 100 independent-error tests became `eligible_error`. All 100 systematic
damage-shaped errors were protected as damage because their observable
distribution remains indistinguishable from damage. All 100 variation tests
were protected, with 94 labelled variation protection and six conservatively
labelled damage protection. No damage, misspecified-damage, systematic-error,
or variation test became cleaning-eligible.

The resulting model was applied offline to the existing `SRR32866683` v5
event report. Of 111,399 rows, 1,478 had likelihood scores. The final gates
protected 340 as damage and 93 as variation, marked no event eligible for
cleaning, and left 110,966 rows insufficient. Before enforcing separate
molecule minima, eleven events passed the likelihood and multiple-testing
gates; all eleven had exactly one alternative molecule. The default minimum of
two alternative and five reference molecules correctly removed these
singletons from eligibility.

This experiment validates safe abstention for the included simulation family,
not a universal production calibration. Models must be retrained and held out
by sample for the relevant library preparation, damage range, community,
coverage, and error process. Generated files remain ignored under
`experiments/validation/generated/event_confidence_validation/`.

## Graph-derived event-confidence stress test

`15_graph_event_calibration_validation.py` moves calibration from generated
observation lists to actual orientation-aware graphs. It reuses exact mutated
k-mer truth for damage and sequencing errors and creates separated C→T/G→A
rare-strain SNPs that are deliberately compatible with the damage channels.
Matched tips, incomplete branches, and individual simple-bubble paths are
scored with a damage profile inferred from a separate high-damage graph for
the same simulated sample. Calibration and validation seed ranges do not
overlap.

The initial resolution-valid run used three 20 kb/20× calibration references
and two unseen validation references at exploratory `alpha=0.1`. Calibration
contained 2,762 scored damage, 24 error, and 12 variation events. The
calibrator now refuses to run when a class has too few examples to attain the
requested empirical p-value resolution; this prevents coarse calibration from
masquerading as conservative evidence.

Validation produced 2,805 damage-truth graph events, of which 1,909 were
scoreable; 150 error-truth events, of which 26 were scoreable; and 39
variation-truth events, of which six were scoreable. No event became
`eligible_error`. Thirteen damage events were protected as damage and eight as
variation; all error and variation events remained insufficient. Only two
independent-error events reached scoring and both had one alternative and one
reference molecule. The remaining scored errors were correlated terminal
errors with insufficient reference support. Weakening those gates would also
weaken the public-data singleton protection, so the correct next step is to
improve evidence recovery rather than tune the threshold.

The test also confirms a topology gap: damage primarily appears in tips and
incomplete branches, while most rare-strain controls appear as bubbles.

The follow-up implementation compares each equal-length non-reference bubble
path directly with its dominant competitor, scores ordinary substitutions
without inventing a damage channel, and counts paired mates by physical
fragment. Experiment 15 now writes `evidence_funnel.tsv` with retained,
detected, matched, and scored truth counts. Repeating the same three
calibration and two validation seeds increased scored calibration errors from
24 to 54 and validation errors from 26 to 47. Damage remained 2,762/1,909 and
variation 12/6 for calibration/validation, respectively. No recovered error
passed the conservative molecule gates; the added observations reveal an
event-ascertainment problem rather than support for weakening those gates.

On the 1,409,072-read `SRR32866683` run, v7 preserved all 111,399 event rows,
the byte-identical 43,447,947-byte assembly, and the byte-identical damage
profile. Scoreable rows increased from 1,478 in v5 to 25,415: 5,692 tips,
19,716 incomplete branches, and seven bubble alternatives. Direct comparison
matched 3,552 of 3,706 non-reference bubble paths. Final rankings were 24,453
variation, 852 damage, and 110 sequencing error; these remain uncalibrated,
report-only explanations.

The first expanded run required 244.46 seconds for event reporting. Profiling
showed repeated beta-binomial damage-model evaluations as the dominant cost.
Aggregating identical `(distance, alternative count, depth)` records and
identical `(allele, quality)` observations is mathematically equivalent and
reduced a controlled profile from 110.13 to 25.50 seconds with byte-identical
matrix outputs. Public event reporting then fell to 135.13 seconds and total
runtime from 567.21 to 446.85 seconds while preserving all 25,415 scores.

The event-conditioned follow-up evaluates all three explanations conditional
on observing at least one consistent alternative fragment and one consistent
reference fragment. On the same held-out matrix, damage-truth rankings improved
from 1,644 to 1,680 of 1,909 scoreable events. Twenty-one of 47 error events
ranked sequencing error instead of zero before conditioning, including 21 of
23 scoreable independent errors. The first version also incorrectly moved four
of six scoreable rare-strain events to sequencing error, demonstrating that a
high-quality singleton is not identifiable as error merely because it survived
event detection.

The refined policy therefore abstains when sequencing error ranks first but
the alternative consists of one fragment with Phred quality above 20. In the
held-out matrix, 17 Q20 independent errors remained error-ranked, four Q23
independent errors became ambiguous, and all four affected rare-strain events
became ambiguous rather than error-ranked. Damage rankings were unchanged at
1,680 damage and 229 variation. Systematic terminal errors remained confounded
with damage or variation, as expected from their deliberately damage-shaped
generating process. No validation event became cleaning-eligible.

On `SRR32866683`, v9 preserved the byte-identical 43,447,947-byte assembly,
byte-identical damage profile, all 111,399 event rows, and all 25,415 scoreable
events. Rankings were 646 damage, 470 sequencing error, 20,428 variation, and
3,871 ambiguous. Every ambiguous event was a high-quality singleton. Of the
470 error-ranked events, 420 were low-quality singletons, 39 had two alternative
fragments, and 11 had three or more. Event reporting fell from 254.33 seconds
in the initial conditioned v8 run to 158.15 seconds in v9; total runtime fell
from 585.55 to 464.59 seconds. These public results establish deterministic
non-regression and plausible candidate behavior, not biological ground truth.

The graph-context phase adds a reusable, report-only extractor for event type,
path length, substitutions, alternative/reference path observations, minimum
edge support, local relative coverage, sequence identity, damage compatibility,
branch degree, terminal fraction and enrichment, strand balance, damage
distance, joint-fragment support, and base quality. Experiment 15 now writes
these fields in `validation_events.tsv`, median values by truth and topology in
`feature_summary.tsv`, and pairwise probability-of-superiority measurements in
`feature_separation.tsv`.

The condition matrix was expanded from nine to 16 graphs per seed. It now uses
three independent error rates, independent terminal errors, systematic paired
terminal errors, three damage strengths, and ordinary plus damage-compatible
rare strains at 1×, 2×, 5×, and 10×. Rare mixtures are held at 20× total depth;
the earlier generator added rare coverage on top of 20× major coverage and
therefore confounded coverage comparisons.

Across three calibration and two unseen validation references, scored error
examples increased from 54/47 to 152/103 and variation examples from 12/6 to
26/17 for calibration/validation. Damage remained 2,762/1,909. Validation
rankings were 1,680 damage and 229 variation for damage truth; 19 error, 65
variation, 12 damage, and seven ambiguous for error truth; and 13 ambiguous,
three damage, one variation, and zero error for variation truth. No event became
cleaning-eligible.

Path length, alternative/reference observations, terminal fraction, terminal
enrichment, and local relative coverage strongly separated variation from the
current damage and error fixtures. Median variation path length was 21 edges
and median terminal fraction 0.33, versus one edge/1.0 for damage and two
edges/1.0 for error. Damage and error overlapped substantially, especially for
terminal error processes, so graph context is not yet a safe deletion model.
The evidence funnel detected and matched 96 validation variation events but
scored only 17, making exact evidence recovery beyond the five-base terminal
window the next implementation target.

The follow-up adds targeted whole-read evidence for ordinary substitutions.
After all event paths are matched once, only changed alternative/reference
edges are indexed during one additional read pass. Repeated occurrences are
collapsed to the best-quality call per physical molecule. Damage-compatible
C→T/G→A events retain the existing terminal-cycle model, so the additional
ordinary evidence cannot create an internal damage explanation.

Using the same three calibration and two validation references, scored errors
increased from 152/103 to 301/202 and variation from 26/17 to 78/53 for
calibration/validation. Validation error rankings were 49 error, 140 variation,
12 damage, and one ambiguous; variation rankings were 44 variation, six
ambiguous, and three damage. Damage results were unchanged at 1,909 scoreable:
1,680 damage and 229 variation. Calibration still produced no
`eligible_error` decisions, so the next bottleneck is calibration of the new
evidence scope rather than event recovery.

On `SRR32866683`, v12 preserved the byte-identical v10 contigs and damage
profile and all 111,399 event rows. Scoreable events increased from 25,415 to
78,663; insufficient-evidence rows fell from 75,463 to 22,215. Rankings changed
from 646 damage, 470 error, 20,428 variation, and 3,871 ambiguous to 645 damage,
2,994 error, 71,021 variation, and 4,003 ambiguous. Event reporting increased
from 171.65 to 547.07 seconds and total runtime from 504.61 to 853.12 seconds.
Peak RSS was 7,011,012 KiB with no swap. This report-only accuracy/performance
tradeoff is documented explicitly; assembly runs without event reporting do
not execute the targeted rescan.

Generated outputs remain ignored under
`experiments/validation/generated/graph_event_calibration_validation/`.

## Damage-aware tip-cleaning validation

### Question

Can the matched-tip classifier drive real graph simplification without deleting
simulated damage or rare-strain paths, and does the resulting assembly preserve
reference-backed accuracy?

The opt-in `--damage-aware-clean-tips` policy applies the existing `2 × k`
length and 20% competing-support limits for at most five rounds. It removes only
tips classified as `error-like` when their observed support is one-sided.
Damage-like, variation-like, ambiguous, unmatched, and bidirectionally
supported tips remain active.

The complete 20-seed experiment-06 matrix covered 360 runs and 36,977 detected
tip candidates. Applying the production gate to those truth-labelled rows gave:

| Truth condition | Candidates | One-sided error-like selected | Protected |
|---|---:|---:|---:|
| Damage | 35,202 | 0 | 35,202 |
| Sequencing error | 1,664 | 82 | 1,582 |
| Rare strain | 111 | 0 | 111 |

An error-like-only rule would have selected 97 error candidates and one
rare-strain candidate. The false rare-strain candidate and 15 error-like error
candidates had bidirectional support, so the final guard retained all 16. This
is a conservative selection audit rather than a claim that every protected
error is biologically real.

Actual graph mutation was then tested against the exact clean bacterial
reference on the existing 5×/20×, two-seed matrix. QUAST compared the frozen
baseline, support-only cleaning, and damage-aware cleaning on identical reads:

| Coverage | Seed | Tips removed | Baseline N50 | Damage-aware N50 | Genome fraction | Misassemblies |
|---|---:|---:|---:|---:|---:|---:|
| 5× | 42 | 0 | 695 bp | 695 bp | 73.523% → 73.523% | 0 → 0 |
| 5× | 43 | 0 | 693 bp | 693 bp | 73.784% → 73.784% | 0 → 0 |
| 20× | 42 | 10 | 4,648 bp | 4,665 bp | 92.947% → 92.952% | 0 → 0 |
| 20× | 43 | 7 | 4,680 bp | 4,703 bp | 92.920% → 92.921% | 0 → 0 |

Mismatches and indels per 100 kbp were unchanged in every condition. The new
policy is therefore high precision but low recall: it preserved measured
accuracy and made only 0.4–0.5% N50 gains at 20×, compared with the much larger
but less biologically guarded support-only cleanup.

The public `SRR32866683` v13 run removed 989 tips and 1,121 edges in 7.93
seconds. Unitigs fell from 740,785 to 739,660 and total assembled bases from
24,239,229 to 24,215,608. N50 remained 32 bp and the largest unitig remained
461 bp. The pre-cleaning event TSV and damage-profile JSON were byte-identical
to v12. Wall time changed from 13:08.86 to 13:28.28 and peak RSS from 7,011,012
to 7,033,856 KiB. Because this public library has no event truth, those results
measure execution, topology change, and performance—not biological accuracy.

The lack of public N50 improvement shows that terminal tips are not the main
fragmentation source in this graph. The next simplification target is the much
larger set of bounded incomplete branches, using the same protection-first
validation sequence before any broader bubble handling.

## K-mer size experiment

### Question

Can a larger fixed k-mer improve contiguity after weak k-mers are filtered?

The same reads and reference were assembled with `min_count=2` and `k` values 21, 31, 41, and 51. MEGAHIT and metaSPAdes outputs from the baseline run were included in the same QUAST comparison.

### Results

| Metric | k=21 | k=31 | k=41 | k=51 | MEGAHIT | metaSPAdes |
|---|---:|---:|---:|---:|---:|---:|
| Graph nodes | 9,902,440 | 9,923,379 | 9,897,805 | 9,827,816 | — | — |
| Branching nodes | 13,665 | 4,669 | 2,883 | 2,066 | — | — |
| Total unitigs | 24,336 | 10,708 | 10,610 | 14,162 | 106 | 196 |
| Contigs ≥200 bp | 8,332 | 4,135 | 5,674 | 8,576 | 106 | 103 |
| Largest contig | 10,197 bp | 28,204 bp | 14,337 bp | 7,836 bp | 448,867 bp | 449,487 bp |
| N50 | 1,698 bp | 3,804 bp | 2,665 bp | 1,667 bp | 151,608 bp | 130,616 bp |
| Genome fraction | 94.940% | 97.309% | 97.779% | 98.113% | 98.601% | 98.347% |
| Duplication ratio | 2.005 | 1.999 | 2.009 | 2.019 | 1.001 | 1.000 |
| Misassemblies | 0 | 0 | 0 | 0 | 2 | 0 |
| Mismatches per 100 kbp | 0.01 | 0.01 | 0.02 | 0.03 | 2.02 | 3.74 |
| Wall time | 3 min 15 s | 3 min 22 s | 3 min 21 s | 3 min 18 s | 38 s | 1 min 20 s |
| Peak memory | 5.83 GiB | 6.13 GiB | 6.41 GiB | 6.38 GiB | 299 MiB | 1.5 GiB |

### Why the results change with k

- `k=21` cannot distinguish enough repeated sequence, so the graph contains more branches and unitigs stop frequently.
- `k=31` adds enough sequence context to resolve many of those branches while retaining sufficient k-mer support, giving the best contiguity in this experiment.
- At `k=41` and `k=51`, longer k-mers become easier to break through sequencing errors or insufficient support. Genome fraction rises slightly, but contiguity declines.
- Graph node counts remain close to 10 million at every k, so k-mer choice does not address the core memory problem.
- Memory rises with longer string k-mers and remains far above the mature assemblers.
- Anvaya outputs unresolved unitigs, while MEGAHIT and metaSPAdes apply graph simplification, multi-k strategies, and additional path resolution to produce much longer contigs.

### Conclusions and next steps

For this clean 150 bp, 20× dataset, `k=31` and `min_count=2` provide the best current accuracy-contiguity trade-off. This is a dataset-specific result, not a universal default for short, low-coverage, damaged, or metagenomic reads.

The approximately 2.0 duplication ratio and graph size near twice the 5.08 Mb reference motivated orientation-aware construction. That representation and its compact storage have now been validated in the experiments above. The next development sequence is:

1. validate conservative tip cleaning on damaged and low-abundance data;
2. add independently tested bubble handling;
3. evaluate paired-read links and multi-k assembly;
4. validate the baseline on controlled mixtures and small communities;
5. introduce damage-aware evidence after the ordinary baseline is stable.

The full report remains ignored under `experiments/baseline/results/k_sweep/quast/`.

## Paired unitig-extension validation

The final paired-extension implementation was validated on the frozen
`GCF_001050915.2` 20× seed-42 ART read pair. Both comparison arms used `k=31`,
`min_count=2`, orientation-aware construction, a five-base end window, and
damage-aware tip cleaning. The treatment additionally enabled bounded reciprocal
paired extension. QUAST 5.3.0 evaluated contigs of at least 200 bp against the
known reference.

| Metric | Damage-aware baseline | 5-pair, 3:1 paired extension |
|---|---:|---:|
| Unitigs/contigs | 13,349 | 12,964 |
| Joined physical links | — | 387 |
| Largest contig | 24,470 bp | 24,470 bp |
| N50 | 4,665 bp | 4,733 bp |
| Genome fraction | 92.952% | 93.014% |
| Duplication ratio | 1.005 | 1.005 |
| Misassemblies | 0 | 0 |
| Mismatches per 100 kbp | 0.00 | 0.00 |
| Indels per 100 kbp | 0.00 | 0.00 |

The paired run mapped 461,916 of 512,820 physical pairs, evaluated 4,325
oriented junctions, work-limited 504, searched 751,210 graph states, and resolved
418 oriented junctions. Wall time was 7:23.96 with 6,775,572 KiB peak RSS,
compared with 3:51.48 and 6,243,052 KiB for the baseline. The continuity gain is
therefore modest and computationally expensive, but passed the accuracy gate.

Two broader policies were rejected. The 3-pair/2:1 rule produced N50 5,296 bp
but 0.43 mismatches and 0.04 indels per 100 kbp. Carrying evidence through a
unique predecessor path produced N50 5,297 bp but one misassembly, 1.50
mismatches, and 0.04 indels per 100 kbp. The production implementation retains
direct junction anchors, strict support, reciprocal agreement, and work-capped
abstention.

On public `SRR32866683` data, the bounded resolver completed in 34.22 seconds
after the earlier unbounded traversal failed to finish after more than 50
minutes. It accepted only 17 `k=21` joins and left N50 at 32 bp. A fixed `k=31`
run raised N50 to 38 bp and reduced branching nodes from 247,669 to 194,244, but
reduced the longest contig from 461 to 366 bp. These mixed fixed-k results make
iterative multi-k construction the next continuity target; paired extension
remains optional and is not presented as the solution to public-data
fragmentation.

## Direct read-thread extension validation

Experiment 16 audits direct source-read crossings of compacted-graph junctions.
It counts independent physical molecules, requires a five-molecule/3:1 strong
winner, and accepts a physical link only when the reverse orientation selects
the reciprocal transition. A boundary-only edge index avoids retaining or
rescanning mappings for every graph edge. Report-only execution is
non-destructive and produced byte-identical FASTA output in validation.

Five seeds were tested under clean, standard bidirectional terminal damage,
0.5% sequencing error, and 15× major plus 5× related-strain conditions. Across
20 runs, all 2,335 reciprocal links were locally present in the known reference
or related-strain truth. Final-contig truth was also measured because locally
valid links can compose through repeated unitigs.

On the clean 20× `GCF_000007145.1` ART dataset, fixed `k=31`, damage-aware tip
cleaning produced 6,342 unitigs. Direct threading initially accepted 2,401
links and produced 3,941 contigs:

| Metric | Damage-aware baseline | Unguarded threading | Repeat-guarded threading |
|---|---:|---:|---:|
| Unitigs/contigs | 6,342 | 3,941 | 3,942 |
| Largest contig | 20,336 bp | 81,176 bp | 81,176 bp |
| N50 | 4,856 bp | 18,021 bp | 18,021 bp |
| Genome fraction | 96.903% | 97.262% | 97.260% |
| Duplication ratio | 1.005 | 1.001 | 1.001 |
| Misassemblies | 0 | 0 | 0 |
| Mismatches per 100 kbp | 0.00 | 0.00 | 0.00 |
| Indels per 100 kbp | 0.00 | 0.02 | 0.00 |

The unguarded indel was a seven-base discrepancy inside a repeated 2,530 bp
unitig, not at an overlap boundary. Its adjacent 37 bp hub had mean edge support
39.43, compared with 15.22 and 19.63 on its flanks. A same-read
adjacent-junction gate was rejected because it removed 19 joins, reduced N50 to
17,474 bp, and did not remove the indel. The final local copy-number guard
requires the internal unitig to have at least 2× the support of both flanks and
then discards the weaker adjacent join. It removed one link and eliminated the
indel without reducing N50 or the largest contig.

Experiment 17 adds short variable fragments, terminal damage, sequencing error,
a related strain, and unrelated contamination. Across 20 runs it observed
6,989 reciprocal links. Exact-reference truth counted 6,987 as correct; the two
failures occurred in one damaged strain-plus-contamination run. All 6,989 were
correct under R/Y topology truth, which treats C→T and G→A damage as compatible,
and no condition gained false contigs. The experiment reports exact and R/Y
truth separately rather than treating retained damage as an assembly error.

Public `SRR32866683` results at `k=31` were:

| Metric | Fixed-k baseline | Unguarded v20 | Guarded v21 |
|---|---:|---:|---:|
| Unitigs/contigs | 639,369 | 606,700 | 606,706 |
| Joined links | — | 32,669 | 32,663 |
| Largest contig | 366 bp | 435 bp | 435 bp |
| N50 | 38 bp | 40 bp | 40 bp |
| Peak RSS | 6,627,640 KiB | 6,619,336 KiB | 6,601,744 KiB |

The guard rejected six links, or 0.018% of reciprocal candidates, and retained
the continuity gain. The measured v21 wall time was 6:55.61 versus 5:38.61 for
v20, but graph construction alone varied by 35.58 seconds and v21 wrote a
larger coverage-annotated audit report. Repeated support-object allocation found
during this comparison was removed afterward; the optimized arithmetic path
has unit-test coverage but was not re-timed on the public dataset. No biological
accuracy claim is made for this unlabeled library.

## Public-data fragmentation attribution

The report-only fragmentation diagnostic was run on the paired `SRR32866683`
fixed-`k=31` graph before topology-changing extension. It classified both
physical ends of all 639,936 unitigs and left the assembled FASTA unchanged.

| Unitig topology | Count | Fraction |
|---|---:|---:|
| Isolated | 135,162 | 21.1% |
| One-sided | 295,792 | 46.2% |
| Connected | 208,982 | 32.7% |

| Physical end category | Count | Fraction |
|---|---:|---:|
| Dead end | 566,116 | 44.2% |
| Unique continuation | 501,026 | 39.1% |
| Ambiguous branch | 212,730 | 16.6% |

This establishes a reference-free baseline for where continuity stops. Nearly
half of all unitigs are one-sided, while isolated unitigs account for more than
one fifth; ambiguous branches are important but are not the only fragmentation
source. These counts describe graph topology, not biological correctness, and
must not be used alone to justify joining or deleting paths. The generated TSV
and run outputs remain excluded from version control.

## P1 dead-end attribution and P2 correction audit

The P1 second pass was added to determine why retained unitig paths stop rather
than inferring a cause from topology alone. Across the frozen 30-run ladder,
uniquely filtered extensions were the largest attributable class in every
perturbed scenario:

| Scenario | Uniquely filtered | Dead ends | Additional P2 candidates |
|---|---:|---:|---:|
| Damage | 64 | 97 | 0 |
| Damage plus error | 101 | 135 | 1,677 |
| Sequencing error | 42 | 48 | 1,868 |
| Rare strain | 146 | 203 | — |
| Contamination | 303 | 339 | — |

Damage-only reads generated no low-quality P2 candidates because the simulated
postmortem substitutions deliberately retained high Phred scores. This is a
useful negative result: ordinary quality-based correction cannot represent the
damage model.

The corrected P1 public-data run used all 1,409,072 `SRR32866683` reads at
`k=31`, `min_count=2`. It completed in 235.25 seconds, including 69.68 seconds
for attribution:

| P1 cause | Count | Fraction of 566,116 dead ends |
|---|---:|---:|
| Read boundary | 85,190 | 15.0% |
| Coverage gap | 0 | 0.0% |
| Uniquely filtered | 254,761 | 45.0% |
| Filtered conflict | 226,159 | 40.0% |
| Retained context elsewhere | 6 | <0.1% |

Within the uniquely filtered class, 8,299 dead ends were predominantly
low-quality and 20,902 were terminal-only. Thus approximately 85% of dead ends
had an observed continuation removed by `min_count`, but almost half of all
dead ends had conflicting raw alternatives and cannot be restored safely from
abundance alone.

The public P2 stage built its spectrum from all reads and audited 20,000 evenly
spaced reads in 92.64 seconds. It reported 2,926 unique complete solid-k-mer
rescues, 50 protected terminal damage-compatible reversals, 58 ambiguous
rescues, and 10,863 incomplete rescues. Its combined invocation later failed
inside the then-unfixed P1 ambiguous-base path, so these are valid P2-stage
counts rather than a completed combined-run benchmark.

Neither stage changed N50 or contig sequences. Their positive result is
localization of the bottleneck: uniquely filtered continuations are a large
candidate pool, while Phred-supported candidates form only a small subset and
filtered conflicts require additional context such as multi-k or molecule-path
evidence.

## Oracle filtered-edge rescue experiment

Experiment 18 now includes two truth-only upper-bound comparators. The
`oracle_unique` arm starts from retained compacted-graph dead ends and rescans
both orientations of the source reads. It selects a continuation only when
exactly one base was observed below `min_count` and its canonical k-mer occurs
in a known generating reference. The broader `oracle_all_truth` arm selects
every observed sub-threshold canonical k-mer present in generating truth.
Both supply exactly enough additional observations for selected edges to cross
the existing threshold. Neither corrects reads, removes damage/error k-mers, or
rescues an unobserved edge.

Each oracle reruns graph construction and reports baseline-versus-oracle
unitigs, total bases, largest contig, N50, primary-reference recovery,
rare-strain recovery, false contigs, and explicit deltas. A truth-valid edge can
still participate in an incorrect graph path, so false-contig changes remain a
required safety outcome rather than being assumed to be zero.

Together the comparators separate the P1 unique-dead-end hypothesis from the
general recovery of weak low-abundance sequence. A small `oracle_unique` gain
is a no-go signal for production dead-end rescue even when `oracle_all_truth`
improves recovery elsewhere. A large, accuracy-preserving unique-arm gain would
justify developing a conservative classifier using evidence available without
a reference.

The completed 30-run ladder produced the following `oracle_unique` means:

| Scenario | Rescued k-mers | Baseline N50 | Oracle N50 | N50 delta |
|---|---:|---:|---:|---:|
| Clean | 1.8 | 4,979.8 | 4,981.6 | +1.8 |
| Damage | 1.8 | 472.4 | 472.6 | +0.2 |
| Sequencing error | 1.8 | 1,314.6 | 1,314.6 | 0.0 |
| Damage plus error | 1.8 | 291.6 | 291.6 | 0.0 |
| Rare strain | 14.6 | 68.8 | 68.6 | -0.2 |
| Contamination | 49.2 | 377.4 | 389.6 | +12.2 |

Across all scenarios, 355 truth-valid unique continuations were rescued. N50
was unchanged in 19 runs, increased in ten, and decreased in one. The
false-contig count did not change, but exact primary-reference recovery fell by
1,230 bp in aggregate. One sequencing-error run lost 1,103 recovered bases
after two individually truth-valid edges connected into a path that was no
longer exactly reference-contained. Edge-level truth therefore did not imply
contig-level correctness.

This is a no-go result for production singleton dead-end rescue. It provided
effectively no continuity improvement in damage or sequencing-error scenarios;
the broader oracle's earlier contamination gain primarily measured recovery of
weak sequence across the separate low-coverage genome. The next upper-bound
experiment should independently normalize known simulated damage and sequencing
errors before deciding between damage-aware consensus projection and multi-k
construction.

## Damage-aware read consensus

Experiment 18 now includes an opt-in pre-graph consensus comparator derived
from the conservative overlap principles used by CarpeDeam and local
quality-aware clustering used by BayesHammer. Exact internal anchors propose
read overlaps; full overlaps must pass 90% DNA and 99% R/Y identity. Only
high-quality nonterminal observations from independent molecules vote on 5′
T-to-C or 3′ A-to-G reversals. Three votes with strict greater-than-4:1
dominance are required, and tied or repetitive placements abstain.

The ladder records corrected bases, true and false corrections, correction
precision, unitigs, largest contig, N50, primary and rare-strain recovery,
false contigs, and baseline deltas. The initial acceptance gate is at least 99%
correction precision, at least 25% N50 improvement in both damage and combined
damage/error scenarios, no clean regression or false-contig increase, no
primary-reference loss, and less than 1% rare-strain recovery loss. Failure of
any accuracy gate rejects topology use; weak continuity gain closes this design
and advances multi-k construction instead of threshold tuning.

The first exhaustive-anchor implementation passed the clean and precision
gates and produced large mean N50 gains: 472.4 to 4,979.8 bp for damage and
291.6 to 1,314.6 bp for combined damage and sequencing error. It nevertheless
failed the strain-safety gate: exact rare-strain recovery fell 9.11%. On the
1,409,072-read `SRR32866683` dataset, its all-read, both-orientation Python
index exceeded the 7.6 GiB WSL memory limit and was killed at 7.17 GiB
anonymous RSS before graph construction.

The scalable revision retains at most eight bottom-ranked canonical internal
anchors per read, packs occurrences into integers, reconstructs reverse
orientation only for shortlisted overlaps, reuses unchanged read objects, and
streams audit rows instead of retaining them. The optimized 30-run ladder was:

| Scenario | Baseline mean N50 | Exhaustive consensus | Bounded consensus |
|---|---:|---:|---:|
| Clean | 4,979.8 | 4,979.8 | 4,979.8 |
| Damage | 472.4 | 4,979.8 | 4,803.6 |
| Damage plus error | 291.6 | 1,314.6 | 982.0 |
| Contamination | 377.4 | 683.0 | 683.0 |
| Rare strain | 68.8 | 93.2 | 93.2 |
| Sequencing error | 1,314.6 | 1,314.6 | 1,314.6 |

The bounded version retained 100% pooled correction precision outside the
rare-strain scenario and 99.74% within it, but rare-strain recovery still fell
8.98%. It therefore remains unsafe for strain mixtures despite passing the
continuity threshold in the damage scenarios.

The optimized public-data run completed in 511.35 seconds with 6,690,384 KiB
peak RSS. Consensus occupied 326.39 seconds, evaluated 4,168,135 terminal
candidates, corrected 157 bases in 150 reads, and wrote a 476,025,873-byte TSV.
Compared with the fixed-`k=31` baseline, contigs decreased from 639,936 to
639,886 and total assembled bases from 27,517,907 to 27,516,329; largest contig
remained 366 bp and N50 remained 38 bp. This closes damage correction alone as
the public-data continuity mechanism. The reusable overlap and R/Y validation
remain opt-in research infrastructure; the next architecture experiment must
test iterative whole-fragment cluster-and-extend assembly rather than feeding
corrected reads back into the same fixed-k graph.
