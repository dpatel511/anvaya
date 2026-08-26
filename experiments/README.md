# Experiments

Generated FASTQ files and result directories must not be committed. Publication-level experiments will use established simulators with recorded versions, commands, configurations, seeds, and reference checksums.

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
