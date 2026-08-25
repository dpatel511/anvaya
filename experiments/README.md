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

The approximately 2.0 duplication ratio after filtering is consistent with forward and reverse-complement paths being assembled separately. Orientation-aware graph construction is therefore the next design priority. The remaining runtime and memory gap also shows that compact k-mer and graph representations will still be required.

The full comparison remains ignored under `experiments/baseline/results/min_count_comparison/quast/`.
