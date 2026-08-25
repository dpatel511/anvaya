# Baseline accuracy experiment

This experiment checks the unmodified assembly baseline against a known reference before damage-aware graph decisions are introduced.

## Requirements

- ART (`art_illumina`)
- MEGAHIT (`megahit`)
- SPAdes (`spades.py`)
- QUAST (`quast.py`)
- `/usr/bin/time`

## Preview

```bash
python3 experiments/02_baseline_accuracy.py \
  --reference reference.fasta \
  --output-dir experiments/baseline/results/run_01 \
  --dry-run
```

Remove `--dry-run` to execute the benchmark. The default simulation uses paired 150 bp HiSeq 2500 reads, 20× coverage, a 300 bp mean fragment length, and seed 42.

Anvaya currently pools both mate files into the graph but does not use pairing links. This limitation must be retained when interpreting comparisons with MEGAHIT and metaSPAdes. QUAST evaluates contigs of at least 200 bp against the supplied reference.

Generated reads, assembler outputs, logs, and QUAST results must remain inside an ignored `results/` directory and must not be committed.
