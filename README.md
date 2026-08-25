# Anvaya

Anvaya is a research project exploring damage-aware de Bruijn graph assembly for ancient bacterial and archaeal metagenomes.

The central idea is to retain evidence from the original DNA molecules—such as position within the fragment, orientation, and base quality—and use it to distinguish postmortem-damage-induced graph paths from genuine biological variation.

## Current status

Early Python prototype. Anvaya supports sequence-file input, directed and experimental orientation-aware de Bruijn graph assembly, topology analysis, unitig extraction, and controlled graph-disturbance tests.

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
- linear-time node degree calculation;
- source, sink, and branch identification;
- maximal non-branching unitig extraction;
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

## Testing

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

## Documentation

- [Research question](docs/research_question.md)
- [Initial design](docs/design.md)
- [Benchmark plan](docs/benchmark_plan.md)
- [Literature notes](notes/literature.md)
- [Simulation experiments](experiments/README.md)
