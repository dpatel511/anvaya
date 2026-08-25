# Anvaya

Anvaya is a research project exploring damage-aware de Bruijn graph assembly for ancient bacterial and archaeal metagenomes.

The central idea is to retain evidence from the original DNA molecules—such as position within the fragment, orientation, and base quality—and use it to distinguish postmortem-damage-induced graph paths from genuine biological variation.

## Current status

Early Python prototype. Anvaya supports sequence-file input, basic de Bruijn graph assembly, topology analysis, unitig extraction, and controlled graph-disturbance tests.

Correctness and scientific validation are the priorities before speed or memory optimization.

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
- node degrees and source, sink, and branch identification;
- maximal non-branching unitig extraction;
- graph topology and support summaries;
- end-to-end file-to-unitigs assembly;
- deterministic test-only simulation helpers;
- unit tests for assembler primitives.

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
