# Anvaya

Anvaya is a research project exploring damage-aware de Bruijn graph assembly for ancient bacterial and archaeal metagenomes.

The central idea is to retain evidence from the original DNA molecules—such as position within the fragment, orientation, and base quality—and use it to distinguish postmortem-damage-induced graph paths from genuine biological variation.

## Current status

Early Python prototype. Anvaya currently supports validated k-mer extraction, basic de Bruijn graph construction, topology analysis, unitig extraction, and an end-to-end reads-to-unitigs workflow.

Correctness and scientific validation are the priorities before speed or memory optimization.

## Initial scope

- ancient bacterial and archaeal metagenomes;
- Illumina short reads;
- untreated double-stranded libraries initially;
- reference-free assembly;
- no damage profile required as input.

## Implemented

- overlapping k-mer extraction;
- DNA sequence and k-mer size validation;
- directed de Bruijn graph construction;
- repeated-edge counting and branch representation;
- incoming and outgoing node degrees;
- source, sink, and branching-node identification;
- maximal non-branching unitig extraction;
- graph topology and support summaries;
- end-to-end reads-to-unitigs assembly;
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
