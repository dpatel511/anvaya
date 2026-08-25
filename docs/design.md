# Initial design

Anvaya is a reference-free, damage-aware de Bruijn graph assembler for ancient metagenomic short reads.

## Initial scope

- bacterial and archaeal metagenomes;
- Illumina short reads;
- untreated double-stranded libraries;
- single-end or merged fragments initially;
- fixed k-mer size initially;
- no reference genome or damage profile required as input;
- Python research implementation before optimization.

## Core idea

Anvaya will construct one graph containing competing paths. It will preserve evidence from original molecules instead of reducing each k-mer immediately to a count.

Initial evidence includes molecule count, distance from both molecule ends, read orientation, base quality, and read identity. Damage-compatible terminal support will be evaluated differently from internal, orientation-balanced support. Ambiguous paths will remain.

## Intended workflow

1. Read and preprocess sequences.
2. Extract k-mers with provenance.
3. Construct the evidence-annotated graph.
4. Compact non-branching paths.
5. Identify tips, bubbles, and competing paths.
6. Compare damage and genuine-variation explanations.
7. Simplify conservatively and output contigs.

Internal damage calibration will be part of the assembler. PyDamage and similar tools are comparators, not required components.

## Orientation-aware graph

The experimental graph stores each canonical `(k-1)`-mer once and uses oriented integer handles to traverse either strand. A physical k-mer edge therefore represents both its forward traversal and the reverse-complement traversal. Forward, reverse, and palindromic observations are retained separately, while `min_count` is applied to their combined support.

This representation prevents the two DNA strands from producing duplicate assemblies and establishes the structure needed for later orientation-dependent damage evidence. It does not yet store read-end position, base quality, molecule identity, or library-level damage evidence.
## Evidence from the current baseline

Minimum-support filtering showed that unsupported k-mers cause most initial graph fragmentation. On the clean single-genome benchmark, `min_count=2` increased genome recovery from 21.444% to 94.940% without introducing a reported misassembly.

A fixed-k sweep showed that `k=31` gives the best current contiguity for 150 bp reads at 20× coverage. Larger k values recover slightly more of the reference but fragment valid paths. This result is dataset-specific and does not establish a default for ancient metagenomes.

All filtered assemblies have a duplication ratio near 2.0, and their graph node counts are close to twice the reference length. This is evidence that forward and reverse-complement paths are currently represented separately. The current output is also a set of unitigs that stops at unresolved branches, rather than fully resolved contigs comparable to mature assemblers.

## Near-term development sequence

1. Design an orientation-aware graph that preserves strand evidence; naïve canonicalization must not discard direction needed for assembly or later damage inference.
2. Represent k-mers compactly as encoded integers and replace nested Python objects with flatter graph structures.
3. Add conservative tip and bubble handling as optional, independently benchmarked operations.
4. Evaluate paired-read links and a controlled multi-k strategy for resolving remaining branches.
5. Freeze and validate the ordinary baseline across clean isolates, controlled strain mixtures, and small communities.
6. Add terminal position, orientation, quality, and molecule evidence for damage-aware branch decisions.

Every algorithmic change will be compared with the frozen baseline for genome recovery, contiguity, misassemblies, mismatch rate, low-abundance retention, runtime, and peak memory.
