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

The compact implementation represents canonical k-mers with rolling 2-bit integers. Physical edge endpoints and strand-support counts are stored in typed arrays, while node degrees and traversal state use byte arrays. This replaces per-node counters, string k-mers, tuple edge keys, and per-edge support objects without changing graph decisions.

## Conservative tip cleaning

The experimental cleaner considers only dead-end paths no longer than `2 × k`. A path is removed when its strongest edge has at most 20% of the support of a competing edge at the attachment branch. Long tips, isolated paths, and strongly supported alternatives are preserved.

Edges are deactivated directly in the compact adjacency arrays, and degrees and active graph totals are updated in place. This avoids copying the graph or allocating per-edge cleaning objects. The operation remains opt-in until it is validated on damaged and low-abundance data.

## Evidence from the current baseline

Minimum-support filtering showed that unsupported k-mers cause most initial graph fragmentation. On the clean single-genome benchmark, `min_count=2` increased genome recovery from 21.444% to 94.940% without introducing a reported misassembly.

A fixed-k sweep showed that `k=31` gives the best current contiguity for 150 bp reads at 20× coverage. Larger k values recover slightly more of the reference but fragment valid paths. This result is dataset-specific and does not establish a default for ancient metagenomes.

The directed filtered assemblies have a duplication ratio near 2.0. The orientation-aware graph corrected this to 1.005 and reduced physical nodes from 9.92 million to 5.01 million. Its compact representation preserved the same unitigs while reducing peak memory from 6.92 GiB to 2.78 GiB and wall time from 193.90 seconds to 113.99 seconds.

Across two clean bacterial genomes, two seeds, and 5×/20× coverage, conservative tip cleaning introduced no reported misassemblies. At 20×, N50 improved by 56–67% and the largest contig by 52–128%. At 5×, cleaning made few changes and did not reduce genome fraction. This supports the current rule on clean isolates but does not establish safety for damaged molecules or rare strains.

## Near-term development sequence

1. Validate tip cleaning on damaged reads and low-abundance related strains.
2. Add conservative bubble handling as an optional, independently benchmarked operation.
3. Evaluate paired-read links and a controlled multi-k strategy for resolving remaining branches.
4. Freeze and validate the ordinary baseline across clean isolates, controlled strain mixtures, and small communities.
5. Add terminal position, orientation, quality, and molecule evidence for damage-aware branch decisions.
6. Profile stable hot loops and consider Cython only where pure Python remains a bottleneck.

Every algorithmic change will be compared with the frozen baseline for genome recovery, contiguity, misassemblies, mismatch rate, low-abundance retention, runtime, and peak memory.
