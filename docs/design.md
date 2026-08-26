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

This representation prevents the two DNA strands from producing duplicate assemblies. With a positive end window, each physical edge now stores compact distance-binned support from both canonical ends, internal support, and ambiguous palindromic terminal support. Reverse observations swap their end distances into canonical orientation.

Base quality, molecule identity, paired-fragment links, and an inferred library-level damage model are not yet stored.

The compact implementation represents canonical k-mers with rolling 2-bit integers. Physical edge endpoints and strand-support counts are stored in typed arrays, while node degrees and traversal state use byte arrays. This replaces per-node counters, string k-mers, tuple edge keys, and per-edge support objects without changing graph decisions.

## Compacted unitig graph

After k-mer-level cleaning and event reporting, each active maximal non-branching path is represented as one physical unitig. Oriented unitig handles retain links in both directions without duplicating reverse-complement sequences. Each unitig stores its sequence, path length, total, minimum, and mean edge support, plus aggregated terminal and internal observations.

The CLI now writes sequences from this representation. No simplification decision is made at the unitig level yet, so compaction must preserve the previous FASTA exactly. The connected unitig graph is the intended substrate for local-coverage bubble and branch decisions because those operations can inspect thousands of unitigs rather than millions of k-mer edges.

## Conservative tip cleaning

The experimental cleaner considers only dead-end paths no longer than `2 × k`. A path is removed when its strongest edge has at most 20% of the support of a competing edge at the attachment branch. Long tips, isolated paths, and strongly supported alternatives are preserved.

Edges are deactivated directly in the compact adjacency arrays, and degrees and active graph totals are updated in place. This avoids copying the graph or allocating per-edge cleaning objects. The operation remains opt-in until it is validated on damaged and low-abundance data.

## Unitig-level bubble scoring

The compacted graph detector follows bounded, internally linear unitig alternatives that leave and rejoin the same oriented unitigs. It skips cycles, complex tangles, overlapping paths, and alternatives beyond the configured search bounds.

Each path is scored by nucleotide length, edge count, total and minimum support, mean edge support, local coverage relative to the strongest path, and sequence similarity to that path. Scores can be written to TSV, but no unitig or link is removed. This separation allows cleaning thresholds to be validated against sequencing errors, damage, and genuine rare-strain variation before they can affect N50 or accuracy.

The first five-seed validation found substantial coverage overlap between detected error paths and genuine 5× rare-strain paths. A maximum relative coverage of 0.20 selected 13 of 17 detected error paths but also 2 of 39 rare paths; a cutoff of 0.30 selected 16 errors but 11 rare paths. Both classes had approximately 0.976 sequence similarity for single-base alternatives. Coverage and sequence similarity alone therefore do not support safe automatic removal.

An expanded 20-seed matrix added bidirectional reads, three error rates, three damage intensities, and rare strains at 1×, 2×, 5×, and 10×. Relative terminal enrichment reduced biological-path selection from 84 of 451 with coverage alone to 20 of 451, while selecting 111 of 182 detected error paths. However, 15 of those 20 biological losses occurred at 1×–2× rare-strain coverage. No tested rule selected any error path while guaranteeing zero biological loss across the complete matrix.

The result supports graded classification rather than a universal hard threshold. Error-like paths can be simplified only under an explicitly chosen sensitivity/retention policy; ambiguous low-coverage alternatives should remain represented. Terminal damage produced no bubbles at any tested intensity and still requires separate tip and incomplete-branch logic.

Removing a false alternative would reconnect the surrounding unitigs and can increase contig length and N50 after recompaction. Removing a genuine strain path would produce the same topological improvement while destroying biological variation, so contiguity cannot be optimized independently of rare-strain retention.

## Explainable path classification

Detected unitig alternatives are now labelled `error-like`, `damage-like`, `variation-like`, or `ambiguous`; the locally strongest path is labelled `dominant`. Every decision records its evidence reasons and sequence substitutions. Classification is non-destructive and does not alter graph topology, contigs, or N50.

The conservative error-like rule combines low local coverage, high sequence similarity, and terminal depletion. A damage-like label additionally requires damage-compatible C→T/G→A substitutions, low local coverage, strand skew, and terminal enrichment. Supported, strand-balanced alternatives are variation-like; conflicting or insufficient evidence remains ambiguous.

Across the 220-condition matrix, the error-like label selected 111 of 182 detected error paths and 20 of 451 biological paths. Tightening the damage rule eliminated damage-like calls among the biological controls. However, terminal-damage simulations produced no detectable unitig bubbles, so the current results establish specificity only—not damage sensitivity. Damage detection must next operate on tips and incomplete branches before cleaning removes them.

## Simple bubble detection

The experimental detector reports bounded alternatives that leave one oriented handle, remain linear and edge-disjoint internally, and rejoin at the same handle. Reverse-complement observations are deduplicated through their shared physical edge paths.

Detection is non-destructive: it records edge paths and support summaries without deactivating edges or changing unitigs. Cycles, long alternatives, and complex tangles are skipped rather than guessed. The detector is a prerequisite for later evidence-aware decisions, not a bubble-removal rule.

## Terminal-evidence reporting

Weak tips and simple bubble paths are collected before graph cleaning. Their sequences, edge support, terminal/internal evidence, and substitutions are streamed to TSV. Equal-length bubble alternatives are marked damage-compatible only when C→T or G→A substitutions have support at the corresponding oriented molecule end.

Reporting does not change the graph. Tips currently remain unclassified when no competing reference path exists, and no evidence-based simplification decision is implemented yet.

## Evidence from the current baseline

Minimum-support filtering showed that unsupported k-mers cause most initial graph fragmentation. On the clean single-genome benchmark, `min_count=2` increased genome recovery from 21.444% to 94.940% without introducing a reported misassembly.

A fixed-k sweep showed that `k=31` gives the best current contiguity for 150 bp reads at 20× coverage. Larger k values recover slightly more of the reference but fragment valid paths. This result is dataset-specific and does not establish a default for ancient metagenomes.

The directed filtered assemblies have a duplication ratio near 2.0. The orientation-aware graph corrected this to 1.005 and reduced physical nodes from 9.92 million to 5.01 million. Its compact representation preserved the same unitigs while reducing peak memory from 6.92 GiB to 2.78 GiB and wall time from 193.90 seconds to 113.99 seconds.

Across two clean bacterial genomes, two seeds, and 5×/20× coverage, conservative tip cleaning introduced no reported misassemblies. At 20×, N50 improved by 56–67% and the largest contig by 52–128%. At 5×, cleaning made few changes and did not reduce genome fraction. This supports the current rule on clean isolates but does not establish safety for damaged molecules or rare strains.

Controlled bubble tests showed zero detections in clean conditions. Sequencing errors and low-abundance strain SNPs both produced simple bubbles, confirming that topology and abundance alone cannot distinguish them. Recovery of 12 introduced SNPs rose from a mean of 1.0 bubble at 1× minor-strain coverage to 11.6 at 10×.

Terminal damage produced hundreds of weak tips but no simple bubbles under the tested `k=21`, `min_count=2` conditions. Damage-aware decisions must therefore evaluate tips and incomplete branches as well as bubbles, and must operate before irreversible tip removal. The current count threshold also limits recovery of very low-coverage strain alternatives.

In a ten-seed controlled validation, terminal-only candidate edges recovered 75.24% of retained damage-derived edges with 98.52% pooled precision against clean and sequencing-error controls. All 5,460 detected damage candidates were terminal-only; 82 terminal-only candidates occurred in error controls. An empirical five-replicate comparison of matched non-UDG and full-UDG R1 subsets showed significantly more normalized bubbles and branching in non-UDG data, but terminal-event enrichment was not significant. These results validate the evidence primitive, not a final damage classifier.

## Near-term development sequence

1. Add and validate incomplete-branch and tip classification for damage topologies not represented by bubbles.
2. Add molecule/link and quality evidence to improve low-coverage error-versus-variation separation.
3. Implement optional conservative simplification only for high-confidence error-like paths.
4. Validate N50, accuracy, and strain retention on simulated and empirical mixtures.
5. Evaluate controlled multi-k assembly and profiled Cython optimization where needed.

Every algorithmic change will be compared with the frozen baseline for genome recovery, contiguity, misassemblies, mismatch rate, low-abundance retention, runtime, and peak memory.
