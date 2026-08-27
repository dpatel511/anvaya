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

Base quality and paired-fragment links are not yet stored. Source-read identity is retained for terminal observations when event reporting is enabled, and a candidate-locus damage profile can now be inferred from matched alternatives. This profile is not yet a whole-library calibrated damage model.

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

Weak tips are also matched non-destructively to an equal-length locally competing linear backbone when one exists. Candidate selection prioritizes RY identity, then DNA identity and backbone support. The report records both sequences, substitutions, relative coverage, oriented damage compatibility, terminal/internal support, end distance, strand balance, three evidence scores, and an auditable classification.

In the original 20-seed matrix, 30,960 of 34,762 matched damage tips were classified damage-like (89.06%). Thirty-two of 1,641 matched error tips and none of 110 genuine rare-strain tips were called damage-like. Error-like calls additionally require terminal depletion: this protects low-abundance variation but deliberately leaves most errors ambiguous. These are heuristic evidence scores, not calibrated probabilities or an authorization to remove graph paths.

Reporting does not change the graph. Tips without a suitable linear competitor remain unmatched, and no evidence-based simplification decision is implemented yet. When an event report is requested, retained terminal edge observations are now stored in compact edge-indexed arrays with their source-read index, oriented end, and end distance. For paths with multiple substitutions, the classifier intersects source-read identities and adds linkage evidence only when the same reads support every expected C→T/G→A change. A single substitution receives no linkage score bonus because its exact edge already associates it with a read.

## Matched-tip evidence classification

Matched weak tips receive normalized heuristic scores for damage, sequencing error, and biological variation. The damage score combines oriented C→T/G→A compatibility, RY identity, support for the substitution on its exact terminal edge, distance from the expected molecule end, and low local coverage. Error and variation scores reuse coverage, terminal/internal support, terminal enrichment, sequence identity, and strand balance. Threshold, evidence, and margin gates keep conflicting cases ambiguous.

In a 20-seed, 460-condition paired ablation, molecule linkage increased overall damage-like recall from 30,960/34,762 (89.06%) to 31,018/34,762 (89.23%). For the 1,155 multi-substitution damage tips, recall increased from 1,093 (94.63%) to 1,151 (99.65%). All 58 label changes were true damage candidates moving from ambiguous to damage-like, distributed across 18 of 20 seeds; no candidate moved in the opposite direction. The paired-error and rare-strain labels did not change. Build time was 2.04× higher while linkage was enabled, with about 672 KB mean compact-array storage per 20 kb/20× graph.

This validates read-level co-occurrence, not biochemical causality. The expanded controls contained 166/1,810 error tips and 1/124 rare tips called damage-like under both scoring modes. In particular, 134/169 deliberately correlated terminal-error tips were already damage-like. Those errors and the molecule-spanning rare SNP pair were represented as separate single-substitution tips, so they did not exercise the multi-substitution linkage gate. A controlled graph fixture or broader path matcher that produces linked and unlinked multi-substitution non-damage alternatives is still required before linkage can justify graph simplification.

## Candidate-locus damage profile

The optional `--damage-profile-report` output aggregates molecule-linked evidence across unique matched weak-tip substitution loci. Incomplete branches are excluded until profile-shape recovery is validated for that detector. Alternative and reference edge evidence are deduplicated independently so a shared backbone edge cannot be weighted once per competing alternative. For each terminal distance the report records C→T/5′ and G→A/3′ alternative observations, matched backbone observations, and their fractions. Because a bidirected path can be represented in either reverse-complement orientation, the primary `damage_fraction` combines both equivalent channels; separate raw channels remain in the JSON for auditing.

Schema version 4 maps edge-boundary observations to the exact appended nucleotide cycle before accepting them. This matters for 5′ C→T evidence because a k-mer may overlap the terminal window while its changed base lies outside it. The report distinguishes `matched_loci` from `observed_loci`; the backward-compatible `eligible_loci` field now equals the latter. It retains the alternative and reference distance vectors for every unique edge-pair/substitution locus and reports observation-weighted, equal-locus, and 20-observation-capped summaries. Contributor and largest-edge fields expose concentration directly. These summaries are diagnostics for selecting and validating a later overdispersed decay model; the capped and equal-locus curves are not used as biological estimates or graph-cleaning inputs. Bins without observations are reported as `null` rather than zero.

Across 20 random references for each low, standard, and high bidirectional damage profile, mean Pearson correlation with the known simulated profile was 0.908, 0.966, and 0.974 respectively. Every run had a higher terminal than interior endpoint, but exact monotonic decay occurred in 12/20 low-damage runs and all standard/high runs. Low-damage profiles averaged 0.4 upward steps and 46.8 observed loci; standard and high profiles had no upward steps and averaged 351.4 and 686.0 loci. These results validate recovery of profile shape from matched candidates when signal is sufficient, not absolute whole-library damage frequency. Unmatched genomic opportunities, base quality, locus dependence, and sampling uncertainty are not yet included.

Schema version 5 adds a report-only candidate-conditioned likelihood model. Each nonempty locus-distance count is modeled with a beta-binomial distribution whose mean follows `background + amplitude * decay**distance`; a constant-background beta-binomial model is fitted as the null. The report includes both log likelihoods, their likelihood ratio, bounded conditional likelihood-support intervals, and explicit insufficiency reasons. Fits require at least 10 observed loci, three distances, three alternative observations, and 20 total observations. These are prototype evidence gates rather than calibrated biological thresholds.

All 60 low/standard/high damage simulations fitted successfully. Mean correlations between the fitted and generating curve shapes were 0.9994, 0.9988, and 0.9994, while mean likelihood-ratio statistics were 26.40, 178.73, and 431.39. Clean controls had no observed loci, and ordinary terminal-error controls averaged 2.3 loci, so all 40 remained insufficient rather than producing fitted evidence. Deliberately correlated terminal C→T/G→A errors did fit in all 20 runs and had a mean statistic of 5.62 (range 0.00–12.97), narrowly overlapping the low-damage range of 12.78–44.52. A universal decision threshold is therefore not justified, and the fit is not yet used by the classifier.

The classifier does not yet change topology or improve N50 directly. It creates an auditable decision layer that can later protect damage-like and ambiguous paths while a separately validated simplifier targets only high-confidence errors.

## Incomplete-branch matching

The event detector now also follows bounded weak paths leaving a branch until their next graph boundary. It excludes paths already represented as short weak tips, alternatives that rejoin as simple bubbles, cycles, paths exceeding 64 edges, and candidates without an equal-length linear backbone comparison. Because detection is non-destructive, it admits paths supported at up to 50% of the matched backbone rather than reusing the cleaner's 20% deletion threshold.

Across 20 clean, standard-damage, 1% sequencing-error, and 2× rare-strain replicates, existing tips and bubbles covered 10,584 of 14,194 retained damage-derived edges (74.57%). Incomplete branches added 2,394 distinct truth edges and raised combined coverage to 12,978 edges (91.43%). All 1,482 damage-condition branch candidates overlapped damage truth; 557 were damage-like and 925 remained ambiguous. No candidate appeared in clean controls. None of 72 error branches or two genuine rare-strain branches was damage-like.

These results validate broader event recovery, not removal. The matrix uses one damage profile, one error rate, and one rare-strain coverage. Incomplete branches now receive molecule-linked evidence during event reporting, but their linkage-specific sensitivity and safety still need the expanded adversarial matrix used for weak tips.

## Evidence from the current baseline

Minimum-support filtering showed that unsupported k-mers cause most initial graph fragmentation. On the clean single-genome benchmark, `min_count=2` increased genome recovery from 21.444% to 94.940% without introducing a reported misassembly.

A fixed-k sweep showed that `k=31` gives the best current contiguity for 150 bp reads at 20× coverage. Larger k values recover slightly more of the reference but fragment valid paths. This result is dataset-specific and does not establish a default for ancient metagenomes.

The directed filtered assemblies have a duplication ratio near 2.0. The orientation-aware graph corrected this to 1.005 and reduced physical nodes from 9.92 million to 5.01 million. Its compact representation preserved the same unitigs while reducing peak memory from 6.92 GiB to 2.78 GiB and wall time from 193.90 seconds to 113.99 seconds.

Across two clean bacterial genomes, two seeds, and 5×/20× coverage, conservative tip cleaning introduced no reported misassemblies. At 20×, N50 improved by 56–67% and the largest contig by 52–128%. At 5×, cleaning made few changes and did not reduce genome fraction. This supports the current rule on clean isolates but does not establish safety for damaged molecules or rare strains.

Controlled bubble tests showed zero detections in clean conditions. Sequencing errors and low-abundance strain SNPs both produced simple bubbles, confirming that topology and abundance alone cannot distinguish them. Recovery of 12 introduced SNPs rose from a mean of 1.0 bubble at 1× minor-strain coverage to 11.6 at 10×.

Terminal damage produced hundreds of weak tips but no simple bubbles under the tested `k=21`, `min_count=2` conditions. Damage-aware decisions must therefore evaluate tips and incomplete branches as well as bubbles, and must operate before irreversible tip removal. The current count threshold also limits recovery of very low-coverage strain alternatives.

In a ten-seed controlled validation, terminal-only candidate edges recovered 75.24% of retained damage-derived edges with 98.52% pooled precision against clean and sequencing-error controls. All 5,460 detected damage candidates were terminal-only; 82 terminal-only candidates occurred in error controls. An empirical five-replicate comparison of matched non-UDG and full-UDG R1 subsets showed significantly more normalized bubbles and branching in non-UDG data, but terminal-event enrichment was not significant. These results validate the evidence primitive, not a final damage classifier.

## Near-term development sequence

1. Expand the candidate-conditioned likelihood model with whole-library opportunities and fully profiled or bootstrap uncertainty.
2. Replace heuristic scores with calibrated damage, sequencing-error, and variation likelihoods using held-out profile fit.
3. Extend incomplete-branch matching to unequal-length and locally non-linear backbones.
4. Implement optional conservative simplification only for high-confidence error-like paths.
5. Validate N50, accuracy, and strain retention on simulated and empirical mixtures, including direct comparison with MEGAHIT, metaSPAdes, and CarpeDeam.
6. Evaluate paired-read linkage, controlled multi-k assembly, and profiled native-code optimization where needed.

Every algorithmic change will be compared with the frozen baseline for genome recovery, contiguity, misassemblies, mismatch rate, low-abundance retention, runtime, and peak memory.
