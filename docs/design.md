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

When event reporting is enabled, source-read identity and side-specific Phred quality are retained for terminal observations in compact edge-indexed arrays. The quality corresponds to the nucleotide appended in the relevant oriented traversal, rather than an arbitrary endpoint of the supporting k-mer. Paired-fragment links are not yet stored. A candidate-locus damage profile can now be inferred from matched alternatives, but it is not yet a whole-library calibrated damage model.

The compact implementation represents canonical k-mers with rolling 2-bit integers. Physical edge endpoints and strand-support counts are stored in typed arrays, while node degrees and traversal state use byte arrays. This replaces per-node counters, string k-mers, tuple edge keys, and per-edge support objects without changing graph decisions.

## Compacted unitig graph

After k-mer-level cleaning and event reporting, each active maximal non-branching path is represented as one physical unitig. Oriented unitig handles retain links in both directions without duplicating reverse-complement sequences. Each unitig stores its sequence, path length, total, minimum, and mean edge support, plus aggregated terminal and internal observations.

The CLI now writes sequences from this representation. No simplification decision is made at the unitig level yet, so compaction must preserve the previous FASTA exactly. The connected unitig graph is the intended substrate for local-coverage bubble and branch decisions because those operations can inspect thousands of unitigs rather than millions of k-mer edges.

## Conservative tip cleaning

The experimental support-only cleaner considers only dead-end paths no longer than `2 × k`. A path is removed when its strongest edge has at most 20% of the support of a competing edge at the attachment branch. Long tips, isolated paths, and strongly supported alternatives are preserved.

The damage-aware mode applies the same topology and support boundary, then
matches each candidate to its local backbone and classifies the evidence before
deactivation. Only an `error-like` tip with no support from the opposite read
orientation is removable. Damage-like, variation-like, ambiguous, unmatched,
and bidirectionally supported candidates are cached as protected. Both modes
run for at most five rounds so newly exposed tips can be considered without an
unbounded cleanup pass.

Edges are deactivated directly in the compact adjacency arrays, and degrees and active graph totals are updated in place. This avoids copying the graph or allocating per-edge cleaning objects. Both policies remain opt-in: the support-only mode is a useful contiguity comparison, while `--damage-aware-clean-tips` is the conservative prototype policy.

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

Reporting itself does not change the graph. Tips without a suitable linear
competitor remain unmatched and are protected by the separate damage-aware
cleaner. When an event report is requested, retained terminal edge observations are now stored in compact edge-indexed arrays with their source-read index, oriented end, and end distance. For paths with multiple substitutions, the classifier intersects source-read identities and adds linkage evidence only when the same reads support every expected C→T/G→A change. A single substitution receives no linkage score bonus because its exact edge already associates it with a read.

FASTQ observations additionally retain the Phred quality of the exact nucleotide represented by each oriented terminal-edge observation. Event reports include observed and missing quality counts, mean quality, and the sum of `log(10^(-Q/10) / 3)` over damage-compatible alternative observations. This is the likelihood of the specific alternative calls under an independent, symmetric sequencing-error model. FASTA input remains supported and produces explicit missing-quality evidence instead of an assumed quality. The value is report-only: it does not include damage or polymorphism likelihoods, does not assume independence between overlapping graph events, and does not alter the heuristic classification.

On the paired `SRR32866683` public dataset, this path retained 7,561 qualities across 2,504 candidate rows with no missing values and weighted mean Q38.01. All topology and event counts matched the previous run, and the independently computed candidate damage fit remained unchanged. Most event rows correctly have no quality likelihood because only an exact, damage-compatible changed-base observation is currently eligible; ordinary substitutions, unmatched paths, and bubble-only reporting are not silently assigned evidence.

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

## Held-out event likelihoods

Damage-profile schema version 6 adds deterministic five-fold cross-fitting. Unique weak-tip loci are sorted by their alternative edge, reference edge, and damage channel, then assigned round-robin to folds. Each locus is scored with the geometric damage curve fitted on the other four folds. Incomplete branches are absent from profile training and use the full tip-derived profile, which is independent of those branch events. If a training fold fails the existing evidence gates, affected events remain explicitly insufficient.

For each exact C→T/5′ or G→A/3′ locus, the scorer uses the same terminal alternative and backbone molecules under three explanations. Following the latent-base emission used by CarpeDeam, a damage probability `p` and Phred error `e` give `P(alt) = p(1-e) + (1-p)e/3`; the corresponding reference probability reverses those terms. The sequencing-error model fixes `p=0`. The variation model replaces the distance-dependent curve with one fitted constant alternative frequency bounded at 0.5 and subtracts `0.5 log(n)` from its log likelihood, the one-parameter BIC penalty. All models therefore see the same observations and qualities.

The output records raw log likelihoods, the penalized variation score, fitted variation frequency, best-ranked explanation, score margin, observation counts, model scope, and insufficiency reasons. It is not a posterior: priors, correlated errors, locus overdispersion in predictive scoring, internal-read quality evidence, and a calibrated abstention margin remain unresolved. Existing heuristic labels and graph topology are intentionally unchanged.

In the 80-run controlled matrix, all 20 damage-generated events ranked damage, all 20 independent Phred-error events ranked sequencing error, and all 20 constant-frequency events ranked variation. All 20 high-quality errors deliberately generated with the same positional curve as damage ranked damage. This negative control demonstrates identifiability rather than an implementation defect: position, substitution type, and base quality alone cannot tell a damage process from a systematic error engineered to have the same distribution.

Five-fold fitting of the public `SRR32866683` candidate loci succeeded in every fold. Each model trained on 3,613–3,614 of 4,517 loci and held out 903–904. Mean fitted probabilities over cycles 0–4 were 0.3367, 0.2691, 0.2167, 0.1761, and 0.1444; fold LR statistics ranged from 33.34 to 44.42. The close fold curves support stable held-out scoring on this sample, while remaining candidate-conditioned rather than whole-library calibrated.

The full schema-version-6 public rerun preserved the assembly byte-for-byte and preserved all 45 earlier event fields across 111,399 rows. It scored 1,478 events: 968 ranked damage, 24 sequencing error, and 486 variation; 95,848 remained explicitly insufficient. Of the scored events, 71.31% had a margin below 1, and none of the damage-ranked events reached a margin of 5. These results validate non-regression, provenance, and conservative abstention, but not a classification boundary. Reporting time increased from 66.24 to 131.18 seconds and the TSV grew by 50.30%, making batched scoring and a more compact optional report the main performance follow-ups.

## Calibrated event confidence

Calibration is a separate report-processing stage rather than part of graph construction. A schema-version-1 calibration model stores class-conditional nonconformity distributions learned from independently identified samples. For a proposed class, nonconformity is the strongest competing log score minus that class's score, divided by the square root of the molecule count. Damage uses the lower cross-fit likelihood when asking whether damage is supported and the upper likelihood when asking whether damage can be rejected. This envelope is deliberately conservative and exposes profile instability instead of treating one fitted curve as exact.

Each report event receives conformal p-values for damage, sequencing error, and variation. Damage and variation p-values are additionally corrected across the report with the Benjamini–Hochberg procedure. An event becomes `eligible_error` only when sequencing error remains conforming while both damage and variation are rejected, cross-fit damage instability is below its gate, and independent molecule minima are met. The defaults require five total observations, two alternative observations, and five reference observations. Plausible damage is `protect_damage`, plausible variation is `protect_variation`, and missing, unstable, or unresolved evidence is `insufficient`. These are still report-only decisions.

Experiment 14 calibrated on 5,000 independently seeded samples spanning fitted damage, a misspecified damage curve, independent Phred errors, damage-shaped systematic errors, and constant-frequency variation. On 100 new samples per condition, all fitted-damage and 99 misspecified-damage cases were protected as damage; the remaining misspecified case was protected as variation. All independent-error cases were eligible, all systematic damage-shaped errors were protected as damage, and all variation cases were protected, although six received the more conservative damage protection label. No protected condition became cleaning-eligible.

Applying the model offline to the public schema-version-6 report scored the same 1,478 events but made none cleaning-eligible after the molecule gates. It protected 340 as damage and 93 as variation; 1,045 scored events and all unscored rows remained insufficient. Before the alternative-molecule gate, eleven rows passed the likelihood tests, but every one was supported by exactly one alternative molecule. This confirms why total coverage cannot substitute for replicated alternative evidence.

The first graph-derived adversarial matrix then separated calibration and validation references and reused exact introduced-event k-mer truth. Three calibration seeds supplied 2,762 scored damage events, 24 sequencing-error events, and 12 damage-compatible rare-strain events, sufficient for an exploratory `alpha=0.1` conformal resolution. Across two unseen references, 1,909 damage, 26 error, and six variation truth events were scoreable. No event became `eligible_error`; 21 damage events were protected and the rest remained insufficient. Only two independent-error events reached likelihood scoring, both singleton alternatives. The remaining scored errors were deliberately systematic terminal errors and lacked the five reference molecules required for a decision. This result is safe but has zero useful error sensitivity: replicated graph molecule evidence and broader bubble scoring are now the limiting factors, not the calibration threshold.

Direct bubble-path comparison, ordinary-substitution scoring, and physical-fragment identities improve evidence recovery without changing graph topology. In the repeated three-seed/two-seed matrix, calibration retained 2,762 scored damage and 12 variation events while sequencing-error examples increased from 24 to 54. On unseen references, damage and variation remained unchanged at 1,909 and six scoreable events, while sequencing errors increased from 26 to 47. The evidence funnel recorded 4,812 retained damage-truth edges → 2,837 detected events → 2,805 matched events → 1,909 scored events; the corresponding error funnel was 2,065 → 153 → 150 → 47, and the rare-strain funnel was 798 → 39 → 39 → six. The added errors still had too little replicated support for cleaning and generally ranked as variation. This exposes event ascertainment as the next statistical problem: after conditioning on a graph branch being detectable, a raw per-base error probability understates the expected alternative fraction.

On `SRR32866683`, the expanded scorer preserved all 111,399 report rows and increased scoreable events from 1,478 to 25,415: 5,692 tips, 19,716 incomplete branches, and seven bubble alternatives. It directly matched 3,552 of 3,706 non-reference bubble paths, although only seven contained both alternative and reference evidence inside the five-base terminal window. The final rankings were 24,453 variation, 852 damage, and 110 sequencing error. These are diagnostic rankings, not cleaning decisions; 14,537 scored rows had only one or two alternative fragments.

The v5, v6, and optimized v7 public contigs and damage profiles were byte-identical. Expanding likelihood coverage initially increased event reporting from 131.18 to 244.46 seconds. Exact aggregation of repeated beta-binomial locus-count patterns and repeated allele-quality observations reduced it to 135.13 seconds without changing scores. Total runtime fell from 567.21 to 446.85 seconds, leaving the expanded report only 3.0% slower than the earlier 1,478-event report.

The event-conditioned scorer evaluates each model given that at least one consistent alternative fragment and one consistent reference fragment were observed. This corrects the raw error model's failure to account for graph-event ascertainment. A high-quality singleton remains non-identifiable as error versus rare variation, so the report now abstains instead of forcing that case into an error rank. Phred quality at or below 20 can still supply direct singleton-error evidence; this is an explicit conservative reporting boundary, not a cleaning threshold.

On the held-out graph matrix, conditioning improved damage-truth rankings from 1,644 to 1,680 of 1,909 and recovered 21 sequencing-error rankings. Conservative singleton abstention retained 17 Q20 independent-error rankings, changed four Q23 independent errors to ambiguous, and changed all four affected rare-strain events from erroneous error rankings to ambiguous. The six scoreable variation events then contained four ambiguous and two damage rankings, with none ranked error. No event became cleaning-eligible, and deliberately damage-shaped systematic errors remained unresolved.

The `SRR32866683` v9 run preserved byte-identical contigs and damage profile, 111,399 rows, and 25,415 scoreable events. It reported 646 damage, 470 error, 20,428 variation, and 3,871 ambiguous rankings. All ambiguous rankings were high-quality singletons. Event reporting took 158.15 seconds versus 254.33 for the first conditioned implementation and 135.13 for v7; total runtime was 464.59 seconds versus 585.55 and 446.85, respectively.

The graph-context follow-up extracts existing evidence without changing the graph or production classifier: event topology, path length, substitution count, alternative and reference path observations, minimum edge support, relative coverage, branch degree, sequence identity, terminal fraction and enrichment, strand balance, damage distance, joint-fragment support, and mean base quality. Experiment 15 writes these fields per truth-labelled validation event plus median summaries and pairwise probability-of-superiority measurements.

The expanded matrix holds rare-strain mixtures at 20× total depth and includes ordinary and damage-compatible variants at 1×, 2×, 5×, and 10×, three independent error rates, independent terminal errors, systematic paired terminal errors, and three damage strengths. Scored calibration examples increased from 54 to 152 for error and from 12 to 26 for variation; validation examples increased from 47 to 103 and from six to 17, respectively, while damage remained 2,762/1,909. No variation event ranked sequencing error: 13 were ambiguous, three ranked damage, and one ranked variation.

Variation was strongly separated by path length, alternative/reference path observations, terminal fraction, terminal enrichment, and relative coverage. Median variation paths contained 21 edges and terminal fraction 0.33, compared with one edge/1.0 for damage and two edges/1.0 for error. Damage versus error remained less separable, and independent terminal errors continued to resemble damage or variation. Graph context is therefore suitable first as variation-protection evidence, not as a standalone deletion rule. The funnel also exposes the next bottleneck: 96 variation events were detected and matched but only 17 were scoreable inside the five-base terminal likelihood window.

Targeted whole-read recovery closes much of that evidence gap for ordinary
substitutions. Event matches are computed once, their changed allele edges are
collected, and the loaded reads are rescanned once. Evidence is aggregated to
one best-quality call per physical molecule and edge, rather than storing a
Python object for every k-mer occurrence. Damage-compatible C→T/G→A scoring
continues to use only orientation-correct terminal cycles.

In the repeated three-calibration/two-validation matrix, scoreable error events
increased from 152/103 to 301/202 and variation from 26/17 to 78/53 for
calibration/validation. Correct validation rankings increased from 19 to 49
for error and from one to 44 for variation. Damage remained unchanged at
2,762/1,909 scored events and 1,680 damage versus 229 variation rankings.

The `SRR32866683` v12 run preserved byte-identical v10 contigs and damage
profile and all 111,399 event rows. Scoreable rows increased from 25,415 to
78,663; error rankings increased from 470 to 2,994 and variation rankings from
20,428 to 71,021, while damage remained effectively stable at 646 versus 645.
Because the public library has no event truth, these are evidence-coverage and
non-regression results rather than biological-accuracy estimates. The targeted
read pass increased event reporting from 171.65 to 547.07 seconds and total
runtime from 504.61 to 853.12 seconds, with 7,011,012 KiB peak RSS and no swap.

Matched-tip heuristic classification can now change topology only when
`--damage-aware-clean-tips` is explicitly enabled. The event-likelihood and
conformal-calibration layers remain report-only and are not used as deletion
probabilities. Simulation truth validates selection behavior; reference-backed
assemblies validate sequence accuracy; an unlabeled public library cannot
establish biological accuracy. That still requires independently characterized
untreated and UDG-treated libraries plus expanded held-out community mixtures.

The complete 20-seed tip-selection matrix contained 36,977 candidates across
360 runs. The initial error-like-only rule selected 97 simulated errors but also
one rare-strain candidate. That rare candidate had balanced support from both
read orientations. Requiring one-sided support retained 82 simulated errors,
protected the 15 bidirectional error-like candidates, and selected zero of
35,202 damage candidates and zero of 111 rare-strain candidates.

Four clean reference-backed assemblies then exercised actual graph mutation at
5× and 20× coverage. The damage-aware policy removed no tips at 5× and 10 and
seven tips in the two 20× replicates. It introduced no misassemblies, mismatches,
indels, or genome-fraction loss. N50 was unchanged at 5× and increased only from
4,648 to 4,665 bp and from 4,680 to 4,703 bp at 20×. The policy is therefore
safe but intentionally low-recall in this validation; it is not a replacement
for the larger contiguity gains of support-only cleaning.

On `SRR32866683`, v13 removed 989 tips containing 1,121 edges in 7.93 seconds.
Unitigs fell from 740,785 to 739,660 and assembled bases from 24,239,229 to
24,215,608, while N50 remained 32 bp and the largest unitig remained 461 bp.
The pre-cleaning event report and damage profile were byte-identical to v12, as
expected from their pipeline position. Total wall time increased from 13:08.86
to 13:28.28 and peak RSS from 7,011,012 to 7,033,856 KiB. This public run proves
execution and deterministic pre-cleaning evidence, not removal accuracy.

## Incomplete-branch matching

The event detector now also follows bounded weak paths leaving a branch until their next graph boundary. It excludes paths already represented as short weak tips, alternatives that rejoin as simple bubbles, cycles, paths exceeding 64 edges, and candidates without an equal-length linear backbone comparison. Because detection is non-destructive, it admits paths supported at up to 50% of the matched backbone rather than reusing the cleaner's 20% deletion threshold.

Across 20 clean, standard-damage, 1% sequencing-error, and 2× rare-strain replicates, existing tips and bubbles covered 10,584 of 14,194 retained damage-derived edges (74.57%). Incomplete branches added 2,394 distinct truth edges and raised combined coverage to 12,978 edges (91.43%). All 1,482 damage-condition branch candidates overlapped damage truth; 557 were damage-like and 925 remained ambiguous. No candidate appeared in clean controls. None of 72 error branches or two genuine rare-strain branches was damage-like.

These results validate broader event recovery, not removal. The matrix uses one damage profile, one error rate, and one rare-strain coverage. Incomplete branches now receive molecule-linked evidence during event reporting, but their linkage-specific sensitivity and safety still need the expanded adversarial matrix used for weak tips.

## Evidence from the current baseline

Minimum-support filtering showed that unsupported k-mers cause most initial graph fragmentation. On the clean single-genome benchmark, `min_count=2` increased genome recovery from 21.444% to 94.940% without introducing a reported misassembly.

A fixed-k sweep showed that `k=31` gives the best current contiguity for 150 bp reads at 20× coverage. Larger k values recover slightly more of the reference but fragment valid paths. This result is dataset-specific and does not establish a default for ancient metagenomes.

The directed filtered assemblies have a duplication ratio near 2.0. The orientation-aware graph corrected this to 1.005 and reduced physical nodes from 9.92 million to 5.01 million. Its compact representation preserved the same unitigs while reducing peak memory from 6.92 GiB to 2.78 GiB and wall time from 193.90 seconds to 113.99 seconds.

Across two clean bacterial genomes, two seeds, and 5×/20× coverage, support-only conservative tip cleaning introduced no reported misassemblies. At 20×, N50 improved by 56–67% and the largest contig by 52–128%. At 5×, cleaning made few changes and did not reduce genome fraction. The new damage-aware policy was much more selective: zero removals at 5× and 17 total at 20×, with 0.4–0.5% N50 gains and no accuracy regression. The 20-seed truth-labelled selection audit additionally observed no damage or rare-strain selections after the bidirectional guard.

Controlled bubble tests showed zero detections in clean conditions. Sequencing errors and low-abundance strain SNPs both produced simple bubbles, confirming that topology and abundance alone cannot distinguish them. Recovery of 12 introduced SNPs rose from a mean of 1.0 bubble at 1× minor-strain coverage to 11.6 at 10×.

Terminal damage produced hundreds of weak tips but no simple bubbles under the tested `k=21`, `min_count=2` conditions. Damage-aware decisions must therefore evaluate tips and incomplete branches as well as bubbles, and must operate before irreversible tip removal. The current count threshold also limits recovery of very low-coverage strain alternatives.

In a ten-seed controlled validation, terminal-only candidate edges recovered 75.24% of retained damage-derived edges with 98.52% pooled precision against clean and sequencing-error controls. All 5,460 detected damage candidates were terminal-only; 82 terminal-only candidates occurred in error controls. An empirical five-replicate comparison of matched non-UDG and full-UDG R1 subsets showed significantly more normalized bubbles and branching in non-UDG data, but terminal-event enrichment was not significant. These results validate the evidence primitive, not a final damage classifier.

## Near-term development sequence

1. Extend the validated one-sided error-like protection gate to bounded
   incomplete branches, first as a truth-labelled selection audit and then as
   an opt-in topology change.
2. Validate incomplete-branch removal against damage, independent errors,
   rare strains, and clean reference-backed assemblies before public-data use.
3. Recalibrate ordinary whole-read events separately from terminal
   damage-compatible events and expand held-out references, coverage regimes,
   library protocols, and recurrent non-systematic error fixtures.
4. Validate the damage model on independently characterized untreated,
   partial-UDG, and full-UDG libraries, then add whole-library opportunities
   and bootstrap or profiled predictive uncertainty.
5. Extend incomplete-branch matching to unequal-length and locally non-linear
   backbones only after the bounded equal-length policy is validated.
6. Evaluate calibrated `eligible_error` decisions as an optional stricter
   simplification policy once error sensitivity is adequate.
7. Validate N50, accuracy, and strain retention on simulated and empirical
   mixtures, including direct comparison with MEGAHIT, metaSPAdes, and
   CarpeDeam.
8. Profile a native or compiled whole-read evidence scanner; retain the current
   pure-Python targeted pass as the auditable reference implementation.

Every algorithmic change will be compared with the frozen baseline for genome recovery, contiguity, misassemblies, mismatch rate, low-abundance retention, runtime, and peak memory.

## Reciprocal paired-read unitig extension

The optional paired-extension stage maps each read's nearest retained 3′ k-mer
to an oriented compacted unitig. Each physical pair contributes one
reverse-symmetric unitig observation. At a branching junction, only mate targets
reachable through exactly one outgoing candidate are discriminating; shared
repeat targets abstain. A candidate must win with at least five independent
pairs and 3:1 support over the runner-up, and the reverse orientation must select
the reciprocal link before sequences are joined.

Reachability is target-aware and capped at 1,000 visited states per junction.
Exceeding the cap marks the junction ambiguous. This reduced the public
`SRR32866683` extension stage from more than 50 minutes without completion to
34.22 seconds at `k=21`, with 4,003 of 21,014 evaluated orientations safely
work-limited. That run accepted 17 physical joins and left N50 at 32 bp.

The final accuracy gate used the existing `GCF_001050915.2`, 20×, seed-42 ART
reads at `k=31`. Damage-aware cleaning alone was compared with the identical
pipeline plus paired extension using QUAST 5.3.0:

| Metric | Damage-aware baseline | Strict paired extension |
|---|---:|---:|
| Largest contig | 24,470 bp | 24,470 bp |
| N50 | 4,665 bp | 4,733 bp |
| Genome fraction | 92.952% | 93.014% |
| Duplication ratio | 1.005 | 1.005 |
| Misassemblies | 0 | 0 |
| Mismatches per 100 kbp | 0.00 | 0.00 |
| Indels per 100 kbp | 0.00 | 0.00 |

The strict run accepted 387 physical joins. A 3-pair/2:1 rule increased N50 to
5,296 bp but introduced 0.43 mismatches and 0.04 indels per 100 kbp. A simplified
unique-predecessor path-context experiment increased N50 to 5,297 bp but caused
one misassembly, 1.50 mismatches, and 0.04 indels per 100 kbp. Both were rejected.
Correct exSPAnder-style growing-path context requires insert-size-aware mapping
and is not approximated by propagating evidence through a topologically unique
chain.

## Direct read-thread unitig extension

The direct threading stage indexes only physical graph edges at compacted-graph
junction boundaries. Source reads are streamed again, and consecutive indexed
k-mers provide oriented unitig transitions. Evidence is deduplicated first per
read and then per physical molecule. A branching source is resolvable only when
one target has at least five molecules and 3:1 support over the runner-up; the
reverse orientation must independently choose the reciprocal physical link.
The audit retains zero-support candidates, orientation counts, terminal and
internal support, base quality, local mean edge support, rank, and dominance.

Reciprocal links are converted into reverse-symmetric successors and spelled
with the same overlap-aware path routine as paired extension. Before spelling,
the resolver examines every internal unitig in a proposed joined chain. If its
mean edge support is at least twice that of both flanks, it is treated as a
local multi-copy repeat signal and the weaker adjacent physical link is
discarded. This is a local ratio rather than a global coverage threshold, which
is important for uneven-coverage metagenomes. The ratio is configurable through
`--thread-repeat-coverage-ratio`; the validated default is 2.0.

A same-read adjacent-junction gate was tested first and rejected. It removed 19
links, reduced reference N50 from 18,021 to 17,474 bp, and retained the same
7 bp indel. The discrepancy lay inside a 2,530 bp repeated unitig, while the
short internal hub had mean support 39.43 versus 15.22 and 19.63 on its flanks.
The final coverage guard removed one unsafe link, restored N50 to 18,021 bp,
and reduced QUAST indels from one to zero without changing the 81,176 bp
largest contig or introducing mismatches or misassemblies.

Experiment 17 adds truncated log-normal short fragments, bidirectional terminal
damage, sequencing errors, a related strain, and an unrelated contaminant.
Across 20 runs, 6,989 reciprocal links had 6,987 exact-reference matches. The
two exact failures were C/T- or G/A-compatible paths in one damaged mixture;
all 6,989 links were correct under R/Y topology truth, and extension did not
increase false-contig counts. Exact and R/Y results are reported separately so
validation does not optimize away the damage signal the assembler must retain.

On public `SRR32866683` at `k=31`, the guard rejected six of 32,669 reciprocal
links. Relative to fixed-k damage-aware output, joined contigs fell from 639,369
to 606,706, N50 increased from 38 to 40 bp, and the largest contig increased
from 366 to 435 bp. Relative to unguarded threading, N50 and largest-contig
length were unchanged. Because this sample has no assembly truth, the public
run validates execution, conservative selection, and resource behavior only.
