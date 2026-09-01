# Changelog

All notable project changes will be documented in this file.

## Unreleased

### Changed

- Established reciprocal ranked read extension with contig merging disabled as
  the safe overlap prototype. On EMN001 500k it improved N50 from 104 to 133,
  NA50 from 102 to 132, and aligned length from 3.83 to 4.90 Mb with zero
  MetaQUAST misassemblies. Contig merging remains explicit and experimental
  after adding one 240 bp translocation for only +3 N50.
- Corrected experiment 18's `summary.tsv` output to use tab separators.
- Replaced stale DBG roadmap statements with the iterative whole-fragment
  overlap-backend decision and separated the active experiment index from the
  preserved chronological research archive.
- Removed an unused private encoded-k-mer wrapper.
- Moved terminal-damage correction after overlap layout so polishing operates
  on frozen contigs and cannot change topology or length.

### Added

- Ranked overlap extension, supported extension consensus, candidate-rejection
  diagnostics, configurable candidate discovery, and an opt-in reciprocal-best
  gate that checks all fragments for incompatible near-best opposite paths.
- A separate experimental `overlap-assemble` command implementing conservative
  whole-fragment overlap clustering, iterative extension, reusable-read layout,
  and contig-overlap merging with DNA and R/Y identity gates.
- Opt-in frozen-layout T-to-C/A-to-G terminal-damage polishing, accepted-edit
  TSV reporting, a one-edit reference truth auditor, and exact reconstruction
  of paired unpolished contigs. The EMN001 500k audit resolved 1,596 edits with
  99.69% precision while preserving the complete contig topology.

- Experimental opt-in damage-aware pre-graph read consensus using a bounded
  canonical-anchor sketch, damage-tolerant R/Y overlap identity, independent
  high-quality nonterminal molecule support, conservative dominance, and a
  streamed per-base audit TSV. The frozen ladder records strong synthetic
  continuity gains but an 8.98% rare-strain recovery loss; `SRR32866683`
  completes within 6.38 GiB peak RSS but retains its 38 bp N50, so the feature
  remains research-only and disabled by default.
- Truth-only `oracle_unique` and `oracle_all_truth` validation-ladder arms for
  measuring the upper bound of filtered-edge rescue without changing the
  production assembler. The 30-run result rejects singleton dead-end rescue as
  a continuity strategy for damage and sequencing-error scenarios.
- P1 report-only dead-end attribution that rescans raw reads only for retained
  dead-end contexts and separates boundaries, coverage gaps, uniquely filtered
  continuations, filtered conflicts, and contexts retained elsewhere.
- P2 report-only solid-k-mer correction auditing for low-quality bases, with
  unique complete-rescue decisions, ambiguous and incomplete outcomes, and
  protection of terminal damage-compatible reversals.
- P1/P2 counts in the frozen validation-ladder output and CLI summaries, plus
  bounded evenly spaced read sampling for public-data correction audits.
- A versioned P0 validation ladder with five deterministic seeds across clean,
  terminal-damage, sequencing-error, combined-damage/error, rare-strain, and
  contamination controls, including truth-aware continuity and recovery
  metrics plus configuration and runtime manifests.
- A deterministic, report-only fragmentation TSV that attributes both physical
  ends of every unitig as dead ends, unique continuations, or ambiguous
  branches and summarizes isolated, one-sided, and connected topology.
- A conservative, report-only damage-projection gate that emits an auditable
  candidate decision only when damage classification, likelihood margin,
  variation contrast, and base-quality requirements all pass.
- Reproducibility manifests for baseline experiments, including runtime and Git
  state, resolved tool paths, commands, parameters, and artifact checksums.
- Opt-in direct read-thread extension across compacted-graph junctions, using
  independent molecule support, 3:1 strong-winner selection, reciprocal
  orientation agreement, and a boundary-only edge index.
- A local two-copy coverage guard that breaks the weaker adjacent join when an
  internal threaded unitig has at least twice the coverage of both flanks.
- Thread-audit TSV fields for directional molecule support, terminal and
  internal evidence, base quality, local unitig coverage, candidate rank,
  dominance, resolvability, and reciprocal selection, plus clean and
  ancient-like truth-labelled validation experiments 16 and 17.
- Opt-in reciprocal paired-read extension of compacted unitig paths, with
  orientation-aware molecule anchors, strong-winner scoring, bounded target
  searches, reciprocal joins, and explicit search-work audit counters.
- Conservative paired-extension defaults of five independent pairs and 3:1
  dominance, selected by reference-backed QUAST validation; weaker and
  path-context variants are documented as rejected accuracy experiments.
- Optional iterative damage-aware tip cleaning through
  `--damage-aware-clean-tips`, removing only short, weak, one-sided tips
  classified as error-like while protecting damage-like, variation-like,
  ambiguous, unmatched, and bidirectionally supported alternatives.
- Cleaning summaries with round, removal, and protection counts, plus
  controlled damage/error/rare-strain selection tests and reference-backed
  QUAST validation of the new policy.
- Targeted whole-read exact-base evidence for ordinary substitutions, collected
  in a second read pass only for matched allele edges and deduplicated by
  physical molecule, while damage likelihoods remain restricted to oriented
  terminal cycles.
- Reusable precomputed tip, incomplete-branch, and bubble-path matches during
  event reporting, avoiding repeated path matching during evidence collection
  and TSV generation.
- Report-only graph-context extraction for matched events, including topology,
  path multiplicity, local coverage, fragment support, terminal behavior,
  strand balance, damage distance, and base quality.
- Expanded experiment 15 with fixed-depth ordinary and damage-compatible
  rare-strain mixtures, independent and terminal error processes, per-class
  feature summaries, and pairwise feature-separation reports.
- Event-ascertainment-conditioned damage, sequencing-error, and variation
  likelihoods evaluated over physical fragments, with the conditioning scope
  and per-model selection probability written to event reports.
- Conservative singleton abstention: a high-quality alternative supported by
  one fragment is reported as ambiguous because sequencing error cannot be
  distinguished from genuine rare variation, while low-quality singleton
  evidence may still rank as error.
- Conditioning provenance in experiment 15 validation tables and shared
  fragment grouping across likelihood models to reduce repeated scoring work.
- Direct likelihood comparison of equal-length bubble alternatives against the
  dominant path, including production event-report fields for matched and
  scoreable bubble paths.
- Ordinary-substitution error-versus-variation likelihood scoring and an
  evidence-funnel report for retained, detected, matched, and scored graph
  truth events.
- Separate read indices and physical-fragment identifiers, with paired mates
  sharing one identifier and contradictory fragment-level allele evidence
  excluded from support counts.
- Report-only class-conditional conformal calibration for event likelihoods,
  with held-out sample provenance, BH-corrected damage and variation tests,
  cross-fit damage instability, and explicit `protect_damage`,
  `protect_variation`, `eligible_error`, or `insufficient` decisions.
- A separate `anvaya calibrate-events` command and experiment 14 covering
  independent errors, systematic damage-shaped errors, variation, and damage
  profile misspecification.
- A graph-derived adversarial calibration experiment using exact damage/error
  k-mer truth, damage-compatible rare-strain variants, and matched tip,
  incomplete-branch, and bubble-path observations.
- Independent minimum alternative and reference molecule gates that prevent
  high total coverage from making singleton alternatives cleaning-eligible.
- Deterministic five-fold candidate-locus damage fitting that excludes each
  weak-tip locus from the profile used to score it.
- Report-only damage, Phred sequencing-error, and constant-frequency variation
  likelihood rankings with an explicit variation-model complexity penalty,
  evidence scope, insufficiency reasons, and uncalibrated margins.
- Damage-profile schema version 6 cross-fit provenance plus controlled and
  public-profile validation experiments 12 and 13.
- Compact, side-specific Phred quality storage for molecule-linked terminal
  observations, with explicit missing-quality accounting and a report-only
  symmetric sequencing-error log likelihood in event TSV output.
- A 20-condition quality-evidence validation matrix covering Phred 8, 20, and
  35 support plus FASTA-style missing qualities across one to five molecules.
- Candidate-conditioned beta-binomial geometric damage fitting, a
  constant-background null, likelihood-ratio evidence, bounded conditional
  support intervals, and explicit sufficiency gates in profile schema version 5.
- Clean, ordinary terminal-error, and correlated terminal-error likelihood
  controls in experiment 10.
- Exact changed-base cycle validation, separate endpoint-enrichment and
  monotonicity metrics, and matched-versus-observed locus counts in
  damage-profile schema version 4.
- Per-locus damage evidence and parallel equal-locus and coverage-capped
  summaries in damage-profile schema version 3.
- Per-bin edge-contribution diagnostics in damage-profile schema version 2,
  including contributor counts and the largest single-edge contribution.
- Initial Stage 0 research, design, benchmark, and literature documentation.
- Validated overlapping k-mer extraction.
- Basic directed de Bruijn graph construction and edge multiplicity.
- Node degree calculations and source, sink, and branch identification.
- Maximal non-branching unitig extraction, including cycles and self-loops.
- Graph topology and support summaries.
- End-to-end reads-to-unitigs assembly workflow.
- Plain and gzipped FASTA/FASTQ input with Phred quality scores.
- DNA normalization, reverse complement, and canonical-sequence utilities.
- Ambiguous k-mer exclusion.
- Deterministic test-only helpers for fragments, errors, terminal damage, and SNPs.
- Initial graph-disturbance smoke experiment.
- Unit tests for assembler primitives and simulation helpers.
- `anvaya assemble` command with progress reporting and FASTA output.
- Optional minimum k-mer support filtering through `--min-count`.
- Experimental orientation-aware graph construction through `--orientation-aware`.
- Canonical physical graph nodes with oriented traversal handles.
- Forward, reverse, and palindromic support tracking for canonical k-mers.
- Orientation-aware unitig extraction and graph summaries.
- Experimental conservative tip cleaning through `--clean-tips`.
- In-place compact-graph edge deactivation and active graph accounting.
- Reproducible clean tip-validation matrix with verified stage resumption and TSV summaries.
- Non-destructive bounded simple-bubble detection for the compact bidirected graph.
- Optional bubble reporting through `--detect-bubbles`.
- Bubble tests covering support, bounds, orientation deduplication, determinism, CLI reporting, and unchanged assembly output.
- Optional compact read-end evidence collection through `--end-window`.
- Non-destructive tip and bubble-path TSV reporting through `--event-report`.
- Sequence, support, substitution, and damage-compatibility fields for graph alternatives.
- Controlled ten-seed terminal-damage validation with reproducible parameter, per-run, and summary tables.
- Unit tests for end-evidence orientation, unchanged assembly output, event sequences, and CLI validation.
- Compacted oriented unitig graph with links, coverage summaries, and terminal/internal evidence.
- Unitig-graph tests covering branches, cycles, reverse links, evidence, and sequence equivalence.
- Non-destructive bounded bubble detection and scoring on the compacted unitig graph.
- Optional unitig-path support, local-coverage, and sequence-similarity TSV reporting through `--unitig-bubble-report`.
- Unitig-bubble tests covering bounds, scoring, unchanged graph state, CLI validation, and byte-identical assembly output.
- Public oriented unitig-path spelling for ground-truth sequence validation.
- Reproducible five-seed unitig-bubble matrix covering clean, sequencing-error, terminal-damage, and rare-strain conditions.
- Candidate-level and threshold-grid TSV reports for measuring error selection against rare-strain loss.
- Orientation-correct forward/reverse and left/right terminal evidence on compacted unitigs.
- Terminal fraction, relative internal coverage, terminal enrichment, and strand-balance bubble-path scores.
- Bidirectional controlled-read simulation and an extended 220-condition validation matrix.
- Non-destructive unitig-path classification as dominant, error-like, damage-like, variation-like, or ambiguous.
- Auditable classification reasons, substitutions, and damage compatibility in unitig-bubble TSV reports.
- Classification summaries in the CLI and controlled validation matrix.
- Non-destructive weak-tip matching to equal-length local backbone paths.
- Tip-event report fields for matched backbone sequence and support, relative coverage, DNA identity, RY identity, substitutions, and oriented damage compatibility.
- Controlled ten-seed tip-matching validation covering clean, terminal-damage, and random-error conditions.
- Unit tests for tip matching, damage compatibility, unchanged graph state, event-report fields, and CLI summaries.
- A 20-seed weak-tip matching matrix with separate 5′, 3′, and combined damage profiles, bidirectional reads, ordinary and terminally concentrated errors, and rare strains at four coverages.
- Candidate-level tip-matching reports and regression tests for 3′ G→A and multiple oriented damage substitutions.
- Non-destructive matched-tip evidence scoring using exact oriented substitution edges, end distance, terminal/internal support, coverage, identity, and strand balance.
- Conservative matched-tip classifications and auditable reasons in event reports and validation tables.
- A 20-seed classification matrix measuring damage sensitivity, error specificity, and rare-strain safety without graph modification.
- Non-destructive bounded incomplete-branch detection with standard-tip and simple-bubble deduplication.
- Equal-length incomplete-branch backbone matching, evidence classification, event-report fields, and CLI summaries.
- A 20-seed incomplete-branch matrix covering clean, damage, sequencing-error, and rare-strain controls.
- Compact source-read linkage for retained terminal edge observations and multi-substitution damage evidence.
- Optional candidate-locus damage-profile JSON reporting through `--damage-profile-report`.
- Strand-symmetric terminal damage fractions with separate auditable C→T/5′ and G→A/3′ channels.
- A 20-seed profile-shape validation across low, standard, and high simulated damage.

### Changed

- Compressed repeated beta-binomial locus records and repeated allele-quality
  observations during likelihood evaluation. A controlled profile fell from
  110.13 to 25.50 seconds with byte-identical matrix outputs; public event
  reporting fell from 244.46 to 135.13 seconds while retaining 25,415 scored
  events.
- Ordinary-only substitutions can no longer rank as damage because a
  zero-damage model ties the sequencing-error model.
- Replaced repeated graph-wide degree scans with linear-time degree calculation.
- Avoided duplicate unitig extraction during CLI assembly.
- Replaced string-based orientation-aware k-mers with rolling 2-bit encoding.
- Replaced nested orientation-aware graph objects with compact integer arrays and byte-array traversal state.
- Reduced the orientation-aware benchmark wall time by 41% and peak memory by 60% without changing QUAST accuracy metrics.
- Replaced full node decoding during path spelling with constant-time oriented base extraction.
- Routed orientation-aware CLI output through the compacted unitig graph without changing FASTA output.
