# Changelog

All notable project changes will be documented in this file.

## Unreleased

### Added

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

- Replaced repeated graph-wide degree scans with linear-time degree calculation.
- Avoided duplicate unitig extraction during CLI assembly.
- Replaced string-based orientation-aware k-mers with rolling 2-bit encoding.
- Replaced nested orientation-aware graph objects with compact integer arrays and byte-array traversal state.
- Reduced the orientation-aware benchmark wall time by 41% and peak memory by 60% without changing QUAST accuracy metrics.
- Replaced full node decoding during path spelling with constant-time oriented base extraction.
- Routed orientation-aware CLI output through the compacted unitig graph without changing FASTA output.
