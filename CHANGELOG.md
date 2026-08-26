# Changelog

All notable project changes will be documented in this file.

## Unreleased

### Added

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

### Changed

- Replaced repeated graph-wide degree scans with linear-time degree calculation.
- Avoided duplicate unitig extraction during CLI assembly.
- Replaced string-based orientation-aware k-mers with rolling 2-bit encoding.
- Replaced nested orientation-aware graph objects with compact integer arrays and byte-array traversal state.
- Reduced the orientation-aware benchmark wall time by 41% and peak memory by 60% without changing QUAST accuracy metrics.
