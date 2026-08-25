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

### Changed

- Replaced repeated graph-wide degree scans with linear-time degree calculation.
- Avoided duplicate unitig extraction during CLI assembly.
