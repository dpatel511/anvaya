# Evidence and implementation priorities

Review date: 2026-09-06. Scope: reference-free assembly of short ancient bacterial
and archaeal fragments, initially untreated double-stranded libraries with merged
reads. References are permitted for evaluation, not for directing assembly.
The objectives are unique aligned recovery, accurate sequence, preservation of
low-abundance strain alleles, and practical resource use. N50 alone is not success.

## What the current experiments establish

The user ran the EMN001 100k experiments. These are observations from one subset,
not independent replications or held-out validation:

| Configuration | Aligned bases | Mismatches / 100 kbp | Genome fraction (%) |
|---|---:|---:|---:|
| Neither new quality filter | 1,320,332 | 988.31 | 0.478 |
| Whole-seed rejection only | 1,199,082 | 1,020.61 | 0.443 |
| Later-recruit filtering only | 1,319,481 | 987.21 | 0.478 |
| Recruit filtering plus dynamic trimming projection | 1,321,565 | 985.95 | 0.478 |

Whole-seed rejection loses substantial aligned recovery. Trimming removed 276
bases from 259 candidate contigs, but the dynamic projection admitted 27 more
contigs, so its apparent gain cannot be attributed to base removal alone.
The fixed-membership comparison selects untrimmed contigs once and applies
trimming to those same records. Length exclusions apply to both paired outputs
and are reported. Separate total-recovery comparisons remain necessary.
Even fixed membership does not fix reference alignment placement or denominators.

## Literature and repository survey

This is a targeted survey of primary papers, repository documentation and selected
source routines, not an exhaustive review or a benchmark of downloaded tools.
No third-party code was vendored or new runtime dependency installed.

### CarpeDeam: principal damage-aware comparator

The paper separates raw-fragment extension, where end-specific damage information
is usable, from contig merging. It uses a supplied damage profile in correction
and extension, and describes a safe-mode consensus step. This supports preserving
original fragment positions rather than assigning damage probabilities to derived
contig ends. It does not establish that our fixed penalties implement the same
model. [Paper, Methods](https://doi.org/10.1186/s13059-025-03839-5).

Source inspected at commit `ab0a2c083a9d6b11566ff2154589af011edecefc`:
[ancientReadsResults.cpp](https://github.com/LouisPwr/CarpeDeam/blob/ab0a2c083a9d6b11566ff2154589af011edecefc/src/assembler/ancientReadsResults.cpp),
[correction.cpp](https://github.com/LouisPwr/CarpeDeam/blob/ab0a2c083a9d6b11566ff2154589af011edecefc/src/assembler/correction.cpp), and
[nuclassembleUtil.cpp](https://github.com/LouisPwr/CarpeDeam/blob/ab0a2c083a9d6b11566ff2154589af011edecefc/src/assembler/nuclassembleUtil.cpp).
The examined extension path calls `calcLikelihoodConsensus`; the utility composes
damage transitions with a sequencing-error matrix. `getSeqErrorProf` receives a
scalar error rate. Do not equate this routine with per-observation FASTQ-quality
weighting. Reuse the statistical decomposition and comparator infrastructure;
do not translate its entire assembler into Python. Its public repository reports
GPL-3.0 licensing. [Repository](https://github.com/LouisPwr/CarpeDeam).

### PenguiN / PLASS: candidate discovery and extension ranking

PenguiN combines overlap assembly and Bayesian extension selection, with protein-
guided and nucleotide-only workflows. Protein guidance is not a default fit for
our ultrashort fragments; nucleotide overlap discovery is the relevant comparison.
[Paper](https://pubmed.ncbi.nlm.nih.gov/39354646/),
[repository workflows](https://github.com/soedinglab/plass).

At commit `ac83d8fa2fd315738afe675609a3dea05837cfea`, the examined
[nuclassembleresult.cpp](https://github.com/soedinglab/plass/blob/ac83d8fa2fd315738afe675609a3dea05837cfea/src/assembler/nuclassembleresult.cpp)
compares overlap mismatch/length evidence through beta-distribution calculations
and uses a priority queue for extensions. Anvaya's mean-minus-standard-deviation
ranking is not that comparison. Evaluate an existing clustering/search backend
against an exhaustive small-fixture overlap oracle before replacing our index.
The repository reports GPL licensing; no implementation was copied.

### StrainXpress: preserve alternative paths

The method clusters reads, assembles locally, then joins through a master overlap
graph while stopping at ambiguous paths. Its published preprocessing retains reads
longer than 70 bp and its global joining thresholds are stringent. Those settings
would exclude much of our target data. Reuse the separation of local assembly and
strain-safe global joining as a design principle, not those numerical thresholds.
[Paper](https://doi.org/10.1093/nar/gkac543),
[repository](https://github.com/HaploKit/StrainXpress).

### mapDamage and PyDamage: inference and authentication

mapDamage provides a statistical damage model and damage-aware quality rescaling.
Its reference-aligned observations are not interchangeable with our selected de
novo overlaps. Start with an externally supplied profile whose library type and
provenance are explicit, before attempting joint assembly/profile estimation.
[Paper](https://doi.org/10.1093/bioinformatics/btt193),
[repository](https://github.com/ginolhac/mapDamage).

PyDamage authenticates contigs from mapped-read damage patterns. Its accuracy
depends on coverage, contig length and damage level; it does not test assembly
correctness. Use it as an optional downstream check when data support inference,
not as a hard gate on our mostly very short contigs. If many contigs are tested,
report adjusted q-values and reliability alongside results rather than treating
every nominal significance result as reliable.
[Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8323603/),
[repository](https://github.com/maxibor/pydamage),
[output documentation](https://pydamage.readthedocs.io/en/0.80/output.html).

### Existing alignment and simulation components

- [Edlib](https://github.com/Martinsos/edlib) provides fast edit-distance paths and
  prefix/infix/global alignment. It is a candidate geometry verifier, not a
  position- and quality-dependent damage scorer. Do not make C/T universally
  equivalent to gain sensitivity: that also erases strain evidence.
- [Parasail](https://github.com/jeffdaily/parasail) supplies SIMD pairwise alignment.
  Evaluate it if affine-gap alignment is needed; a conventional substitution
  matrix alone does not encode fragment-end position and individual base quality.
- [Minimap2](https://github.com/lh3/minimap2) is useful established alignment
  infrastructure, but successful short-read mapping is not proof of sufficient
  all-vs-all overlap sensitivity for 31–100 bp damaged fragments. Measure recall
  on that range before adopting it as our candidate generator.
- [gargammel](https://github.com/grenaud/gargammel) separates fragmentation, damage,
  adapter addition and sequencing errors. Prefer it for biologically realistic
  external challenges. Keep tiny in-house fixtures for exact unit-test truth,
  preserving source IDs, orientation and fragment coordinates in evaluation data.

## Implementation sequence and acceptance checks

1. **Finish the fixed-membership trimming audit.** Run matched before/after
   evaluation; reconcile all candidate exclusions. Report mismatches and aligned
   bases together, and separately report sequences below the evaluator's useful
   alignment-length range. Stop tuning trimming if effects remain negligible.
2. **Track raw evidence per position.** Extend the current seed-only mapping into
   read placements: read/molecule ID, orientation, original read interval, contig
   offset and alignment path if gapped. Keep raw qualities immutable. Use 0-based,
   half-open internal intervals; derive original read-end distances through the
   mapping. Read ends equal physical ends only when preprocessing establishes it.
   Test reverse strands, merged reads, duplicates and successive extensions.
3. **Implement and test one observation likelihood.** Initially accept a validated
   external damage profile. A proposed decomposition is
   `P(observed | true, position, quality) = sum_b D(true,b,position) E(b,observed,quality)`.
   Damage and sequencing error are separate processes. Test normalization, zero
   damage, direction-specific terminal transitions, high-quality strain variants,
   missing qualities and numerical stability. This is a proposed model, not a
   claim that nominal probabilities are calibrated on selected overlaps.
4. **Apply it to consensus, then extension, separately.** Count independent
   molecules rather than duplicate evidence; preserve unresolved alleles. Do not
   assume the growing center is always correct. Compare overlap scores against
   competing origins and length/background effects; calibrate acceptance on held-
   out truth rather than calling a heuristic margin a correctness probability.
5. **Measure candidate recall and optimize the measured bottleneck.** Use exhaustive
   tiny fixtures to quantify losses from anchors, occurrence caps and greedy
   ownership. Benchmark existing search/alignment libraries before writing a new
   aligner. Profile pool copying and alignment work; batch updates or replace a
   backend only when profiling supports the change, keeping outputs controlled.

Validation must cross damage, ordinary error, fragment length, coverage, repeats,
strain distance and abundance. Separate tuning and evaluation seeds, references
and datasets. Measure unique reference coverage, allele retention, false joins,
paired base errors, runtime and RSS. Genome fraction must use unions within each
reference ID; aggregate aligned length may include duplicated sequence. Report
uncertainty across independent simulation replicates, not p-values based on
treating correlated bases as independent replicates. Do not fit a damage model
and claim independent validation on the same accepted overlap evidence.

The [first raw-placement consensus projection](raw_consensus.md) now implements
a supplied-profile observation likelihood and conservative calling rules on the
fixed progressive layout. It is opt-in and awaits actual-data evaluation. Complete
strain-mixture inference, duplicate identification, gapped layouts and calibration
remain open. A new trim threshold or a longer contig alone is not that milestone.
