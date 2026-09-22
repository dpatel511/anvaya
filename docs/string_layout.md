# Exact overlap path control

Sources: [Myers 2005, The fragment assembly string graph](https://doi.org/10.1093/bioinformatics/bti1114)
describes constructing string graphs and transitive reduction from read overlaps.
[Li 2016, Minimap and miniasm](https://pmc.ncbi.nlm.nih.gov/articles/PMC4937194/)
separates overlap layout from base correction and reduces graph redundancy.
Their results do not establish performance for Anvaya's short ancient reads.

The synthetic-only `--string-layout` mode reuses Anvaya's bidirected edge helper
and `_project_master_edges`. It builds exact read dovetails, collapses identical
sequences up to reverse complement while retaining every raw observation, removes
two-edge transitive overlaps with matching shifts, and projects reciprocal
nonbranching paths. It does not pop bubbles or resolve repeats by guessing.
This is a bounded control, not a complete Myers algorithm or miniasm port.
Equal-length reads and at most 1,000 reads are required. Candidate discovery uses
the existing k15/eight-anchor finder; it is not exhaustive overlap enumeration.

The existing projector now optionally returns each output's node orientations
and offsets. The experiment transforms these into original raw placements
(0-based, half-open raw intervals), preserving reverse orientation and duplicate
read evidence for one consensus pass. No simulation origins enter layout.
Singletons are retained, unlike the progressive baseline's cluster admission,
so total contig counts and coverage require that qualification. Progressive
diagnostics in the case summary describe the separately computed baseline;
the `graph` field and output scores describe the graph experiment.

Verified matrix: `results/controlled-string-v1`, 36 cases. Clean unique 20x
controls become one contig of 1,198/1,200 bp with zero base errors, matching their
truth-connected spans. Clean strain controls retain 84/85 contigs (N50 99/99)
with no conflicting-origin contigs. Repeat controls have four contigs and one
conflicting-origin contig each; that remains unresolved evidence, not proof of
a false join. The damaged unique controls have N50 80/80, 163/164 contigs and
128/124 resolved mismatches, unchanged by consensus. Exact graph placement
fragments the evidence needed for correction.

This is the first layout control that reconstructs the high-coverage unique
reference while abstaining at strain branches. It is not a damage-aware
assembler result. Next, validate damage-compatible edges while preserving
competing branches and placement provenance; do not relax all mismatches or
add a polishing loop. New seeds and realistic fragment lengths are required
before promoting an edge policy. The v1 hashes precede a storage-allocation-only
change in the experiment; the verified path/branch unit tests cover that change.

## Phase 1 ambiguity and provenance hardening

Verified in `results/controlled-string-phase1-final` on 2026-09-08. All 36
tracked assembly/evaluation cases are metric-identical to `controlled-string-v1`.
Clean unique 20x remains one 1,198/1,200-base contig with zero resolved errors;
damaged unique 20x remains fragmented at N50 80 and is intentionally unresolved
by this phase. The final matrix records source hashes, including:

- `experiments/string_layout.py`: `550125fa535a4e43e74933a0fa6e0e7eaaedadad3b4ee664e1355c18f08e1d8e`
- `src/anvaya/overlap_graph.py`: `fcc00377785d032cf316b6ef752a86d4931fa5f8109139a2076cbda5c03433e6`
- `src/anvaya/overlap_progressive_links.py`: `c1e708ec5241a9f6effc0e4bbb31e76eccd6e280d3cbc83788659c2d03f10e97`

The bidirected-edge helper now has an opt-in ambiguity set. When two different
edge records compete for one physical oriented pair, it removes both directions
and blacklists that physical pair instead of allowing the later offset to
overwrite the first. Existing production callers omit the set and retain their
previous behavior; the bounded string-layout experiment enables it and reports
`ambiguous_physical_edges`.

When separate physical paths spell the same sequence, path projection now unions
their physical nodes and coordinate-normalized layouts instead of retaining only
the first path's provenance. Reverse-complement paths are transformed into the
stored output orientation using `contig_length - offset - read_length`. Identical
layout tuples are deduplicated. The corrected-base diagnostic takes the maximum
across identical spelled paths rather than counting the same output correction
twice.

The experiment accepts optional molecule IDs and records contributing molecules
from those IDs, not raw read indices. Palindromic raw sequences retain both
possible orientations, causing the existing consensus layer to treat the
molecule placement as ambiguous instead of assigning an arbitrary physical end.

Focused tests cover clean/transitive chains, branches, cycles, conflicting pair
offsets, identical spelled paths with distinct physical provenance, reverse
duplicates, duplicate molecule IDs, and palindromic orientation. The full suite
ran 239 tests with no failures or skips using pysam 0.24.0 in a temporary Python
3.13 environment. No actual FASTQ data was run. Phase 1 meets its acceptance
gate. Phase 2 adds and tests allele-compatible transitive reduction below.

## Phase 2 damage-compatible edges

The synthetic-only `--damage-string-layout` mode keeps exact candidates and
adds an approximate edge only when every overlap mismatch has one directional
damage explanation under the supplied profile. Evidence is taken from original
FASTQ observations aligned across the site, with original qualities and raw-end
coordinates. A latent base must have both an exact observation and a terminal
C>T or G>A observation with nonzero profile probability, and at least three
distinct unambiguous molecules must cover the site. Internal transitions,
ordinary substitutions, zero-rate damage events, insufficient evidence and
missing or low-quality-only evidence do not create an approximate edge.

Coordinates are 0-based. For an oriented raw alignment, original raw position
is `p` forward and `L - 1 - p` reverse. Damage is tested against raw 5' and 3'
distances, never distance to a graph-node or contig end. Each molecule supplies
at most one observation at a site; conflicting representations are excluded.

Approximate edges carry their inferred overlap alleles. A direct transitive
edge is removed only when a two-edge path has the same cumulative shift and no
conflicting inferred allele. Reciprocal correction tuples are coordinate-sorted
so discovering the same physical edge from opposite orientations cannot create
a false ambiguity. Path placements are retained for the existing single
consensus pass.

Development results are in `results/controlled-damage-string-phase2-v1`, using
the already-examined seeds 93001 and 93002. Against
`results/controlled-string-phase1-final`, all 18 clean cases are unchanged. Of
18 damaged cases, N50 improves in 14 and is unchanged in four; longest contig
improves in 16. Every damaged case has fewer absolute resolved mismatches and a
lower resolved-base mismatch rate. Aggregate resolved mismatches fall from
1,698 to 498. Repeat-panel conflicting-origin contigs do not increase.

The frozen rule was then evaluated without tuning on predeclared seeds 94011
and 94012. Results are in `results/controlled-damage-string-phase2-heldout` and
the exact comparison is in `results/controlled-string-phase2-heldout-control`.
All clean cases are unchanged. Damaged N50 improves in 15 of 18 cases and is
unchanged in three; longest contig improves in all 18. Aggregate resolved
mismatches fall from 1,883 to 554, no case worsens its resolved-base mismatch
rate, and repeat conflict totals fall from 17 to 13.

These are 80-base, error-free mechanism controls with a matched supplied damage
profile. Phase 3 variable-length results and their additional limitations are
recorded below.

Fifteen focused graph tests cover the damage directions, both orientations,
asymmetric and zero profiles, internal and ordinary mismatches, missing and low
quality, duplicate molecules, reciprocal edge representation, branches, cycles,
provenance and incompatible transitive alleles. The repository suite ran 247
tests successfully with pysam 0.24.0; one external wrapper smoke test was skipped
because minimap2 and samtools were unavailable in that temporary environment.
`git diff --check` passed apart from existing line-ending conversion warnings.

## Phase 3 variable-length prototype status

The implementation now lives in `src/anvaya/damage_string_graph.py`; the old
experiment file is a compatibility import. Dovetail shifts and overlaps use the
two participating read lengths. Raw placements retain each original `[0, L)`
interval. Unique exact containments are retired into their single longest
parent while preserving every molecule's offset and orientation. Equal-length
duplicates remain grouped. Repeated containments with multiple equally long
parents remain separate. Pool records are built in one batch rather than copied
once per read.

The exact 80-base oracle was rerun in
`results/controlled-string-phase3-source-regression`. All tracked metrics in all
36 cases match `controlled-string-phase1-final`. Twenty-one focused graph tests
cover variable-length chains, variable-length damage, containment coordinates,
ambiguous containment and the earlier Phase 1/2 invariants.
The repository suite passes 256 tests with one external remapping integration
skipped because `minimap2` and `samtools` are unavailable in the test environment.

Profiling a 560-read damaged case attributed 0.424 of 0.846 seconds to the two
candidate-alignment passes and 0.211 seconds to damage-edge classification. A
1,000-read 31–200 bp synthetic case took 1.03 seconds without profiler overhead.
A 10,000-read case took 14.25 seconds and reached 141 MiB maximum RSS. These are
scaling observations, not a 100k resource guarantee.

The package now has a dedicated `damage-assemble` command. It requires FASTQ
qualities and a supplied CarpeDeam-format profile, then emits pre-consensus and
single-pass consensus FASTA from one layout plus JSON, placement and decision
audits. Its default input cap remains 1,000 because the variable-length safety
gate below did not pass. The reader stops at the first record beyond the cap.

Development seeds 95031 and 95032 cover 31–200 bp fragments, matched damage,
Q10 substitution errors, a half-strength supplied profile, repeats and an 80:20
two-strain mixture. With the validated Phase 2 edge policy, damage-aware layout
reduces consensus mismatches in every damaged case and does not add repeat or
strain origin conflicts. Contiguity remains inconsistent: N50 rises in six of
ten damaged comparisons, is unchanged once and falls three times. Seed 95032's
strain case also loses 18 uniquely recovered reference bases.

An exact-backbone plus reciprocal unique-best bridge experiment was evaluated
and reverted. Although it improved the unique panels, it raised equal-length
repeat origin conflicts from 15 to 24 across the established 36-case oracle. It
also produced a 229-base cross-strain path through a locally identical 48-base
overlap whose flanking strain markers were not spanned by any read. Local edge
thresholds cannot establish that phase, so the source module retains the Phase 2
edge policy.

Two further branch-selection mechanisms were rejected rather than promoted.
Reciprocal longest-overlap filtering improved N50 and mismatches on every
held-out unique case, but seed 95103 gained a repeat-origin conflict and seed
95104 gained a strain-origin conflict while losing 84 uniquely recovered bases.
Requiring uniquely junction-spanning molecules directly between raw-read nodes
made supported damage edges largely redundant. Applying the repository's
raw-supported linker after exact unitig construction still increased strain
conflicts from 2 to 3 on seed 95031 and from 5 to 9 on seed 95104. Junction
support is therefore not haplotype support when competing strains share the
whole local overlap.

Phase 3 therefore remains incomplete. Do not run the 100k FASTQ yet. The next
implementation must carry branch or haplotype context across locally identical
overlap bridges. The next unused seed set, 95201–95204, is declared in
`experiments/variable_length_validation.py` and must remain unexamined until a
multi-site rule is frozen. Raising the CLI cap and issuing the real-data command
happen only after that gate.
