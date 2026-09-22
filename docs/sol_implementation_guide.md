# Sol implementation guide — active contract, 2026-09-20

This is the ONLY active plan. Read `AGENTS.md`, this document and the
[2026-09-20 review](astra_review_20260920.md). The
[previous guide](archive/sol_implementation_guide_pre_20260920.md) is historical,
including all its appended next actions. Its relative links were written from
`docs/`. Implementer observations belong in [findings](sol_implementation_findings.md).

## Goal and reviewer decision

Deliver an independent overlap assembler with correct variable-length geometry,
immutable raw molecule provenance, damage-aware decisions and honest output
accounting. CarpeDeam remains a comparator, never the production engine. No
particular improvement or release date is guaranteed.

This review authorizes ONE architectural candidate: a conventional exact unitig
core separated from raw evidence, followed by one constrained damage integration.
It does NOT reopen failed containment-suppression patches, post-layout linking,
greedy branch projections or paired opportunity sweeps. Do not lower thresholds,
restart diagnostic variants or select favorable seeds.

The reviewer reproduced 138 correct clean mate merges with downstream N50
2,399→99; the merged graph retains 64 directed full-read overlap edges and 265
ambiguous containments. This disproves the merger-failure interpretation but does
not prove a simple edge filter solves the graph. See
`experiments/review_20260920.py` and `results/reviewer-20260920/`. That consumed
fixture is a regression, NOT validation.

## Execution rules

1. Preserve unrelated dirty work; standard branch prefixes only, never `codex/`.
   Do not stage/commit unless requested. Keep ignored local notes out of Git.
2. Inspect referenced functions before edits. Proposed files below are NEW files.
   Reuse existing records, orientation helpers, raw projection and runners.
3. One implementation, minimal reproductions, no new audit framework. Fix ordinary
   bugs against declared invariants; failed biological gates close the candidate.
4. Agent runs synthetic/unit tests. User runs actual FASTQ benchmarks. Existing
   results may be inspected; do not launch actual datasets yourself.
5. Coordinates: 0-based half-open, orientation explicit. Both mates represent ONE
   molecule. Damage distances refer to original raw ends, not assembled ends.

## Phase A — verify exact geometry first

Deliver a small exact unitig core outside the default production route. Proposed
module: `src/anvaya/exact_string_graph.py`; proposed tests:
`tests/test_exact_string_graph.py`. No public backend flag yet.

Read `damage_string_graph.py::assemble`,
`overlap_graph.py::_project_master_edges`, and candidate/index helpers. Consult
[Readjoiner](https://link.springer.com/article/10.1186/1471-2105-13-82) and
[SGA's official README](https://github.com/jts/sga/blob/master/src/README) for
standard overlap/unitig semantics. Pin and inspect any specific implementation
used as reference; check licensing before copying. Independently implement the
needed algorithm, not a wrapper around another assembler.

Required invariants:

- Separate identical groups, containment, proper dovetails and unrelated matches.
  For oriented A starting at 0 and B starting at s, a proper right dovetail has
  `0 < s < len(A)` and `s + len(B) > len(A)`. Overlap `len(A)-s` is shorter than
  BOTH reads. Full-read overlap is containment, never an extension edge. Emit
  geometrically consistent reverse-complement mirror edges.
- Separate topology from immutable raw evidence. Retain exact containment
  placements needed to identify maximal strings and compose offsets/orientations.
  Multiple containers are not automatically distinct loci. Never allocate an
  ambiguous raw placement to an arbitrary parent.
- Contained strings need not be branch vertices in this CLEAN exact core. This
  does not authorize dropping their alleles or claiming damage/strain safety.
  Evidence allocation remains a mandatory, separately tested Phase B contract.
- Reduce only sequence/geometry-consistent transitive edges. Traverse maximal
  unbranched unitigs symmetrically, including isolated nodes and cycles. Do not
  greedily choose longest paths through branches. Preserve layout relationships
  when identical output sequences arise; string equality is not common origin.

Verification:

1. Add a simple exhaustive exact-overlap oracle in TESTS only, <=100 synthetic
   reads. Compare indexed candidates, proper-edge geometry and reduction/path
   invariants. This independent oracle justifies a small test-only implementation;
   do not use all-pairs scanning as the production candidate index.
2. Cover variable-length containment chains, two overlapping containers at one
   locus, repeated containers, prefix/suffix containment, duplicates, reverse
   orientation, branches, cycles and disconnected components. Check permutation
   and reverse-complement invariance of sequences and layouts.
3. Run the consumed clean merge fixture. Target recovery of the unmerged 2,399-bp
   sequence with correct layouts and no introduced sequence. This is a target,
   not a promised result. If it fails, use the oracle to identify the discrepancy;
   do not tune thresholds to that seed. No truth enters the assembler. The review
   script asserts OLD defects intentionally: retain it as historical evidence,
   and add candidate tests separately.
4. Report edge types and sequence/layout checks, not only N50. Measure synthetic
   scaling. Index once; avoid global scans per edge/transition.

Decision: correct oracle agreement and recovered clean target → Phase B. Oracle
mismatch → fix the implementation. Correct oracle also fails target → document
one concrete obstruction and STOP for review; no mate/support threshold search.
No real-file run in Phase A.

## Phase B — integrate damage evidence once

Deliver one candidate using verified geometry and current raw likelihood and
consensus, keeping the safe baseline reproducible. Read `raw_consensus.py`, raw
placement types, current damage edge scoring and
`tests/test_containment_survival_fixture.py` first.

[CarpeDeam](https://link.springer.com/article/10.1186/s13059-025-03839-5) motivates
terminal-damage-aware extension, not guaranteed strain safety.
[PEAR](https://pmc.ncbi.nlm.nih.gov/articles/PMC3933873/) and
[AdapterRemoval](https://adapterremoval.readthedocs.io/en/2.3.x/manpage.html)
inform mate overlap/quality treatment; retain BOTH original mate coordinate maps.

- Reuse existing likelihood and molecule support gates. Mismatches do not turn
  containment into a dovetail. Exact sequence agreement is not proof of a correct
  allele. Do not count merged consensus as two independent observations.
- Preserve contained reads' correction dependencies at exposed primary sites.
  Compose raw placements and offsets without arbitrary allocation or duplicated
  votes. Unresolved placements remain unresolved. A sidecar cannot compensate
  for missing correct bases in the scored primary output.
- Keep the baseline four-read regression unchanged. Add candidate biological
  invariant tests under permutation/reverse complement; different representations
  need not preserve old edge/path counts. Reproduce consumed strain-interval and
  reference-position 1876/2333 failures and demonstrate evidence preservation.
- Compare exact and damage-aware layout with common consensus, and report truly
  damage-disabled controls separately. Exact layout with damage consensus is not
  a wholly damage-disabled control.

### Frozen development gates

Use existing runners, consumed seeds 199001/199002, depths 2x/5x/10x/20x, existing
six main and three error-only scenarios, both modes, layout and consensus.
Freeze scenario names/settings and ALL imported source hashes before the matrix.
48 main combinations/192 rows and 24 error-only combinations/96 rows are correlated
cases, not independent biological replicates.

EVERY row must have no lost previously correct coordinate, new discordance,
newly evaluated incorrect coordinate, increased mosaic/incompatible/unresolved
burden, reduced evaluated span or increased mismatch fraction. Sidecars do not
satisfy primary gates. Check ID accounting and duplicate keys before joining.

Benefit: clean unique 20x primary N50 >=1.5x baseline for EACH seed in both modes.
EACH seed must improve damage-aware layout N50 in at least one damaged scenario
without lower consensus N50. Report all other N50 regressions. Error-only controls
may be safe no-ops. Assembly time <=2x matched baseline; evaluator costs separate.
These retained engineering criteria are not literature-derived guarantees.

Decision: all gates pass → Phase C. Safety failure or insufficient benefit →
reject candidate and retain baseline. No second architecture, relaxed threshold,
fresh-seed replacement or new opportunity sweep without review. A conservative
alpha may still be packaged, explicitly without claiming the rejected benefit.

## Phase C — validate once, package, request one real run

1. Register two uninspected seeds in the existing ledger, excluding every consumed
   seed. Freeze code/settings and run ONCE with SAME gates. Failure ends promotion;
   do not tune and reuse those seeds as held-out validation.
2. Retain the sensitivity requirement before readiness claims: preregister
   realistic lengths/qualities and errors not perfectly tagged Q10. Prefer an
   existing versioned simulator, record assumptions; no simulator framework or
   fitting to test outcomes. Idealized passes remain limited evidence.
3. Harden the existing CLI: reject input/output collisions and existing
   destinations before assembly; stage outputs and publish only after validation.
   Record partial multi-file publication as failure. Preserve assembled/unresolved
   accounting. Test malformed inputs, interruptions and collisions with small
   fixtures. Incrementally hash large inputs and record full source/config hashes
   and stage time/RSS. No speculative CLI modes or unvalidated paired claims.
4. Give the user exact commands for the SAME frozen EMN001 50k cohort, baseline
   and candidate; reuse comparator results only if input/settings match. Verify
   paths and installed help; never invent them. Separate output directories and
   explicit time limits. No TAF016 audit, 100k expansion or parameter grid yet.

Retain 50k operational gates: <=120 seconds, <=1 GiB RSS, all-record output/input
N50 >=1.10 and record count/input read count <=0.90. Also report assembled-only
N50/auN/longest/length, unresolved count, aligned length/coverage, mismatch/indels
and misassemblies with identical filters/reference. Never mix all-record and
assembled-only denominators. Zero short-contig misassemblies is weak evidence;
reference mapping is a proxy, not truth.

Decision: runtime failure → profile that stage and permit only semantics-preserving
optimization checked by sequence/layout equality. Accuracy failure → no promotion.
Safe but no useful real contiguity/coverage gain → report conservative prototype
status, not a competitive breakthrough. A pass supports a documented minimal
prototype; broader claims need additional evidence.

## Findings and reporting contract

Append one section per increment to `sol_implementation_findings.md`:
`Implementer findings — model/author — date — phase`. Preserve reviewer entries.
Include hypothesis, changed files, commands/environment, source/input hashes,
tests RUN/PASSED/SKIPPED, per-case numerators/denominators, failing coordinates,
resources, PASS/REJECTED/NO-OP/NOT-RUN disposition and permitted next action.
Never rewrite a criterion after observing results.

User updates state the measured change and whether actual-file commands are
needed. Phase A: NO actual-file test yet. Begin with exact-core regressions and
the minimal kernel. Do not append another future phase to this guide.
