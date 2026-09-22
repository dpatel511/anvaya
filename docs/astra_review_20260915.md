# Assembly progress review — Astra, 2026-09-15

## Verdict

The project has made useful progress in diagnosis, reproducibility and runtime,
but has not demonstrated the requested damaged-data contiguity improvement with
held-out safety. The immediate problem is not lack of another consensus tweak.
Two containment experiments changed which primary reads could participate in
correction, and the second experiment repeated that risk despite keeping removed
reads separately. The guide accumulated conflicting active instructions instead
of ending the mechanism. This is also a reviewer/GPT-6 workflow failure, not just
Sol underperforming.

The [replacement guide](sol_implementation_guide.md) allows one bounded attempt
at a different stage: merging baseline contigs with existing raw support logic,
without modifying raw-read containment. It explicitly closes that route if the
existing evidence cannot support useful joins. It is a testable hypothesis, not
a promise of improvement. Prototype readiness remains unproven; a calendar or
percentage-complete estimate would be invented.

## Evidence and attribution

Reviewed current assembly, projection, linking, provenance, simulator, paired-gate
and regression-test paths; frozen Phase A/B/C reports; the current guide/findings;
and primary literature below. Graphify supplied navigation (1,200-token query
budget); current files were checked directly because index line numbers can be
stale. This is a targeted audit of the active mechanism and evidence, not a claim
that every repository line or every old output was independently revalidated.

No actual FASTQ assembly was run. Saved results were read, gate-contract probes
executed and the complete unit suite run. Pre-edit documentation, probe source,
probe output, dependency log and unit log are in
`results/astra-review-20260915/`. The current graph hash matches the restored
baseline `5bae98bb7dc5fb4fcedc3619907db82860c63be5d3e331c1f762395c094a828f`.
Dirty/untracked files predate this review and cannot all be attributed to one
model. Findings explicitly identify recent work as GPT-6 continuations.

| Implementation | Evidence | Assessment |
| --- | --- | --- |
| Baseline 50k operation | 43,478 contigs; input/output N50 62/65; saved run 62.65 s, 326,504 KiB | Operational/resource success; fails 1.10 N50 response criterion |
| Selected mapping proxy | 578/50,000 eligible reads; 65/70 directed expected overlaps selected | Biased 1.156% subset; not a whole-dataset recall or connectivity ceiling |
| Boundary-edge deletion | Development clean N50 132 to 131 | Failed simplification; no contiguity evidence |
| Longer-path redundancy diagnostic | No sequence-equivalent witness for the 690 residual strict physical branches examined | No-op on this case, despite a valid small algorithmic counterexample |
| Maximal nonbranching diagnostic | N50 166 versus 132, but 732 outputs from 416 reads and 314 reused nodes | Output reuse prevents accepting the N50 increase alone |
| Containment attempt 1 | Clean N50 132/135 to 2398/2382; error-only strain case loses 51 bp | Rejected; coordinate ownership did not ensure biological phase |
| Containment attempt 2 | Development passes; validation main fails 4/192 rows, error-only passes 96/96 | Correctly rejected; clean gain does not compensate for primary errors |
| Four-read reduction | Removed clean partner loses incident correction edge; damaged primary survives with T | Useful causal regression, derived from consumed validation |

The four held-out failures are two correlated scenarios at two stages, not four
independent biological failures. They nevertheless violate the declared safety
contract. In seed 20261401 at 20x, reference-0 positions 1876 and 2333 change from
concordant C to C/T discordance, while N50 falls 137 to 136. Coordinates here are
0-based; the attempt-1 lost interval is half-open `[2069, 2120)`. The main/error
gate statuses were reread directly by this review's probe script.

## Findings that affect the next implementation

### P1 — the executable acceptance gate is weaker than the written contract

`experiments/paired_coordinate_gate.py::compare_summaries` matches case keys and
selected metadata, but lacks mandatory content/evaluator identity checks and an
external expected matrix. Both sides may omit the same required case. Its benefit
condition admits any damaged-mode gain, including a one-base N50 change; it does
not enforce the guide's clean 20x 50% requirement.

The review probes return true for all three demonstrations: a 100-to-101 gain
without clean controls, unequal declared input hashes, and the same missing
declared case on both sides. These demonstrate missing enforcement, not evidence
that historical matrices actually used mismatched inputs. The old coordinate
gate DID detect the containment failures and should not be weakened.

Fix the manifest/hash/benefit contract prospectively before the next candidate.
Archive old schemas unchanged. Call the 48/24 values case combinations, not
independent cases: only two reference seeds underlie each panel.

### P1 — retaining a read is not the same as retaining its correction role

`damage_string_graph.py` performs correction-sensitive overlap decisions before
projection. Suppressing a child from topology removes relationships even if its
sequence remains in another output collection. The frozen traces show primary T
observations and unresolved C observations at the failed loci. The four-read
case isolates loss of a clean correction partner. Do not attempt another variant
of containment suppression or treat sidecar coverage as primary preservation.

The existing paired evaluator's union of primary and unresolved evidence was a
useful detector, but a primary-only N50 beside union-evidence accuracy is a mixed
reporting basis. Any new mechanism must report primary correctness and all output
categories explicitly; moving records between categories is not an improvement.

### P1 — the proposed linker is reusable but not integration-ready

`overlap_progressive_links.py::audit_raw_supported_progressive_links` already
provides exact overlap proposals, raw spanning support and rejection of molecules
supporting multiple discovered edges. Rewriting those components is unnecessary.
However, it returns projected sequences/diagnostics without a composed provenance
pool. `overlap_graph.py::_project_master_edges` has layout support that can be
reused through a narrow adapter. Preserve existing shared-call behavior.

Never reinterpret a derived pool record's molecule ID as all of its raw support.
Composition must carry raw lengths, intervals, strands, qualities and molecule
identity into the merged coordinate system. Otherwise subsequent consensus can
apply damage rates at the wrong ends or count duplicated evidence.

Before implementing the adapter, measure bridge feasibility. `_spanning_molecules`
requires reads that cross both overlap boundaries; long contig overlaps may be
impossible for short raw fragments to bridge. Raising overlap length or support
counts can make this worse. Its docstring's "unique sequence" refers to outside-
overlap geometry, not proven genomic uniqueness. Candidate-local uniqueness is
limited by sparse discovery. Use tiny exhaustive tests and abstain on competing
placements; do not claim a genome-wide repeat guarantee.

### P2 — regression evidence and independence were overstated

`tests/test_containment_survival_fixture.py` is deterministic at runtime, but its
literal reads were selected from the failed validation seed. It is not independent
validation. Delta debugging supports a deletion-minimal counterexample under the
chosen predicate, not proof of a globally smallest possible fixture. Only one of
the two failed positions was reduced this way.

The test asserts baseline path/edge counts as well as biology. Keep it as historical
protection, but add candidate tests of sequence and provenance rather than forcing
safe longer layouts to keep the same topology. Its dictionary keyed by placement
sets could overwrite duplicate groups: check uniqueness explicitly. Whole-fixture
reverse-complement and order-permutation invariance are additional missing checks.

### P2 — model, simulator and runtime claims need narrower interpretation

Graph damage compatibility uses positive profile rates categorically; changing a
positive profile to half its value does not test probabilistic graph calibration.
Consensus does use numerical likelihoods. Exact-layout runs followed by that same
consensus are not a completely damage-disabled pipeline control.

The variable-length simulator uses uniform 31–200 bp reads and marks simulated
sequencing errors at Q10 versus Q35 otherwise. This is useful for controlled tests,
but favorable to quality-aware inference and insufficient as realism validation.
Depth/scenario rows from the same reference seed are correlated. No precision or
generalization claim should use their row count as independent sample size.

The rich main evaluator consumed about 1.46 GiB in saved runs; do not attribute
that to the assembler. Keep case evidence on disk and aggregate compact summaries
if memory becomes blocking. Profile before optimizing: existing measured 25k
candidate search dominated runtime, with phase auditing another material stage.
Do not delete audit work or change graph behavior merely to improve a timing line.

## Literature: what it does and does not justify

[CarpeDeam](https://link.springer.com/article/10.1186/s13059-025-03839-5) separates
fragment-based extension from later corrected-contig merging. Its safe mode uses
extension-region evidence, and it acknowledges repeat/strain misassembly risks.
These principles motivate trying a separate merge stage while preserving raw
fragment damage coordinates; they do not validate our support rule or promise a
gain on this dataset. The paper and official repository documentation were read;
this review did not perform a pinned-commit source audit of CarpeDeam internals.
Sol must pin the [official implementation](https://github.com/LouisPwr/CarpeDeam)
and inspect its [analysis workflows](https://github.com/LouisPwr/CarpeDeamAnalysis)
before using them for matched commands. Check GPL-3.0 obligations before copying
code; conceptual inspiration and a CLI benchmark do not require a code transplant.

[RAFT](https://pubmed.ncbi.nlm.nih.gov/39406502/) addresses haplotype information
lost through contained-read removal in long-read assembly. It reinforces caution
about deleting evidence, but its operating regime differs from ultrashort damaged
reads. The [official code](https://github.com/at-cg/RAFT) is a reference, not a
validated solution to transplant here.

## Decision and verification

The new guide replaces accumulated active plans with: one contract fix; one
bounded post-layout opportunity check; one provenance-preserving candidate if
supported; one frozen development/validation sequence; then user-run matched
50k comparison. Failure redirects to a comparator-backed architecture decision,
not another correction loop. Existing correction thresholds remain fixed.

This review ran 326 tests: **325 passed, one skipped**, with no failures, under
WSL Python 3.12 and `pysam 0.24.0`. The project declares Python >=3.13; this run
does not verify that declared runtime. Historical wording "326 passed with one
skip" is incorrect when unittest reports "Ran 326 ... OK (skipped=1)". Interpret
each historical test log similarly; do not invent a corrected count without its
log. Current suite success does not prove assembler improvement.

Production code was not changed in this review. Updated directives and the audit
are documents; small gate probes live under ignored review results. No actual
FASTQ was run and nothing was staged or committed.
