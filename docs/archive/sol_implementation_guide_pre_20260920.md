# Sol implementation guide — active contract, 2026-09-15

This is the ONLY active plan. It supersedes every next-action instruction in the
[archived guide](archive/sol_implementation_guide_pre_20260915.md) and older
findings. Read `AGENTS.md`, this document and the latest reviewer entry in
[findings](sol_implementation_findings.md). See the
[review](astra_review_20260915.md) for evidence and literature limitations.
GPT-6 continuations share responsibility for the previous loop; do not attribute
all dirty files or all failures to Sol.

## Objective and standing

Improve contiguity on damaged reads while preserving correct sequence and raw
molecule provenance. A unit-suite pass, clean-data N50 gain or real-file execution
alone does not meet this objective. No next mechanism is guaranteed to succeed.

The restored `damage_string_graph.py` SHA256 is
`5bae98bb7dc5fb4fcedc3619907db82860c63be5d3e331c1f762395c094a828f`.
Freeze/hash ALL imported sources, not just this file, before experiments.
Preserve unrelated dirty workspace changes.

Accepted: 50k operation, containment performance work, coordinate safety checks,
reproducible failed candidates and a four-read regression. Unproven: improved
damaged-data contiguity passing independent safety validation. The real 50k N50
response gate failed: input 62, output 65.

Both ambiguous-containment suppression attempts are CLOSED. Attempt 1 lost a
51-bp strain interval; attempt 2 exposed incorrect primary bases at reference-0
positions 1876 and 2333 despite retaining removed children in a sidecar. Do not
revive suppression, move reads into unscored categories, or relax the safety gate.

### Current execution checkpoint — 2026-09-15

Increment 1 is COMPLETE. Increment 2 is CLOSED as a measured no-go. On seed
199001, clean unique 20x had 2,512 exact, geometrically span-able dovetails but
zero with two unique spanning molecules; N50 stayed 132. Damage unique 10x had
three supported paths and N50 145→146. This cannot meet the clean >=1.5x gate.

Do NOT execute Increment 3 or implement the provenance adapter for this linker.
Those sections remain the frozen protocol explaining why it was not authorized.
The active next action is the matched comparator and architecture decision under
"Stop and redirect." Do not rerun the opportunity check, lower support, enable
near-exact links or substitute new seeds.

### Comparator decision — 2026-09-16

The architecture review selected CarpeDeam safe mode for a matched external
comparison. It is not Anvaya's backbone or production dependency. The unchanged
internal string graph remains Anvaya's baseline. Do not add another custom
containment, linking, correction or graph-selection heuristic while the matched
comparison is open.

Use `anvaya carpedeam-assemble`, implemented in `carpedeam_backend.py`. It must:

- call a separately installed, pinned CarpeDeam executable; do not vendor or copy
  its GPL-3.0 source into Anvaya;
- verify the documented `ancient_assemble` interface, stream-validate the input,
  hash the executable/input/profile files and preserve stdout/stderr;
- force safe mode, minimum contig length 31, five raw-read iterations, ten total
  iterations, merge identity 0.99 and safe coverage five;
- refuse existing output, temporary, diagnostics and log paths; and
- treat missing, empty or malformed output FASTA as failure while preserving an
  `invalid_output` diagnostics manifest.

The next actual-data action is ONE 50k run using the already frozen cohort
`results/damage-string-graph-emn001-smoke-50000-v1/reads-50000.fq.gz` and the
same `results/raw-consensus-100k/profile` prefix used by its internal baseline.
The user runs it; the agent does not. Do not request 100k or tune the pinned
settings first.

First enforce the existing operational response gates: wall time <=120 seconds,
RSS <=1 GiB, output/input N50 >=1.10 and contig/read count <=0.90. A response pass
authorizes a matched reference evaluation; it is not an accuracy pass. Compare
the backbone output and frozen internal baseline with identical reference,
minimum-contig and MetaQUAST settings. Require no new misassembly and no worse
mismatches, indels, aligned length or evaluated coverage before calling the
backbone an improvement. Record every denominator. If operational response fails,
close this pinned candidate. If accuracy fails, retain it only as a comparator;
do not repair it with an unregistered Anvaya consensus pass.

CarpeDeam output currently exposes no raw-read placement sidecar. Therefore do
not claim Anvaya molecule provenance, coordinate safety, strain retention or a
validated post-assembly damage correction from its FASTA. A later adapter is
not authorized as the product direction: Anvaya is intended to remain an
independent competing assembler. Use comparator stage counts and output behavior
to specify the next internal architecture, without copying implementation code or
reconstructing placement provenance from comparator FASTA.

### Matched-result decision — 2026-09-16

The 50k reference comparison is COMPLETE. CarpeDeam is not the next Anvaya
architecture. Its longer N50 came from a much smaller extended-only output and
did not improve reference recovery or sequence accuracy. Do not integrate it,
copy its implementation, tune Anvaya against its FASTA, or add another consensus
pass.

The next bounded task is an output-contract diagnostic using existing immutable
Anvaya placements. Split the current output by number of DISTINCT raw molecule
IDs: assembled records have at least two; singleton/unresolved records go to a
separate FASTA. Every original output record must appear exactly once across the
two files and sequences must remain byte-identical. Molecule counting uses IDs
only; placement coordinates in the TSV are 0-based with half-open read intervals
but are not transformed by this split.

Use `experiments/split_assembled_output.py`; do not recreate this partition with
ad hoc FASTA filtering. It rejects ungrouped placement records, missing/duplicate
FASTA IDs and reused output paths, and records input/output hashes and accounting.

On the frozen 50k output, the preread diagnostic gives 5,125 assembled records,
478,252 bases, N50 101 and longest 287, versus CarpeDeam's 4,962 records,
510,650 bases, N50 112 and longest 346. Run one matched MetaQUAST comparison of
these similarly selected outputs. Require the split Anvaya set to retain its
lower mismatch/indel rates and competitive aligned length/genome fraction. This
is evaluation normalization, not a contiguity improvement.

If that gate passes, implement the output split as an explicit CLI contract:
primary assembled FASTA plus an unassembled-fragment sidecar, while preserving
the current all-record output for compatibility. Then return to internal layout
work. The only authorized algorithmic target is branch resolution/path emission
inside `_project_master_edges`; candidate discovery already supplied 4,609
physical edges at 50k and the graph rejected 628 physical edges at ambiguous
ends. Do not revisit containment suppression, exact post-layout linking, anchor
thresholds, damage consensus or CarpeDeam integration.

Any later branch-resolution candidate must be based on direct raw molecule paths,
preserve alternative branches instead of forcing one winner, and pass the existing
per-coordinate synthetic safety gate before another actual-data run. If short
reads do not span a branch with independent evidence, abstain; no algorithm can
recover linkage absent from the data.

The matched assembled-only gate PASSED and the explicit CLI contract is now
implemented. `damage-assemble` retains required all-record `--output-before` and
`--output-consensus`; requesting `--output-assembled` requires the paired
`--output-unresolved`. The two optional FASTAs partition consensus records using
at least two distinct, uniquely placed raw molecules and retain the original
`unitig_N` identifiers. Diagnostics record counts/bases for both partitions.

Before branch work, run one 50k contract check and compare both optional FASTAs
byte-for-byte with the already validated diagnostic split. This rerun validates
plumbing only; do not repeat MetaQUAST if both hashes match. A mismatch blocks
branch work and must be explained from eligibility semantics rather than hidden.

The 50k contract check PASSED byte-for-byte on 2026-09-16. All-record before and
consensus FASTAs retained their historical hash; assembled and unresolved FASTAs
matched the diagnostic hashes. Output normalization is CLOSED. Do not rerun this
dataset or MetaQUAST for the next feasibility audit.

The next phase is diagnostic only: measure branch-transition evidence after
transitive reduction and before `_project_master_edges` rejects ambiguous ends.
For each ambiguous oriented node, enumerate incoming/outgoing edge pairs and test
whether an independent raw molecule spans both junctions with one consistent
orientation and coordinate system. Endpoint molecules and duplicate placements
do not count as independent support. Conflicting pairings abstain. Report branch
nodes, possible transitions, transitions with >=1 and >=2 independent molecules,
uniquely resolved transitions, affected physical nodes and projection-only N50.

Use the already consumed development seeds first; do not touch validation seeds or
the 50k FASTQ. If no transition has two independent spanning molecules, close this
branch-resolution route without production code. If opportunities exist, freeze
one path-emission candidate that never duplicates a raw molecule across primary
outputs and preserves every unsupported alternative as unresolved evidence.

### Branch-transition decision — 2026-09-18

The diagnostic is COMPLETE and this direct branch-transition route is CLOSED.
On consumed development seed 199001, clean unique 20x had 1,661 possible
transitions, 41 with at least two independent non-endpoint molecules and 41
uniquely resolved transitions. Damage unique 10x had 514 possible transitions,
14 with at least two molecules and 14 uniquely resolved transitions. Conflicting
molecule pairings and duplicate placements abstained. The audit's zero-transition
projection reproduced baseline N50 exactly.

An optimistic no-node-reuse projection improved clean N50 only 132 to 152
(1.15x) and damaged N50 only 145 to 156 (1.08x). It therefore cannot meet the
preregistered clean 20x >=1.5x benefit gate even before coordinate-safety
validation. Do not implement, tune or lower support for this path-emission rule,
and do not run it on validation seeds or actual FASTQ.

The next decision must be architectural rather than another local correction or
graph heuristic. Preserve the current internal assembler and assembled/unresolved
output contract as the prototype baseline. Before authorizing more code, compare
standard overlap-layout options that can use short-read pairing or reference-free
long-range linkage while retaining molecule provenance. Any proposal must state
where information beyond these single-end short reads comes from; without added
linkage, the measured branch ceiling is the stopping evidence.

### Paired-end architecture check — 2026-09-18

Two existing paired mechanisms are CLOSED. Do not rerun or tune them:

- Merging overlapping mates before the current string graph failed the clean
  consumed-seed control. On seed 199001 at 20x, the unmerged exact paired reads
  produced N50 2,399 while 138/356 merged pairs produced N50 99. This is a hard
  clean-data regression, not a damage-threshold issue.
- Post-contig paired scaffolding found no usable linkage on the existing empirical
  TAF016 runs. At 100k pairs it had 591 uniquely mapped pairs but zero cross-contig
  pairs and zero supported links; at 500k it again emitted zero supported links.
  Its insert calibration was invalid on this short-fragment library.

The only paired hypothesis still open is diagnostic graph-level mate threading:
measure whether independent mate pairs select one transition at the already
identified reduced-graph branches without first merging mates or requiring both
mates to map to completed contigs. This requires genuine synchronized R1/R2 data;
the EMN001 merged-fragment benchmark cannot test it. Before production code, run
one bounded opportunity audit on an ordered TAF016 subset. Require at least two
independent molecules per transition, conflicting pairings to abstain, exact
pair-name/count synchronization, no use of reference truth for selection, and a
projection-only clean benefit capable of meeting the existing 1.5x gate. A no-op
or sub-gate projection closes paired threading without threshold variants.

The bounded audit is implemented in `experiments/paired_branch_opportunity.py`.
It is exact/reference-free and changes no assembly output. It validates synchronized
pair names and counts, excludes endpoint molecules and duplicate placements,
requires inward opposite orientations across both junctions, and removes any
molecule supporting competing transitions. Eight deterministic anchors per read
bound the index; candidate placements are then verified base-for-base. Focused
and full synthetic tests pass. The next action is one user-run 25,000-pair ordered
TAF016 audit. Do not increase the prefix or add damage-tolerant placement until
this exact opportunity result is reviewed.

### Paired-threading decision — 2026-09-19

The 25,000-pair TAF016 audit is COMPLETE and graph-level paired threading is
CLOSED. The unchanged exact graph had 2,769 physical branch nodes and 30,114
possible transitions, but zero transitions had even one eligible independent
pair. Consequently zero reached two-pair support, zero were uniquely resolved,
and projection-only N50 remained 76. The zero-transition projector independently
reproduced baseline N50 76.

This is not permission to relax placement geometry, include endpoint molecules,
lower support, add damage-tolerant mapping or run larger prefixes. The synthetic
positive control demonstrated that the audit detects a uniquely supported
two-pair transition, while input synchronization and the projection invariant
passed on the empirical run. The empirical zero therefore closes this frozen
mechanism rather than starting a parameter search.

The audit also failed the engineering budget: 44:25 wall time and 1,056,460 KiB
peak RSS for 50,000 reads. Do not optimize a biologically empty mechanism.
Retain the current single-fragment assembler and explicit assembled/unresolved
outputs as the prototype baseline. Further contiguity claims require new physical
linkage data (for example capture, linked-read or long-read evidence), not another
rule over this short-read graph.

## Non-negotiable execution rules

1. One hypothesis, one production mechanism, one frozen candidate. No simultaneous
   anchor, overlap, support, consensus and graph changes; no extra correction pass.
2. Record the predicted measurable change and rejection condition before coding.
   Diagnostic truth may enter evaluation only, never assembly decisions.
3. Coordinates are 0-based, intervals half-open. Preserve original molecule ends,
   strands and qualities. Corrected contig ends are not raw damage coordinates.
4. A merged output may replace its inputs only with checked provenance composition.
   No silent evidence loss, duplicate-path inflation, or N50 gain by changing the
   output denominator. Other evidence categories never compensate for primary loss.
5. Preserve graph projection defaults, containment decisions and existing damage,
   quality and search thresholds. Do not opportunistically refactor shared code.
6. One bounded feasibility check, one development candidate, at most one held-out
   run. No-op closes this approach; safety failure rejects the candidate. No seed
   replacement, threshold sweep or automatic new diagnostic phase.
7. Agent runs unit/synthetic tests; user runs actual FASTQ. Provide verified exact
   commands when ready. Do not stage/commit without a fresh request. Keep local
   notes out of commits; use standard branch prefixes, never `codex/`.

## Increment 1 — make the acceptance contract executable

Limit scope to `experiments/paired_coordinate_gate.py`,
`experiments/variable_length_validation.py` and focused tests. This is one small
increment, not a new benchmark framework. Do not rerun rejected matrices.

Reviewer probes in `results/astra-review-20260915/gate-probes.json` show that the
current gate accepts a one-case 1% N50 gain with no clean controls, ignores unequal
input hashes, and ignores a declared case omitted from both summaries.

- Require a versioned expected manifest of seeds, scenarios, depths, modes and
  stages declared BEFORE execution. Reject missing/extra/duplicate members against
  it, including a case missing from both baseline and candidate.
- Hash ordered read IDs, sequences, qualities, molecule IDs, reference sequences
  and supplied damage profiles with deterministic unambiguous encoding. Compare
  actual content hashes, not only counts/configuration. Require generator and
  evaluator identity and record all source hashes. Assembler hashes may differ;
  inputs, reference, generator and evaluator must match. Fail closed on omissions.
- Preserve historical artifacts and label them legacy-schema; do not retroactively
  claim they passed the stronger contract. Validate reference identity before
  comparing integer coordinate keys.
- Replace `independent_cases` with `case_combinations`; separately report distinct
  reference seeds. Depths/scenarios/modes/stages sharing a seed are correlated.
- Implement ALL benefit requirements below, including clean 20x >=1.5x. Keep the
  existing per-coordinate safety failures unchanged.

Verification: tests reject equal counts with unequal content, missing hashes,
wrong reference/evaluator, jointly omitted cases, duplicates and insufficient
clean gain. Include a fully declared passing fixture and old regression tests.
Record tests RUN, PASSED and SKIPPED separately. Then move on.

## Increment 2 — one post-layout merging feasibility check

Hypothesis: compatible supported extensions remain between baseline contigs even
though changing raw-read containment topology was unsafe. This differs from old
progressive-link experiments only if this new substrate has measured opportunities.

Reuse `audit_raw_supported_progressive_links` in
`src/anvaya/overlap_progressive_links.py` on baseline damage-string-graph pools.
Use exact overlaps only (`allow_near_exact=False`), existing overlap settings and
minimum raw support 2. No RY rescue or threshold search. Literature motivates a
separate stage, not a guarantee for these defaults.

Run ONLY development seed 199001 clean unique 20x and damage unique 10x initially.
Report input contig lengths, candidate physical joins, exact joins, overlap
lengths, geometrically span-able joins, raw-supported joins before/after competing
edge rejection, accepted paths and projection-only N50. Mark coordinate safety
`not_evaluable_without_composed_raw_placements` until the adapter exists; never
score projected FASTA as if its records retained the baseline pool's placements.

The existing bridge rule requires a raw read to extend past BOTH overlap
boundaries. A read cannot bridge an overlap as long as itself. Quantify this
before blaming support thresholds. Flanks outside an overlap are not necessarily
genomically unique; uniqueness among sparse discovered candidates is not proof
that all repetitive alternatives were considered.

If no useful joins survive, or this exact mechanism cannot meet the clean gain
criterion, CLOSE it and follow the redirect below. Do not add more diagnostic
phases. If opportunities exist, implement the provenance adapter. Counts alone
do not establish safe joins or successful assembly.

## Increment 3 — compose provenance and evaluate the frozen candidate

The existing linker returns projected sequences/diagnostics, NOT a complete pool
with composed raw placements. Do not plug its FASTA into consensus. Implement a
narrow adapter using existing projection layout support; preserve current callers'
behavior by default.

Order: baseline raw graph/layout -> exact raw-supported contig merging -> ONE
raw-evidence consensus pass. Compare both layout and consensus using the same
evaluator. No pre-polishing to manufacture exact joins in this candidate.

- Carry source-contig identity, orientation and path offset. Compose each raw
  interval/strand through the layout, including reverse complements and unequal
  lengths. Validate bounds and retained qualities against immutable raw evidence.
- Do not use a derived pool slot's `molecule_id` as all its raw support. Count each
  raw molecule once per locus. Copies of a derived contig are not new molecules.
- Conflicting offsets/orientations or competing joins require abstention. Do not
  infer strain phase from one shared matching read or discard competing evidence.
- Emit unmerged inputs unchanged. Reject cycles/ambiguous paths. Audit ownership
  before/after: no disappearing molecules or newly conflicting primary placements.

Required focused tests: forward/reverse joins, unequal lengths, permutations,
nested provenance, duplicate molecules, conflicting offsets, cycles, repeat and
shared-strain alternatives, no spanning read, and raw-end damage after merging.
For tiny fixtures compare indexed discovery with exhaustive pair enumeration to
expose hidden alternatives. The exhaustive oracle is test-only.

Keep the four-read survival regression as historical protection. It was reduced
from consumed validation: deterministic and RNG-free, NOT independent validation.
Add biological-invariant tests under reverse complement/permutation. Do not force
a safe candidate to preserve exactly two paths or baseline edge counts. Check
placement-group uniqueness before constructing dictionaries that could overwrite
duplicate groups.

## Development gates — safety AND contiguity

Freeze the candidate before the full matrix: seeds 199001/199002; depths
2x/5x/10x/20x; existing six main scenarios and three error-only scenarios; exact
and damage-aware modes, layout and consensus stages. Declare exact scenario names
in the manifest. Totals: 48 main combinations/192 rows and 24 error-only
combinations/96 rows, but only TWO reference seeds per panel.

EVERY row must pass current coordinate safety: no lost previously correct
coordinate, new discordance, newly evaluated incorrect coordinate, increased
mosaic/incompatible/unresolved burden, reduced evaluated span or increased
mismatch fraction. Fewer mismatches over less evaluated sequence are not an
accuracy gain. No averaging away strain-specific losses.

Prospective benefit requirements for this CONTIGUITY candidate:

1. Clean unique 20x primary N50 >=1.5x baseline for EACH seed in both layout modes.
2. EACH seed has at least one damaged scenario with strictly improved damage-aware
   layout N50 and no lower consensus N50. Correction-only improvement is insufficient.
   Report all N50 losses elsewhere, even when they are not safety failures.
3. Error-only controls may be safe no-ops. Exact layout followed by common damage
   consensus is not a fully damage-disabled control. Within-mode improvement does
   not prove the damage model superior to exact assembly.
4. Record assembler time/RSS separately from evaluator costs. Engineering budget
   for this increment: candidate assembly time <=2x matched baseline. This limit
   is not a biological or literature-derived criterion.

A failed gate rejects the submitted candidate. Preserve useful failure evidence,
restore only candidate edits and redirect. No tuning on failures.

## Validation and actual FASTQ checkpoint

After development passes, register two uninspected seeds in the existing ledger,
excluding every development/diagnostic/consumed seed. Freeze all code, manifests
and settings; run once with the SAME gates. No replacement validation seeds.

The current uniform length distribution and Q10 errors/Q35 other bases are limited
controls. Before claiming readiness, preregister a separate sensitivity panel with
empirical lengths/qualities and errors not perfectly marked by low quality. Prefer
an existing simulator, record its version/assumptions, and do not build a simulator
framework or tune on this panel. A pass on idealized data is not real-data validation.

After validation, supply user commands for the SAME 50k FASTQ cohort on baseline,
candidate and a versioned CarpeDeam safe-mode comparator. Verify installed CLI
help first. Include input/source hashes, environment, separate output directories,
time/RSS, diagnostics, damage-profile conventions and identical evaluation
filters/reference/settings. Reuse verified runners; never invent input paths.

Keep the existing 50k limits: 120 seconds, 1 GiB RSS, output/input N50 >=1.10
and output contig count/input read count <=0.90.
Also compare matched baseline/comparator N50, auN, longest, aligned length,
coverage, mismatch/indel and misassembly metrics. Operational response is not
correctness; zero QUAST misassemblies in tiny contigs is weak evidence. The mapping
proxy is not truth. No new 100k request before assessing this checkpoint.

## Stop and redirect

This mechanism gets ONE development candidate and at most ONE held-out attempt.
Insufficient linking opportunities or safety failure closes it. Retain baseline
and complete a matched comparator assessment. Identify discovery, support
geometry, branching or error/phase preservation as measured limitations; do not
infer the full dataset's connectivity ceiling from the selected mapping proxy.

The next review must choose a concrete architecture: retain the custom assembler
only with measured benefits; integrate a vetted assembly backbone while retaining
separately validated aDNA handling; or narrow the prototype scope. Another custom
heuristic is NOT automatically authorized. Check licensing before copying code.

## Required primary literature

- [CarpeDeam paper](https://link.springer.com/article/10.1186/s13059-025-03839-5):
  staged extension/contig merging and raw-end damage treatment motivate this
  experiment. Safe mode examines extension evidence; it does not guarantee
  strain-safe joins. Do not blindly copy its support thresholds.
- [Official implementation](https://github.com/LouisPwr/CarpeDeam) and
  [benchmark workflows](https://github.com/LouisPwr/CarpeDeamAnalysis): inspect a
  pinned commit before relying on code behavior or commands. CarpeDeam is GPL-3.0.
- [RAFT paper](https://pubmed.ncbi.nlm.nih.gov/39406502/) and
  [code](https://github.com/at-cg/RAFT): contained reads can carry haplotype
  information. Its long-read evidence warrants caution, not direct short-aDNA
  algorithm transplantation.

## Findings contract

Reviewer instructions live here; implementer observations append under an
author/model/date heading in `sol_implementation_findings.md`. Never rewrite a
criterion after seeing results. Every increment records hypothesis, source/input
hashes, commands/environment, tests RUN/PASSED/SKIPPED, per-case deltas, failing
coordinates and denominator changes, assembler resources, and disposition
PASS/REJECTED/NO-OP/NOT-RUN. Select the next action from this guide.

End user updates with the measured delta and whether actual-file commands are
needed. Rejected mechanisms are useful evidence, not improved contiguity. Do not
say a phase is almost ready repeatedly without running its declared gate.
