# Implementation guide for GPT-5.6 Sol

> ARCHIVED 2026-09-15. Everything below is historical, including statements that
> call themselves active. Follow [the current guide](../sol_implementation_guide.md).

Read the 2026-09-13 review below first. It is the active execution plan.
Older checkpoints and phases are retained as history; they do not authorize the
latest proposed maximum-overlap heuristic. Section 14 remains the findings format.

## Astra review — 2026-09-13: active instructions

**Standing:** a bounded overlap assembler exists and real-input operation reached
50,000 reads. Improved damaged-data contiguity with independent safety validation
is still unproven. Sol and subsequent GPT-6 continuations built useful tooling,
preserved rejected candidates and improved containment runtime. Recent causal and
literature interpretations overstate what those measurements establish. Attribution
must follow each findings entry, not assume every dirty file was written by Sol.

Detailed review: [2026-09-13 audit](astra_review_20260913.md). Reproducible probes,
source hashes, saved-output hashes and pre-edit document snapshots are in
`results/astra-review-20260913/`. No actual FASTQ assembly was run by the reviewer.

### Corrections that govern the next implementation

1. The 50k graph has **4,609 input physical edges, 213 transitive removals,
   4,396 reduced edges, 628 branch rejections and 3,768 reciprocal edges**.
   Survival is 81.75%, not 86.4%. High survival among accepted edges does not
   measure true-overlap recall or establish the cause of fragmentation.
2. The exact-cohort proxy includes only 578/50,000 reads (1.156%). Its 65/70
   selected directed overlaps and 178-bp maximum component apply only to that
   selected subgraph. Excluded reads can connect its components. Neither the full
   dataset's connectivity ceiling nor its taxonomic ordering is established.
   Keep the failed 50k response gate; do not call the dataset intrinsically
   unassemblable or infer that whole-input recall can improve by at most 7.14%.
3. Clean-ladder N50 deterioration is real. Its cause is not isolated to the
   degree-one path rule: boundary containments enter the extension graph, and
   reduction tests two-hop witnesses only. The review's four-read, strict-dovetail
   probe gives 120/100/100-bp outputs when a redundant three-hop-reachable edge
   remains; its chain control gives one 160-bp contig with the SAME path rule.
   This proves a reduction limitation, not its share of the ladder's losses.
4. **Suspend reciprocal maximum-overlap selection.** A longest-overlap choice is
   a heuristic, not a strain-safety guarantee or a general string-graph theorem.
   Do not swap heuristics until residual branches are classified causally.
5. The graph uses damage rates categorically (`rate > 0`) for compatibility;
   half-strength profiles with the same support do not test graph calibration.
   Consensus uses numerical likelihoods, but cannot create long layout paths.
   Keep one consensus pass and current thresholds during topology work.

### Phase A — isolate the failing graph invariant, then stop auditing

Use current source as a frozen baseline. First promote the four-read probe into
a focused regression test. It must reproduce the 160-bp sequence and raw placement
coordinates under a sequence-equivalent path, including reverse complements and
unequal read lengths. All coordinates are 0-based, half-open. Add negatives for
different correction spellings, conflicting placement offsets, real branches,
cycles and repeated sequence. Do not assume any alternate walk is safe to remove.

On existing development seed 199001 clean 20x ONLY, capture the graph immediately
before projection. Classify residual branch edges into strict dovetails,
boundary containments and multi-offset/ownership ambiguity; find a bounded
sequence-equivalent alternate-path witness for redundant strict edges. Report
physical and oriented counts separately. Diagnostic truth may label joins but
must never enter production assembly decisions.

Perform diagnostic-only ablations separately: boundary relations represented
outside the dovetail graph; verified redundant strict edges removed. Preserve
all original reads and raw evidence. Compare lost/recovered correct junctions,
output sequence and coordinates, not only degree counts. Neither ablation is
automatically an accepted implementation. Multi-parent containment is not proof
of cross-strain ambiguity, and deleting contained reads can erase connectivity.

**Decision after this bounded trace:**

| Measured cause | Implement one next change | Stop condition |
|---|---|---|
| Redundant strict edges block correct paths | Sequence/offset-preserving reduction with a valid retained witness | No witness or conflicting spelling: retain edge |
| Boundary relations distort degrees or witnesses | Separate containment relations from strict extensions, preserving evidence and ownership | Damaged/strain regression: reject, do not delete more children |
| Both contribute | Repair the first demonstrated invariant and measure its isolated delta | Do not combine policy changes in one candidate |
| Remaining branches are genuine alternative paths | Keep alternatives; request review of unitig emission/selection semantics | No automatic maximum-overlap tie-break |
| No useful causal delta | Record the no-op and unresolved cause | No new seed matrix or proxy rerun to search for a pass |

Do not replace two-hop reduction with unbounded all-path enumeration. A bounded
oracle is acceptable for tiny tests. Production search must have measured cost,
preserve reachability through retained witnesses, and handle correction-bearing
paths and reverse-edge symmetry. `_project_master_edges` is shared: exercise its
existing callers rather than assuming the damage assembler is its only user.

### Phase B — one topology repair, case-level acceptance

**Phase A result, 2026-09-13:** complete on clean development seed 199001 at 20x.
Boundary removal regressed N50 132→131. No one of 690 residual strict branch
physical edges had a bounded sequence-equivalent witness, so the verified reducer
was a no-op. Diagnostic maximal non-branching paths reached N50 166 (+25.8%, below
the frozen +50% gate) but emitted 732 contigs from 416 reads and reused 314 nodes
across paths. It is rejected because output/provenance expansion is uncontrolled.
Do not implement any of these three candidates in production.

The measured remaining target is narrower: 603/690 residual truth-labelled branch
edges touch ambiguously contained nodes. Earlier containment suppression improved
clean contiguity but failed uneven-strain validation through nonlocal degree changes.
The next Phase B proposal must represent coordinate-compatible containment
alternatives without suppressing a node before reduction, choosing the longest
overlap, or copying every branch path. It must first specify how raw placements and
independent molecules remain unambiguous when one node has several output owners.
If that invariant cannot be expressed in the current fixed-slot sequence pool,
stop and propose the minimal representation change before altering graph policy.
Do not run another ladder until the tiny clean, repeat and two-strain ownership
fixtures pass. Evidence: `results/phase-a-invariant-trace-dev-v3/`.

Freeze baseline/candidate source and imported module paths. Test the minimal
regression first, then all unit tests with pysam present. The review ran 321 tests:
320 passed and one external-aligner integration test skipped under Python 3.12.
The package declares Python >=3.13; a supported-interpreter suite is still required
before a readiness claim. Missing dependencies are a reported gap, not a pass.

Use the existing 199001/199002 coverage ladder as development evidence. Reuse the
same input realization per baseline/candidate. Preserve six scenarios, both modes,
and layout/consensus stages; run the separate error-only panel too. Before scoring,
assert unique/equal keys `(seed, scenario, target_depth, mode, stage)`, matching
input/settings hashes, and complete denominators. Adapt the existing paired gate
to depth-aware identities if needed; do not silently overwrite repeated cases.

Freeze these gates BEFORE the candidate matrix:

- Fix the minimal sequence/placement regression. Improve clean 20x N50 by at
  least 50% on EACH development seed, with no clean incorrect/conflicting output.
- No previously correct uniquely evaluated coordinate becomes incorrect,
  discordant, unresolved or uncovered. No new demonstrated incompatible join or
  diagnostic mosaic; report newly recovered incorrect/conflicting regions too.
- Preserve unique correct recovery. Report unresolved output and duplicate
  coverage separately. Apply the existing paired coordinate contract as well as
  per-case error counts; averages and altered denominators cannot hide failures.
- At least one damaged case must improve truth-compatible contiguity or unique
  correct span, with every repeat/strain case passing safety. A clean-only success
  is a limited graph repair, not demonstrated aDNA superiority.
- A zero-action candidate fails the benefit gate. A failure is preserved with its
  causal trace and rejected; do not relax thresholds or reinterpret the gate.

Report exact-mode versus damage-mode differences separately from candidate versus
baseline differences. Existing damage-aware strain uncertainty is a baseline
limitation; matching it is not proof of strain safety. After at most two focused
mechanism failures, stop and request review with the smallest counterexample.

**Phase B attempt 1, 2026-09-13 — rejected and reverted.** The candidate treated
multiple exact containment parents as one locus when all parents implied one
0-based graph coordinate, the component spelling was compatible, no robust phase
variable was present, and provisional projection left one output coordinate. It
then retired the child and attached its raw placement to that output. The focused
clean/repeat/strain fixtures and all 333 unit tests passed under Python 3.12.3 with
`pysam 0.24.0` (one external-aligner test skipped). The 48-case main development
ladder passed the paired-coordinate gate. Clean 20x N50 changed 132→2398 and
135→2382, with one output contig for each seed; at least one damaged case improved.

The separate depth-matched error-only panel rejected the mechanism. Seed 199001,
2x 80:20 strains lost a previously correct minor-strain interval `[2069, 2120)`
on reference 1 in exact and damage-aware layout; consensus repeated the same loss.
The 51-bp child originated from reference 1 but was exactly contained by two
reference-0 parents at a common graph coordinate. The component had zero robust
phase-variable sites. Provisional projection supplied one output owner, but that
owner was a 360-bp cross-reference candidate with five unexplained mismatches.
Thus graph-coordinate uniqueness and output-coordinate uniqueness do not establish
biological evidence ownership when the discriminating allele lies outside the
contained fragment. No truth-blind threshold can recover absent linkage.

The production baseline is restored at `damage_string_graph.py` SHA256
`5bae98bb7dc5fb4fcedc3619907db82860c63be5d3e331c1f762395c094a828f`.
Evidence and the rejected source are frozen in
`results/phase-b-coordinate-ownership-dev-v1/`. This counts as one Phase B
mechanism failure. Do not restore or tune this candidate.

Before a second topology attempt, make the minimal representation proposal that
separates primary layout from unresolved containment evidence. An ambiguous child
may cease contributing extension-graph degree only if its immutable raw evidence
retains an explicit, nonduplicated unresolved owner; it cannot be assigned to a
primary contig, copied to every parent, or silently discarded. Define how primary
FASTA, unresolved alternatives, consensus eligibility, N50, and paired-coordinate
coverage are reported before changing `SequenceRecord` or graph policy. The
51-bp counterexample above, the four-read clean recovery, repeated offsets, nested
containment, reverse orientation, and a damaged C/T or G/A case are mandatory
representation tests. If primary-versus-unresolved output cannot satisfy both the
clean benefit gate and the existing coordinate-preservation gate without using
truth, stop this containment line rather than introduce a heuristic cutoff.

**Phase B attempt 2 contract, frozen before implementation on 2026-09-14.** Extend
`ProgressiveSequencePool` additively with an unresolved-containment record tuple;
do not change `SequenceRecord` states or the meaning of `active_derived`. Each
unresolved record owns disjoint immutable raw placements exactly once. It is not
placed in a primary layout and is excluded from consensus. Synthetic scoring uses
the union of primary and unresolved records for coordinate/error preservation,
while `contigs`, `bases`, N50 and longest describe primary layouts only and report
unresolved counts/bases separately. `damage-assemble` must write a separate FASTA
when unresolved records exist; absence of that output path is an error, not silent
discard. Minimum primary-contig length filtering does not delete unresolved evidence.

Reuse attempt 1's coordinate/spelling/phase eligibility and provisional unique
output-coordinate test solely to remove a child from topology. Do not attach that
child's raw placements to the primary owner. Unique single-parent containments keep
the established behavior. Require the mandatory representation fixtures above,
the dependency-complete unit suite, then rerun the exact frozen main and error-only
199001/199002 depth ladders. Attempt 2 passes only under the unchanged Phase B
benefit and zero-regression gates. Report primary contiguity and unresolved evidence
separately; never combine their lengths into a favorable N50.

**Phase B attempt 2 result, 2026-09-14 — retained for independent validation.**
The additive pool representation and unresolved sidecar contract were implemented
as frozen above. Coordinate-compatible children affect topology but their raw
placements are absent from primary layouts and consensus. The exact 51-bp
counterexample remains a separate unresolved record and all positions in
`[2069, 2120)` stay concordant correct in the combined evidence audit.

The dependency-complete suite passed 332 tests with one external-aligner skip under
Python 3.12 and `pysam 0.24.0`. The paired 199001/199002 main panel passed all 48
independent cases and 192 mode/stage rows; the error-only panel passed all 24 cases
and 96 rows. Clean 20x primary N50 changed 132→2398 and 135→2382 in both modes,
with 272/279 unresolved records reported separately and no evaluated-coordinate
loss. Damage-aware benefit occurred in low-depth damage and damage-plus-error cases.
Across the main panel, 2,170 coordinate-compatible containments were proposed,
1,750 had one projected owner and became separately owned unresolved evidence, and
420 abstained. Candidate and baseline main runs used 33.66/32.28 seconds and
1,463,332/1,468,644 KiB peak RSS; evaluator evidence dominates memory. Evidence is
frozen in `results/phase-b-unresolved-containment-dev-v1/`.

Phase B development gates are satisfied. Do not tune this mechanism before Phase C.
Register genuinely uninspected validation seeds, freeze this exact source and the
primary-versus-unresolved evaluator contract, and run the same main and error-only
panels once. Report primary assembly and unresolved evidence separately. A fresh
failure rejects the candidate without seed replacement or threshold changes.

**Phase C validation protocol, preregistered before execution on 2026-09-14.**
Seeds `20261401` and `20261402` are reserved as the only fresh validation seeds.
Run the unchanged six-scenario main panel and three-scenario error-only panel at
2x, 5x, 10x and 20x target coverage for the frozen Phase A baseline and Phase B
attempt 2 candidate. Apply the unchanged paired-coordinate gate separately to
both panels. The main panel must have no safety failures and retain at least one
safe damage-aware N50 or correct-coordinate benefit; the error-only panel must
have no safety failures. Execute each matrix once. Any failed gate rejects the
candidate without replacing seeds, changing thresholds, or editing the mechanism.
All coordinates are 0-based and intervals are half-open.

**Phase C result, 2026-09-14 — rejected and restored.** The error-only panel passed
all 24 independent cases and 96 mode/stage rows. The main panel failed four of 192
rows, all for validation seed `20261401`, damage-aware mode, 20x `damage_unique`
and `profile_mismatch_unique`, at both layout and consensus. Previously concordant
correct reference-0 positions 1876 and 2333 became C/T-discordant; evaluated
coordinate coverage decreased, discordant bases and mismatch rate increased, and
N50 changed 137→136. These are 0-based positions.

The clean 20x benefit reproduced: primary N50 changed 135→2381 and 137→2400 for
the two seeds, with unchanged evaluated spans. That benefit cannot override the
per-coordinate safety failure. The failing T observations were primary records,
not unresolved sidecar records, showing that removing ambiguous containments from
topology can change which damaged reads survive as primary outputs even when the
removed evidence is preserved separately. The source, matrices, gates, resource
logs, and reconstructed realized-coverage ledger are frozen in
`results/phase-c-unresolved-containment-validation-v1/`. Seeds `20261401` and
`20261402` are consumed validation seeds. Attempt 2 is rejected without tuning or
replacement seeds, and production was restored to Phase A SHA256
`5bae98bb7dc5fb4fcedc3619907db82860c63be5d3e331c1f762395c094a828f`.

Stop the ambiguous-containment suppression line. A 2026-09-15 delta reduction
produced a four-read seed-independent fixture. One clean 138-bp read is exactly
contained by more than one longer representation and is also the only accepted
damage-correction partner for a 136-bp read carrying terminal C→T damage. The
Phase A baseline retains two physical edges/two paths and corrects one overlap
base. Attempt 2 retires the clean partner, leaving one edge/one path, zero corrected
bases, and the damaged read as a primary T-containing output. The permanent test
`test_containment_survival_fixture.py` locks the safe baseline sequence and raw
placement membership. The literal reads and frozen baseline/candidate comparison
are in `reduced-survival-failure.json`.

This is the second focused Phase B mechanism failure, so the Phase B stop rule now
applies. Do not implement a third topology heuristic or open another seed matrix.
The review question is whether to open a separate phase for diagnostic-only local
path emission that leaves containment membership, primary-output eligibility and
the four-read fixture unchanged. Such a phase must preregister a direct raw-molecule
join-support contract and first measure actionable opportunities on existing
development cases. No sidecar ownership scheme or containment cutoff is allowed.

### Phase C — independent validation, then user-run comparison

Before fresh validation, maintain a seed-role ledger including development
199001/199002 and all inspected historical seeds. Current validation checks do not
exclude every development seed. Record actual sampled bases/coverage per strain;
the 80:20 label is target allocation, not guaranteed realized coverage. Current
strain ladders are not nested across depths because shared RNG consumption shifts
minor-strain sampling. Do not call them paired depth experiments. Do not rewrite
old fixtures or results; version any revised fixture and freeze a new baseline.

Validate once on genuinely uninspected seeds after Phase B passes. Keep repeats,
uneven strains, damage, error-only and mixed-error cases. Uniform 31–200-bp reads,
Q10-labelled substitution errors, and two synthetic strains do not reproduce the
real EMN001 distribution. Add a separate length/quality-matched synthetic panel
for transfer testing; no need to build a new simulator framework. Numerical
profile calibration needs a later dedicated test, not this topology increment.

Only after benefit and safety pass, give the user the exact command for a matched
50k comparison using the already frozen FASTQ, profile, minimum output length and
new output directory. Read current CLI help and verify every path before providing
it; do not invent a command here. The user runs actual files. Compare baseline
versus candidate N50/NA50, aligned span, mismatch/indel numerators and denominators,
unaligned output, duplication, misassemblies, runtime and RSS. Reference mapping is
conditional evidence. Byte-identical FASTA means no extra QUAST rerun is needed.
Do not escalate to 100k merely because the command succeeds.

### Optimization and reporting boundaries

The per-node `list(range(len(reads)))` and node-ID equivalent create avoidable
quadratic allocation work. Hoisting immutable ID lists is the first simple
performance candidate, in a separate change after profiling. Raw evidence search
may be deferred until a nonexact overlap needs it, but this changes diagnostic
counters; distinguish output equivalence from counter equivalence. Measure before
claiming either is the dominant runtime cause. Also profile component traversal
(`min(unseen)` per component) and retained evaluator evidence. Do not mask a graph
regression with runtime improvements.

The 48-case ladder's 1,467,144-KiB process RSS is not isolated assembler RSS.
Stream completed case evidence to disk and retain compact summaries if this is
blocking iteration; preserve coordinate gates and atomic completion status.
Reject nonfinite/nonpositive depths and preflight fixture read counts against the
assembler cap before starting a matrix. These are harness fixes, not N50 gains.

Sol must append findings under **Implementation findings**, with author, source,
hypothesis, frozen gate, exact tests/skips, per-case results, inference and next
decision. Reviewer guidance and reviewer findings remain separate. No new review
checkpoint authored by Sol, no rewriting failed historical evidence. No staging
or committing in this task, and keep local session notes out of version control.

## Historical Astra review checkpoint — 2026-09-10

This was the reviewer-authored plan at that checkpoint. The 2026-09-13 checkpoint
above now supersedes it. It replaced contradictory historical
instructions; do not restart completed phases. Implementation measurements belong
in [Sol implementation findings](sol_implementation_findings.md), with author,
hypothesis, evidence and interpretation separated. Append findings after every
attempt, including rejections and no-ops. Sol must not rewrite this guidance to
make an experiment appear to have passed. Request review when evidence warrants
a change of plan.

**Decision:** defer further phase splitting and consensus tuning. Increment A has
established containment-induced branching as the dominant clean-layout loss. Two
ownership-aware containment candidates restored clean contiguity but failed frozen
safety gates and were reverted. The evaluator now counts each uniquely projected
0-based reference coordinate once and reports disagreeing emitted bases explicitly.
The two development failures were traced to correction-bearing paths that became
damaged singletons. The fourth candidate passed its recorded development gate but
failed one strain-validation case and was reverted. This review reproduced an
unequal-length reverse-edge signature defect in that candidate. Its preservation
check detects coordinate co-occurrence, not necessarily selected traversal or
preserved spelled alleles. Therefore its claimed preservation guarantee is not
established. The latest discordance needs a paired provenance trace before it can
be called a new spelling error rather than newly resolved conflicting output. No
approach guarantees positive results before testing. The assembler remains a
bounded experimental prototype, not a validated strain-safe system.

This review inspected current source, tests, saved results and primary methods.
No separate active Sol task was found; this assesses completed repository work.
There is no clean attribution boundary in the dirty tree: do not attribute every
uncommitted change to Sol. Entries in the findings file identify their authors,
including GPT-6 continuations.

Current review evidence: `results/astra-review-20260910/audit.py`, `audit.json` and
`review.md`. Prior guide/findings were copied into that directory before editing.
The audit checks unique/equal case keys across four saved panels and reproduces
the signature defect using the frozen candidate. No actual FASTQ was run.
Sections 11–14 below specify the next three phases and conditional decisions;
they supersede historical prescriptions to protect every adjacency immediately.

Previous review evidence: `results/astra-guide-review-20260909/review.py` and `review.json`.
The directory preserves the previous guide and findings. The probe checks all
12 case identities across three recent runs, compares complete layout/consensus/
consensus-diagnostic payloads and counts containment multiplicity on existing
clean synthetic seeds. It is not a causal ablation.

## 1. Working agreement and target

- Read `AGENTS.md` every task; inspect branch and dirty files. Reviewed branch:
  `fix/support-two-quality-gates`. Preserve unrelated changes; do not stage,
  commit, reset, stash or delete experiments. Use standard branch prefixes,
  never `codex/`. Keep ignored local notes out of Git.
- User runs actual FASTQ and reference benchmarks. Agent may run unit tests and
  bounded synthetic fixtures. Give verified commands when real-input testing is due.
- One measured mechanism change at a time. No DBG revival, extra polishing loop,
  simultaneous threshold sweep or speculative replacement framework.
- Objective: higher correctly assembled span and aligned contiguity, preserving
  original molecules and plausible alternatives. N50 alone is not success.
- Graphify contains candidate/overlap/project nodes but lacks the new phase modules.
  Verify indexed relationships in current source; `.codegraph/` was absent.

## 2. Findings by implementation

| Work | Evidence and disposition | Implication |
|---|---|---|
| Equal-length exact string layout | Historical clean unique 20x seeds 93001/93002: 1198/1200 bp single contigs, versus progressive N50 92/88 | Easy-model success, not variable-fragment validation |
| Damage-compatible edges and immutable raw provenance | Retained in source, with one consensus pass | Useful foundation; edge damage model remains categorical |
| CLI IDs, conservative containment, transitive allele spelling | Repaired and tested | Preserve regressions; these are not current blockers |
| Candidate-aware evaluator | Legacy scores plus candidate transforms, diagnostic mosaics, damage-event explanations and denominators | Do not report only legacy resolved-origin metrics |
| v5 development | 36 damaged paired cases passed recorded gates | Development evidence only |
| v5 one-time validation | `95202-strains-3-0.4`: mosaic grew 135 to 254 bp; nonambiguous recovery 2189 to 2062 bp | Failed gate despite stable mosaic count |
| Phase helper repairs | Singleton exclusions, nonadjacent molecule links, branching and cycle ambiguity tested | Old helper blockers repaired; blocks are not full haplotypes |
| Component adapter and coordinate repair | Raw projection and narrow repeated-parent containment fix retained; broad alternatives rejected | Never force contradictory coordinates into a spanning-tree placement |
| Pre-layout phase audit | Development v2 damage-aware: 89 components, 87 consistent; 38 variable sites, 9 resolved links | Audit only, using filtered accepted edges and unique containments |
| Unlinked-block splitter | One strain split per seed; 35/83 duplicated overlap bases; no safety/recovery gain; 95032 N50 157 to 155 | Reverted: missing linkage is not incompatibility |
| Positive-contradiction veto | Zero splits; all 12 cases' layout/consensus/consensus-diagnostic payloads match audit v2 | Reverted no-op, not proof that whole-input linkage is absent |
| Clean layout-loss trace | 189/247 and 208/268 surviving branch edges touch ambiguous-contained nodes; truth-only removal gives N50 2358/2359 | Increment A complete; causal target established |
| Coordinate-equivalent containment | Clean contigs 153/158 to 2/1, but strain unresolved bases rise by up to 650 and mismatch-rate gate fails | Rejected and reverted; coordinate agreement is not unique evidence ownership |
| Unique post-projection owner | Clean contigs 153/158 to 2/1 and strain provenance remains stable; legacy mismatch-rate gate fails, including two +1 layout mismatches | Rejected and reverted; exposed the emitted-base denominator bias and two unresolved safety cases |
| Coordinate-deduplicated evaluator | Collapses concordant duplicate observations, abstains on conflicts and handles reverse transforms | Retained evaluator-only; use mismatch rate together with discordance and evaluated-coordinate coverage |
| Exact-edge coordinate owner | Same clean gain, but new discordance remains at 0-based reference positions 2143 and 2101 | Rejected and reverted; corrected graph edges are not the cause |
| Correction-adjacency owner | Development: 72/72 gate rows pass and clean N50 reaches 2358/2359; validation: 70/72 rows pass, one strain case fails at both stages | Rejected and reverted; signature defect reproduced, preservation claim and latest causal diagnosis need repair |

The retained baseline includes the narrow repeated-containment repair and audit.
It is not universally byte-identical to historical v5. Equal score payloads do
not establish FASTA equality. Recent evidence directories are:

- `results/variable-length-prelayout-phase-audit-dev-v2`
- `results/variable-length-phase-split-dev-v1`
- `results/variable-length-positive-phase-dev-v1`

Each contains only `summary.json`. Hashes identify source versions but cannot
recover reverted implementations. Do not claim rejected code is fully reproducible.
Historical equal-length controls and v5 development/validation evidence remain
recorded in the findings and archived previous guide; none were deleted.

The consumed 95202 diagnosis examined a final 254 bp path. Diagnostic sites 45
and 95 are 0-based contig coordinates; no retained path molecule spans both.
That does not establish absence of linkage in the whole input. The later no-op
used 95031/95032, not that failure; its counters do not distinguish every reason
for abstention. Phase work is deferred for lack of demonstrated benefit, not
because impossibility was proved.

## 3. Where progress is blocked

Even the clean variable-length controls are highly fragmented:

| Seed | Nodes | Unique containments retired | Multiple containment placements | Contigs | N50 | Longest | Pre-layout components | Ambiguous ends |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 95031 | 208 | 35 | 105 | 153 | 150 | 200 | 6 | 127 |
| 95032 | 208 | 27 | 108 | 158 | 160 | 250 | 4 | 133 |

These are clean 31–200 bp fragments from existing 2400 bp unique-reference
fixtures. A child retires only if one exact parent/orientation/offset exists
globally. Multiple parent reads can represent the same locus; multiple parents
need not mean multiple genomic placements. But arbitrary parent selection can
lose repeat/strain evidence. These counts justify a containment/branching audit,
not deleting all contained reads. Component counts are not attainable N50 bounds.

Other losses remain unranked: missed anchors, best-placement ties, occurrence caps,
accepted-edge topology, two-hop reduction and reciprocal degree-one selection.
Do not assume candidate recall is dominant either. An audit that leaves layout
unchanged cannot itself repair these losses.

## 4. Current source boundaries

- `damage_string_graph.assemble` canonicalizes identical sequences up to reverse
  complement while retaining raw groups. Graph nodes are not molecules. It
  enumerates exact containments, accepts edges, audits coordinates, retires
  eligible children, projects paths and expands raw provenance.
- `_candidate_alignments` in `overlap_assembly.py` selects best anchor-supported
  placements and suppresses ties. Reuse `stage_trace`. The existing
  `extension_candidate_recall.py` is a small pair/decoy test, not a full graph audit.
- `_project_master_edges` in `overlap_graph.py` uses two-hop reduction and
  reciprocal nonbranching paths, not full Myers reduction. The string-layout
  experiment wraps the same source; it is not an independent oracle.
- The damage-edge rule uses independent, unambiguous raw molecule support and
  directional terminal evidence. Positive rates enable categorical events;
  halving positive rates does not test quantitative edge calibration.
- `phase_graph.py` uses accepted edges and unique containment constraints, not all
  input alternatives. Offset/orientation conflicts abstain. `phase_blocks.py`
  retains supported pairs, molecule IDs, ambiguous links and excluded evidence.
  Two-molecule support is a heuristic, not calibrated error rejection.
- Dense coordinate columns, all-pairs marker enumeration and unconditional
  pre-layout auditing need profiling before scaling. Exact containment search
  also compares sequence pairs. Do not call these production-ready.
- Raw placement creation is in source, not only experiment code.
  `project_raw_consensus` returns reads and diagnostics, not a pool.
- `damage-assemble` has an overridable default cap of 1000 reads. Raising the cap
  is not evidence of 100k readiness.
- Graph outputs retain singletons; progressive centers have different membership.
  Match membership/filtering or stratify comparisons. Some controlled runs compute
  progressive diagnostics before replacing layout: those are not graph counters.

## 5. Increment A — completed: locate one recoverable layout loss

Completed in `results/layout-loss-trace-dev-v3`. The bounded trace and the
four-read forward/reverse reproducer established the following:

- Every selected clean edge was truth-compatible. Candidate discovery missed
  40/1072 and 26/1104 proper dovetails, while 247 and 268 survived reduction as
  branches. Missing candidates were not the dominant measured category.
- Of 140/135 contained nodes, 105/108 had multiple parent representations of one
  truth locus. No clean containment candidate mapped to multiple truth loci.
- Ambiguous-contained nodes touched 189/247 and 208/268 surviving branch edges.
- A truth-only ablation removing all contained nodes produced 2/1 contigs with
  N50 2358/2359. It diagnosed causality but is forbidden as assembly policy.
- The smallest four-read fixture shows one multiply contained read blocking a
  correct 140 bp reciprocal path in both orientations, producing 60, 80 and
  120 bp outputs.

The original procedure is retained below as the audit contract for reproducing
or extending this diagnosis. Do not rerun it as if increment A were unfinished.

1. Freeze current source/settings and seed roles in a new result directory.
   Use clean 95031/95032 first. Truth stays private to evaluation; the assembler
   receives only reads, qualities and molecule IDs.
2. Reuse graph stages and tracing. Enumerate exact overlaps exhaustively only on
   these small fixtures as an independent diagnostic. Separate proper dovetails
   from containment and record oriented offsets. Removing a redundant edge with
   an equivalent path is not a lost assembly opportunity.
3. Follow expected adjacencies through anchors, ties/caps, identity, acceptance,
   containment, reduction and reciprocal selection. Include missed candidates
   and selected extra edges. Count physical edges once, checking reverse mirrors
   instead of treating them as independent replicates.
4. Distinguish multiple parent representations of one coordinate-compatible locus
   from repeated offsets, incompatible parent paths and strain alternatives.
   Capture parent edges and surviving degrees. Truth may diagnose this distinction;
   an eventual decision must use read evidence.
5. Locate actual broken output junctions and lost correctly assembled span.
   Choose a recoverable loss, not the category with the most redundant edges.
   Reproduce it in a tiny clean forward/reverse fixture.

Check unique/equal case keys, complete terminal dispositions, coordinates and raw
membership. Include repeated-placement and uneven-strain negatives before accepting
any simplification. If containment multiplicity is not causal, report that and
follow the demonstrated loss. Do not open fresh validation seeds in this increment.

## 6. Increment B — evaluator repair complete; next causal trace

The first candidate is preserved in
`results/variable-length-coordinate-containment-dev-v1` and was rejected. It
retired a multiply contained node when candidate parents implied one coordinate,
then projected the child through compatible parent layouts. Clean N50 reached
2358/2359, but strain cases gained cross-reference unresolved output (up to
650 bases) and 16/24 mode/case comparisons violated at least one predeclared gate.
Coordinate agreement does not prove that shared evidence has one strain owner.

The second candidate, preserved in
`results/variable-length-unique-owner-containment-dev-v1`, attached a contained
molecule only after one post-projection output coordinate remained. It removed the
strain-provenance regression and kept clean N50 2358/2359, but failed the literal
legacy mismatch-rate gate. A retained evaluator repair then showed that most rate
failures were denominator artifacts: concordant duplicate output had previously
been counted repeatedly. It also exposed one new discordant coordinate in each of
two damage-aware inputs.

The third candidate, preserved in
`results/variable-length-exact-coordinate-containment-dev-v1`, allowed only exact
accepted edges to prove parent-coordinate equivalence. It retained the clean gain
but reproduced both discordances at 0-based reference positions 2143 and 2101.
This falsifies the corrected-edge hypothesis. It was reverted.

Do not restore either implementation unchanged. The completed trace shows that
suppression changes path spelling. Baseline raw reads 21 and 16 belong to corrected
two-read paths of 119 and 111 bp with zero mismatches. The candidate emits them as
97 and 109 bp singletons, losing partners 134 and 146 and retaining a terminal
damage observation at contig positions 3 and 106. Raw evidence is high-quality T at
5-prime distance 3 or 2; reverse orientation produces the discordant A.

The fourth candidate computed the unsuppressed projection first and identified
selected correction-bearing adjacencies. Apply containment suppression provisionally,
then required each such adjacency to remain in one output at the same orientation
and relative 0-based offset. It restored proposed suppressions by component and
recomputed monotonically. Development passed all 72 normal/error-only comparisons.
Frozen validation seeds 96011/96012 were then consumed once. Seed 96012's
damage-aware 80:20 strain case gained a discordant site at reference 1 position
1539 (0-based), despite N50 150→155 and nonambiguous recovery 4178→4341. The
candidate was reverted and these seeds must not be reused as confirmation.

Do not restore this candidate unchanged or immediately implement an all-adjacency
guard. First complete Phase 1 in section 11: paired diagnosis of position 1539,
correct physical-edge identity and a gate that distinguishes common-coordinate
regression from newly evaluated output. Inspect consumed seed 96012 diagnostically;
it can become a regression fixture but can never again confirm generalization.
Develop the next mechanism on separate minimal fixtures and established development
inputs. Keep the historical rejection; do not retroactively change its gate.

Write the failing reproducer first. State the decision rule and observable evidence.
A may identify redundant containment, reduction or candidate loss; none is chosen
in advance as a guaranteed solution. Preserve these invariants:

- All coordinates are 0-based, intervals half-open. Raw position p maps to offset+p
  forward or offset+L-1-p reverse. Reversing a C-base contig maps full-read offset o
  to C-o-L and flips orientation. Preserve clipped bounds in raw orientation.
- Raw bases/qualities and raw-end distances never become corrected-node evidence.
  Count a molecule once per column; conflicting/ambiguous placements stay explicit.
- Reverse edges agree; different offsets do not silently overwrite. Identical
  spelled paths retain compatible provenance; cycles are not silently linearized.
- Retiring a node requires retained connectivity and evidence ownership. Do not
  choose the longest parent, discard all children, or use reference truth to invent
  uniqueness. Preserve real repeat/strain alternatives.
- Transitive replacements agree in coordinates and spelled alleles. Do not loosen
  damage thresholds or add polishing to conceal a layout failure.

Run focused regressions then the full suite. Freeze and evaluate the existing
12-case variable-length development panel with exact and damage-aware controls.
Retain the already implemented separate error-only controls (`--error-only`):
the main scenarios mix error with damage and assign errors Q10. Do not reimplement
this panel or count it as a new contribution.
Profile magnitude sensitivity belongs to consensus; zero versus positive profiles
can test categorical graph behavior.

Historical gate (retain when reproducing old runs): fix the tiny break and improve correct span
on targeted clean cases; no per-case increase in diagnostic-mosaic/incompatible
counts or bases, nonduplicated mismatch rate, discordant reference coordinates or
unresolved bases; no loss of evaluated reference coordinates or sequence-resolved
nonambiguous recovery. Continue reporting the legacy emitted-base rate so the
denominator effect remains visible. Report N50, longest, truth-compatible contiguity,
duplicate observations, time and RSS. A clean topology fix must work in exact mode
too. Aggregate gains cannot hide failed cases. Obtain review before changing a
tolerance after seeing results.

The existing comparison script automates only part of this contract; section 11
requires explicit completeness checks and paired-coordinate diagnostics before a
new gate is frozen. Any revised gate is prospective and recorded separately.

After one failed candidate, preserve source/results and diagnose the failure. Do
not cycle through phase splitters, thresholds or new seeds without a new causal
observation. Stop rules prevent repetition; they do not mean hiding negative results.

## 7. Next increment C — frozen validation and user-run prototype

After B passes, predeclare untouched validation seeds and the same gates, freeze
source and evaluate once. Seeds 93001/93002, 95031/95032 and inspected 95101–95104
are development; 94011/94012 are historical validation; 95201–95204 are consumed
validation; 96011/96012 are now consumed validation too. Never reuse them as
independent confirmation or silently relabel runs. The current driver's `--held-out`
default points to consumed seeds, and explicit `--seeds` labels even validation
as development. Repair metadata and register roles before another validation run.

Profile representative bounded synthetic sizes before increasing real-input caps.
Measure containment, indexing, reduction, phase audit and consensus separately.
Optimize measured costs with sequence/provenance equivalence checks. Wrapper/source
equivalence is a regression check, not independent correctness validation.

A user-run bounded actual FASTQ smoke test may check CLI operation, deterministic
output, provenance joins and resources once those checks pass. Label it operational,
not evidence of better N50 or strain safety. Do not require another phasing invention
before offering that test. A matched 100k comparison additionally needs resource
evidence and the frozen synthetic gate. The agent never runs actual FASTQ itself.

Before giving commands, inspect `damage-assemble --help`, verify current paths and
profiles, and use a new output directory with an explicit cap. Compare layout at
equal consensus status and output-length filtering, then one-pass consensus effects.
Remap changed assemblies; byte-identical FASTAs do not need another MetaQUAST run.
Real-reference agreement remains a proxy, not read-origin truth.

## 8. Literature and repository cross-check

- [Myers 2005](https://doi.org/10.1093/bioinformatics/bti1114) and
  [Simpson and Durbin 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2881401/):
  separate redundant containment/transitive overlaps from unitig extension. This
  motivates a graph audit, not assuming the current reducer is complete.
- [Kamath et al. 2024](https://doi.org/10.1101/gr.279311.124) and
  [RAFT](https://github.com/at-cg/RAFT): contained-read deletion can remove needed
  haplotype coverage. Their long-read setting differs; borrow the preservation
  concern, not thresholds or a blanket keep/delete policy.
- [CarpeDeam paper](https://doi.org/10.1186/s13059-025-03839-5) and
  [official repository](https://github.com/LouisPwr/CarpeDeam): separates read
  extension from later contig merging and documents a contiguity/misassembly
  tradeoff. Staged extension is relevant after basic layout loss is measured;
  lower identity or unsafe merging does not establish this project's safety goal.
- [POLYTE paper](https://doi.org/10.1093/bioinformatics/btz255) and
  [HaploConduct](https://github.com/HaploConduct/HaploConduct): overlap-based
  haplotype preservation is relevant, but known-ploidy/paired-end assumptions
  differ from merged ancient metagenomic fragments. Never invent read pairing
  or molecule IDs to supply missing linkage.
- [Li 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4937194/) separates layout
  and correction; its long-read simplifications need new justification here.
  [Li 2011](https://doi.org/10.1093/bioinformatics/btr509) motivates observation
  likelihoods, not calibrated claims about join safety from nominal posteriors.

Primary search records and repository documentation were checked. Some full-text
publisher/PMC opens were blocked. This is not a source-line audit of every external
implementation. Before borrowing code, inspect its version, implementation and
license; distinguish the borrowed principle from the Anvaya-specific heuristic.

## 9. Commands, reproducibility and reporting

Unit command from the activated project environment; verify dependencies first:

```bash
cd /mnt/e/Projects/anvaya
python -c "import sys, pysam; print(sys.executable); print(pysam.__version__)"
PYTHONPATH=src python -m unittest discover -s tests -q -b
git diff --check
```

The former `/tmp/anvaya-audit-tests/bin/python` is absent in this review environment.
On 2026-09-10, `/usr/bin/python3` ran 300 tests: 278 passed, 21 skipped and one
test-module import error (`pysam` missing). This is incomplete verification, not
evidence of a new assembler regression. Earlier suite counts refer to different
source/environment states and were not reproduced here. Never equate tests run
with tests passed; recheck external tools before claiming integration coverage.

Development matrix command after snapshotting and declaring the gate; use a new,
unused result name:

```bash
PYTHONPATH=src python experiments/variable_length_validation.py \
  --seeds 95031 95032 --output results/variable-length-layout-candidate-v1
```

The driver now hashes itself, the scorer and all `src/anvaya/*.py`, not just four
files. Hashes do not replace recoverable source. Preserve source snapshots, exact
commands/settings, environment versions, time/RSS, stdout/stderr, run status, seed
roles, FASTAs and raw placements before reverting. Partial experiments count too.
Check case-key joins and denominators; keep nonambiguous/repeat-masked recovery
separate from whole-reference coverage. Reference length is not sample depth.

Historical paths, to verify before issuing actual-data commands:

```text
data/carpedeam_p0/simulated/synth-EMN001/upload/synth-EMN001/p0-subsets/EMN001-simulated-100k-min31.fq.gz
data/carpedeam_p0/simulated/synth-EMN001/upload/synth-EMN001/p0-subsets/EMN001-simulated-500k-min31.fq.gz
data/carpedeam_p0/simulated/synth-EMN001/upload/synth-EMN001/refs/all.fa
results/raw-consensus-100k/profile5p.prof
results/raw-consensus-100k/profile3p.prof
```

The profile is supplied, not established as fitted to EMN001. Historical MetaQUAST
settings: `--min-contig 31 --min-identity 90 --threads 2`; verify executable paths.
The old overlap proxy resolved only 1264/97964 reads: its 92.6% recall does not
represent the whole dataset. The 474,607,190-base reference sum is not coverage.

Each findings entry: author/date, hypothesis, exact files/source snapshot, tests
run/passed/skipped, results per seed, failed gates, measurements versus inference,
and one next action. No-ops need attempted/eligible/acted/abstained counts by reason.
Zero actions cannot by itself locate missing information. Never claim byte equality
without comparing bytes. Implementation reports stay separate from reviewer guidance.

## 10. Reconciliation of every previous guide section

| Previous section | Disposition |
|---|---|
| Checkpoint and blockers | Fixed bugs marked historical; impossibility claim removed; one active task |
| Increment 1 repairs | Invariants retained in sections 2 and 6 |
| Evaluation conclusions | Updated scorer/manifest acknowledged; denominator and snapshot requirements retained |
| Retained/rejected work | Consolidated in section 2, historical findings and archived guide |
| Next three increments | Contradictory completed/rejected instructions replaced with A/B/C |
| Objective and working agreement | Overlap goal, surgical scope, user-run actual-data rule retained |
| Established baseline | Easy equal-length success separated from variable-length failure |
| Implementation boundaries | Updated adapter ownership, evaluator, overridable cap and test count |
| Literature | Added containment counterevidence, verified repositories and transfer limits |
| Original Phase 1 | Completed repairs distinguished from remaining candidate/graph risks |
| Original Phase 2 | Infrastructure retained; failed generalization prevents validated status |
| Original Phase 3 | Variable-length CLI exists; validation and scaling remain open |
| Commands and output discipline | Bounded operational smoke distinguished from comparative benchmark |
| Completion report | Separate author-labeled findings, failed gates and abstention reasons required |
| Suggested next instruction | A is complete; follow section 11 before another mechanism change |

The full previous guide remains in local review evidence. Do not resurrect its
superseded instructions as a parallel to-do list.

## 11. Historical Phase 1 — make the next decision trustworthy

**Deliverable:** one paired failure diagnosis, tested preservation semantics and
an executable prospective gate. This is a bounded prerequisite, not another broad
benchmarking or polishing project. No new validation seeds and no default assembler
change in this phase. Keep each of the following as a separate reviewable change.

### 11.1 Reconstruct the failed case by raw provenance

Use the frozen baseline/candidate summaries and source under
`results/variable-length-correction-adjacency-containment-validation-v1` and the
development snapshot. Verify imported module paths and hashes before reproducing
the small synthetic case; do not accidentally compare the current source to itself.

For seed 96012, `damage_error_strains_80_20`, reference index 1 position 1539:

1. Export baseline/candidate FASTA, selected traversal steps, raw molecule groups,
   placements and evaluator candidate transforms. Join by immutable raw identities,
   orientation and coordinates. `unitig_160` is not a stable cross-run identifier.
2. Identify the A-bearing candidate output and locate its supporting raw observations
   in baseline outputs, including previously unresolved outputs. Record raw base,
   quality, raw 5-prime/3-prime distance, multiplicity and eligible/ambiguous status.
3. Decide separately whether spelling changed, projection became resolvable, a new
   path appeared, or alternative ownership changed. Show the baseline and candidate
   local spelled sequences and the exact decision responsible. The coordinate was
   already covered in baseline; calling it simply a newly recovered locus is wrong.
4. Reduce the causal mechanism to a tiny fixture independent of the validation seed.
   If summaries cannot recover the lineage, regenerate this synthetic case from the
   snapshot and save it. Do not guess from aggregate metrics or contig ordinals.

All reference, contig and raw positions are **0-based**; intervals are half-open.
For contig length C, a forward interval [a,b) becomes [C-b,C-a) after reversal.
Reference truth is evaluator-only and must never choose assembly ownership.

**Exit:** an evidence table links one changed output to its raw evidence and causal
operation, or explicitly documents why origin ambiguity prevents a conclusion.
An unresolved trace is not permission to implement a stronger guard by intuition.

### 11.2 Test physical identity and actual spelling preservation

The archived `_edge_signature` must not be copied unchanged. For source length Ls,
target length Lt, overlap O and forward shift s = Ls-O, the mirror shift is Lt-O,
equivalently Lt-Ls+s. For a source-oriented correction position p, the mirrored
position is Lt-1-(p-s), and its base is complemented. Use the established reverse
edge constructor where possible instead of maintaining a second convention.

Specify two distinct checks: physical adjacency identity (orientation/relative
placement) and spelling/evidence preservation. Endpoint co-location alone proves
neither that a traversal was selected nor that its correction reached the output.
If selected traversal information is currently discarded, expose only the minimal
trace needed at path construction; do not build another graph framework.

Required synthetic regressions before testing a new candidate:

- unequal lengths, both orientations: an edge and its constructed mirror identify
  one physical adjacency; distinct offsets stay distinct;
- the existing four-read clean containment break improves in forward and reverse;
- a correction-bearing path becomes an uncorrected singleton after suppression:
  the candidate must preserve the necessary spelling/evidence;
- coalesced layouts with co-located endpoints from different paths cannot falsely
  satisfy a selected-path check;
- same coordinates but conflicting spelled alleles do not count as preserved;
- ambiguous repeat/strain ownership remains explicit; duplicate raw observations
  cannot become independent molecule support;
- rollback stops monotonically, verifies its final invariant and fails closed to
  the baseline projection if the invariant cannot be established;
- a failure in one disconnected component does not cancel a safe change elsewhere.

Do not encode the historical reference base or a seed-specific read number in a
production rule. Keep tiny fixtures synthetic; this phase requires no actual FASTQ.

### 11.3 Make the comparison contract explicit

Extend the existing evaluator/driver only enough to emit a paired table over the
union of `(reference_id, position_0based)` coordinates. States are: concordant
correct, concordant incorrect, discordant, and not uniquely evaluated. The last
state needs an explicit unresolved/uncovered reason. Do not collapse ambiguity
into correctness or treat missing data as zero errors.

Report transitions, not just totals: correct→incorrect, correct→discordant,
correct→unresolved/uncovered, incorrect→correct, newly evaluated correct/incorrect/
discordant, plus unchanged states. Keep duplicate observations and unplaced outputs
separate. A sequence-resolved transform is conditional evidence; retain uniquely
origin-resolved and sequence-resolved strata so sequence-based resolution does not
silently improve the apparent accuracy denominator.

Prospective acceptance contract, frozen before the next candidate matrix:

- Zero new errors/conflicts/losses on previously correct uniquely evaluated
  coordinates; no new demonstrated incompatible joins or diagnostic mosaics.
- Report new-region incorrect/discordant output separately; it never counts as
  correct recovery. A novel conflict needs a causal classification, not automatic
  approval based on aggregate gain. If a tolerance is proposed, obtain reviewer
  agreement before running; default to requiring no unexplained new conflicts.
- Retain old total mismatch/discordance/coverage and legacy-rate diagnostics.
  Report numerator and denominator; zero denominator means not evaluable.
- Measure unique correct recovered bases and truth-compatible contiguity separately.
  Coverage with mismatches is not correct span. Supply known repeat masks in the
  synthetic evaluator or explicitly label results as unmasked coverage.
- Positive benefit is mandatory: fix the target fixture and improve correct span
  or truth-compatible contiguity on targeted development inputs. At least one
  damaged development case must benefit before this is called an aDNA improvement.
  A clean-only gain may be reported as a clean-layout improvement, with limited scope.
- All per-case safety checks must pass. Report independent input-case counts and
  correlated mode/stage rows separately; no significance claim from row totals.

Before comparing, assert unique and equal case keys, expected modes/stages, complete
statuses, matched inputs/settings and evaluated-coordinate denominators. Add tiny
tests for duplicate/missing cases, missing metrics, zero denominator, discordance
hiding an error and all-no-op outputs. A failed gate must return nonzero status.

Repair seed metadata with an explicit role and a consumed-seed manifest. Reject
contradictory flags and an attempt to claim known consumed seeds as fresh validation.
Record historical metadata inconsistencies in new audit records; do not rewrite
old evidence. Source hashes must cover actually imported modules; retain recoverable
source, inputs/settings and environment, not hashes alone.

**Phase 1 stop:** if the same result has merely been reformatted, stop diagnostics.
Write the causal conclusion and select exactly one Phase 2 branch. Do not start
another matrix until that conclusion or a precise unresolved evidence gap exists.

## 12. Historical Phase 2 — one containment change, chosen by evidence

| Phase 1 result | Next implementation | Required evidence to proceed |
|---|---|---|
| Lost selected traversal changes a previously corrected base | Preserve the necessary traversal/spelling or restore the particular responsible suppression | Tiny damaged fixture repaired; common-coordinate regressions absent |
| Spelling unchanged but an existing conflicting output becomes uniquely projected | Repair ownership/projection only if raw evidence shows it is wrong; otherwise classify the evaluation transition | Provenance join identifies old output; no invented unique origin |
| Genuine new alternative path loses allele evidence | Retain the necessary child/path using read-supported constraints | Positive contradiction or indispensable path evidence; missing linkage alone is insufficient |
| Mirror/signature bug causes spurious restoration | Correct canonical identity and remeasure action counts before expanding the guard | Reverse-invariance tests and case-level rollback delta |
| Guard passes by undoing all damaged changes | Localize restoration to demonstrated dependencies; do not relax safety thresholds | Safe independent component still improves and targeted damaged case benefits |
| Output unchanged and no new causal evidence | Stop this mechanism and report the no-op | One bounded next audit at the first measured lost stage, not a fresh threshold sweep |
| Containment no longer dominates remaining losses | Rank candidates/reduction/branch selection with existing stage trace | Lost correct junctions, not simply edge counts, justify the next target |
| Better N50 but incorrect/unresolved output increases | Reject or request review of a predeclared tradeoff | Never accept on N50 alone or silently alter the gate |

Try one minimal mechanism first. Do not preserve every baseline edge automatically:
that can freeze the fragmentation being repaired. Do not import the whole rejected
implementation and change several rules simultaneously. No new phase splitter,
consensus loop, relaxed identity threshold or DBG code belongs in this increment.

Run focused tests, then the dependency-complete unit suite. Freeze the candidate
and baseline, then run the existing 95031/95032 main and separate error-only panels
with exact/damage-aware modes and layout/one-pass consensus stages. Known seeds are
development/regression only. Example commands, after activating a verified Python
environment and choosing unused output directories:

```bash
PYTHONPATH=src python experiments/variable_length_validation.py \
  --seeds 95031 95032 --output results/next-containment-dev-v1
PYTHONPATH=src python experiments/variable_length_validation.py \
  --seeds 95031 95032 --error-only --output results/next-containment-error-dev-v1
```

These are current driver commands, not a complete frozen baseline/candidate runner.
Invoke each source snapshot in isolation and record its imported path. Section 11's
metadata repair must precede any new validation commands.

Measure proposed/eligible/suppressed/restored containments, protected physical edges,
actual selected traversals, rollback reasons/iterations, output changes, correct
span, contiguity, errors and resources. A rate of zero actions is a no-op, even when
every safety row passes. Match exact-mode topology checks to damage-aware ones.

After a failure, preserve source and evidence, classify it once, and request review
if a different invariant is needed. Do not consume another seed panel to search for
a pass. After two mechanism attempts without a demonstrated damaged-data benefit,
stop this line and provide a ranked alternative backed by existing stage-loss data.
That is a review checkpoint, not authorization for a wholesale assembler rewrite.

## 13. Historical Phase 3 — validation and a usable prototype

There are two distinct readiness tracks; report each explicitly.

**Biological improvement:** after Phase 2 passes a predeclared gate, freeze source,
settings and validation manifest, reserve genuinely unused seeds, then evaluate
once. Include clean, damage-only, error-only, mixed error/damage, repeats, uneven
strains and profile mismatch. Existing fixtures have independent substitution errors
at Q10 and omit indels/PCR duplication/real abundance distributions; passing them
does not establish broad metagenome robustness. Add such realism only as a separate
declared panel after the current mechanism works, not as an endless prerequisite.
Failure remains recorded; inspected validation seeds become diagnostic, never fresh.

**Operational use:** can advance without claiming biological superiority. Confirm
CLI/provenance/determinism tests, then profile bounded synthetic sizes up to the
intended smoke cap. Record reads/bases, wall time and peak RSS with containment,
indexing, reduction, phase audit and consensus timings. Set the resource budget
before running and stop when exceeded. Profile first; pairwise containment and
dense audits are known risks, not yet measured dominant runtime at 100k.

Once those checks pass, give the user a verified bounded FASTQ smoke command with
an explicit cap, input/profile paths and a new output directory. The default cap
rejects oversized input; it does not select the first N reads. If a small subset is
needed, provide a streaming FASTQ extraction command and validate its record count.
The user runs it. Report this as CLI/resource validation, not improved N50.

If the smoke succeeds, proceed to a matched larger comparison only within measured
resources and with the frozen synthetic acceptance result for the candidate.
Do not simply raise the 1000-read default to 100k. Optimize only measured hot paths,
preserve sequence/provenance equivalence on fixtures, then repeat the relevant
resource check. No need to rerun every biological panel for a byte-equivalent change.

For actual comparisons, match input molecules, output filtering and consensus stage.
Report N50/NA50, aligned span, unaligned output, mismatch/indel burden, duplication,
misassemblies and resources against the retained baseline and available CarpeDeam
result. Reference alignment is a proxy; zero reported misassemblies in very short
contigs does not demonstrate safe long joins. If only consensus accuracy improves,
do not call it a contiguity gain. If FASTA is identical, record the no-op and do not
ask the user to repeat MetaQUAST for the same sequences.

**Current distance to goal:** the overlap prototype and CLI exist, but three gates
remain open: a retained damaged-data contiguity gain, independent validation, and
resource-backed actual-input evaluation. These are milestones, not three guaranteed
patches or a reliable percentage/ETA. A bounded smoke test is much closer than a
validated 100k ancient-metagenome assembler.

**Phase 3 status, 2026-09-12:** the selected-edge containment candidate failed the
fresh independent zero-regression gate in uneven-strain cases and has been reverted
without tuning. The restored baseline passed the predeclared 1,000-read synthetic
resource and determinism checks. A bounded 1,000-read EMN001 smoke also passed its
operational/resource gate but produced 994 contigs and zero supported consensus
changes, so it is too sparse to evaluate biological benefit. Consumed-seed diagnosis
shows that pre-reduction containment changes graph degrees and enables an unsafe
exact cross-strain path elsewhere in the component. Close this containment branch.
Indexed exact-containment discovery now matches the exhaustive oracle and frozen
outputs, reduces the measured 1,000-read containment stage 24-fold, and passes the
predeclared 10,000-read synthetic resource gate. The exactly 10,000-record EMN001
smoke also passed at 5.57 seconds and 83,924 KiB. It emitted 9,343 contigs with
N50 61 bp from reads with N50 60 bp; 351 paths joined 727 reads. Twenty-four
damage-compatible candidate overlaps were accepted and seven overlap bases were
corrected, but consensus made zero changes because nearly all columns lacked
three independent molecules and all other covered columns already agreed.

This result changes the immediate decision. Candidate discovery is the largest
measured stage (2.241 seconds), but it is not yet an operational blocker and a
runtime-only rewrite cannot improve N50. Do not tune consensus or weaken overlap
identity. A provisional 25,000-record EMN001 run passed resources at 16.56 seconds
and 195,688 KiB but emitted 22,702 contigs with N50 62 bp from reads with N50
60 bp. Both response gates fail. Its diagnostics predate or bypass the current
candidate-stage instrumentation, so rerun that exact frozen FASTQ with
`PYTHONPATH=src` and record the imported module path and SHA256. That rerun is now
complete: it reproduced 22,702 contigs in 19.26 seconds and 195,328 KiB, with both
candidate diagnostic objects present. Of 1,420 input physical graph edges, 1,302
(91.7%) survived transitive and branch filtering, so graph reduction is not the
dominant fragmentation stage. Zero anchors were occurrence-capped. The aggregate
candidate funnel cannot distinguish absent true overlaps from correct filtering.
The current-source reference-proxy opportunity audit now reproduces 565 selected
of 610 expected directed dovetails (92.62%): 23 had no anchor vote and 22 failed
R/Y identity. No DNA-identity or occurrence-cap loss was observed. Only 1,264 of
97,964 mapped query groups were eligible, so this is evidence from an easy,
selection-biased subset rather than read-origin truth. It is sufficient to reject
a broad anchor or R/Y relaxation, which would spend repeat/strain safety for at
most 7.4% proxy recall. The bounded 50k coverage gate is now authorized with no
algorithm or threshold changes.

The predeclared 50k resource stop is 120 seconds or 1,048,576 KiB RSS. Response
requires output/input at most 90% and assembly N50 at least 1.10 times input-read
N50. If resources fail, stop and profile candidate runtime. If either response
criterion fails, do not run 100k; quantify reference-supported physical overlap
opportunity on the frozen 50k subset. If both pass, compare layout and consensus
against the reference before changing assembly policy. The manifest is
`results/damage-string-graph-emn001-smoke-50000-manifest-v1.json`.

The 50k run is complete. Resources pass at 62.65 seconds and 326,504 KiB, and
output/input passes at 43,478/50,000 (86.96%). N50 changes only 62 to 65 bp
(1.048), failing the 1.10 response gate. Stop before 100k. Run the exact-cohort
opportunity audit using `--read-names-fastq` with the frozen 50k FASTQ and existing
query-name BAM. Require requested/matched equality. Use per-reference eligible
depth, overlap-component lengths and expected-dovetail stages to decide whether
the available cohort lacks physical connectivity or the candidate finder loses
reference-supported joins. Do not change consensus, graph reduction or thresholds
until that result is recorded.

The exact-cohort audit matched all 50,000 names. Only 578 reads were eligible and
they formed 535 optimistic >=30-bp overlap components; the longest was 178 bp.
Candidate recall was 65/70 (92.86%), with five no-anchor misses and no identity
loss. This closes EMN001-prefix tuning: measurable physical connectivity is too
sparse to diagnose long-contig behavior, and 100k remains blocked.

A truth-known 2x/5x/10x/20x ladder exposed fragmentation. The following causal
claim and proposed heuristic are historical and superseded by the 2026-09-13
review; do not execute them. On clean
unique references, mean N50 falls 237→180→149→133.5 bp while reads rise
42→104→208→416 and outputs rise 16.5→53→149.5→343. At clean 20x, an average
3,024.5 physical edges yields only 20 reciprocal edges despite zero sequence
mismatches. The next candidate must act only after transitive reduction: choose a
unique best-overlap outgoing and incoming edge, retain their reciprocal
intersection, and abstain on ties. Do not change anchors, identity, damage scoring,
containment, consensus or transitive reduction in the same candidate.

Development acceptance: clean 20x N50 must increase by at least 50% and retain
zero discordant/mosaic/incompatible bases. Across all repeat and 80:20 strain
development rows at 2x/5x/10x/20x, discordant, mosaic, incompatible and unresolved
bases must not increase versus branchless baseline. A failure is preserved and
reverted without threshold tuning; a pass receives fresh independent seeds.
Keep ambiguous contained reads through reduction; path-policy work remains gated
against strain mosaics. Detailed measurements are in the Phase 3 and indexed-
containment findings.

## 14. Required Sol findings format and reviewer boundary

Append to `docs/sol_implementation_findings.md` after each attempt. Never insert
implementation claims into reviewer instructions or silently rewrite a failed gate.
Use this template, including rejected and incomplete runs:

```text
Author/model, date, branch, source snapshot + hashes, imported module paths
Hypothesis and single changed mechanism:
Predeclared success/failure criteria and seed roles:
Tests: run / passed / skipped / failed; missing prerequisites:
Inputs: case IDs, modes/stages, key-join validation, environment:
Results: per-case before/after correctness, contiguity, coverage and resources:
Actions: eligible / changed / restored / abstained, with reasons:
Finding (measured):
Interpretation (inference, with uncertainty):
Disposition: retained / rejected / no-op / incomplete:
Next action: exactly one branch from section 12, or reviewer question:
User-run test needed: yes/no; exact verified command if yes:
```

A later reviewer should be able to identify what changed, why, which criterion
failed and the next evidence needed without rereading the entire conversation.
Keep this guide reviewer-authored; additions from Sol belong in findings first.
