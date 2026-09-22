# Sol implementation findings

> **Active reviewer contract, Astra, 2026-09-20:** only
> [the current guide](sol_implementation_guide.md) authorizes next steps. All older
> next-action instructions below are historical. Implementer findings append under
> author/model/date headings; do not edit reviewer criteria after observing results.

## Reviewer findings — Astra, 2026-09-20 (not an implementation result)

- **Correction:** failed greedy/exact-placement opportunity tests do not establish
  a biological contiguity ceiling. The 44-minute paired audit remains a null under
  its restrictive model; zero reported conflicts means the reproduced global
  conflict bug does not explain that real result.
- **New synthetic reproduction:** all 138 clean merges equal original fragments,
  yet exact assembly N50 drops 2,399→99; 265 ambiguous containments and 64 reduced
  full-containment directed edges remain. Verify variable-length graph semantics
  before rejecting mate merging. This consumed fixture is NOT validation.
- **Evidence:** `experiments/review_20260920.py`, source-hashed outputs and unit log
  in `results/reviewer-20260920/`. Python 3.13.15: 348 tests RUN, 347 PASSED,
  one SKIPPED. No actual FASTQ assembly rerun or production algorithm changes.
- **Decision:** authorize a conventional exact unitig core separated from immutable
  raw evidence, then one damage integration and frozen validation. Previous
  containment patches/linkers and paired sweeps stay closed. See the
  [review](astra_review_20260920.md) and [active guide](sol_implementation_guide.md).
- **Next:** Phase A; no user actual-file test yet. Append implementer results under
  distinct model/author/date headings; never rewrite reviewer criteria.

## Reviewer findings — Astra, 2026-09-15 (historical; superseded above)

- **Assessment:** useful runtime and diagnostic work, but no held-out-safe damaged
  contiguity improvement. GPT-6 continuations and the accumulating guide share
  responsibility for the repeated containment approach. Do not blame every dirty
  file or failed experiment on Sol.
- **Verified evidence:** current graph hash matches restored Phase A; frozen main
  validation fails four rows and error-only passes. Small executable probes show
  the paired gate ignores input hashes/joint omissions and lacks clean-gain
  enforcement. See `results/astra-review-20260915/` and the
  [detailed review](astra_review_20260915.md).
- **Interpretation corrections:** 48/24 are correlated case combinations, not
  independent replicates. The literal four-read reduction is RNG-free at runtime
  but derived from consumed validation; it is not independent validation or a
  proof of a globally smallest counterexample. Historical entries remain intact.
- **Tests this review:** 326 run, 325 passed, one skipped; Python 3.12 with
  pysam 0.24.0. Declared Python >=3.13 remains unverified. Do not count skips as
  passes. No actual FASTQ or rejected assembly matrix was rerun.
- **Decision:** close ambiguous-containment suppression. Fix the experiment
  contract, then one bounded exact raw-supported post-layout merging attempt,
  reusing the existing linker and composing raw provenance. No useful joins or
  a safety failure closes that route and triggers a matched-comparator architecture
  decision. No production changes, staging or commits in this review.
- **User-run test needed now:** no; commands follow the new guide's validation gate.

### Codex (GPT-5) — 2026-09-15 — contract pass and exact-link no-go

- **Increment 1 hypothesis:** the next candidate must be compared on a declared,
  content-identical matrix and must satisfy the written per-seed contiguity gate.
- **Implemented contract:** schema-v2 manifests now declare expected cases before
  execution. Case records hash ordered read/molecule IDs, sequences, qualities,
  references, private truth and supplied damage profiles. The gate requires
  generator/evaluator/source identity, distinct assembly-source hashes, complete
  cases, reference-base agreement and clean/damaged per-seed N50 criteria. It
  reports case combinations and distinct reference seeds separately.
- **Contract verification:** 25 focused tests passed. A 2x synthetic smoke wrote
  six expected and six observed cases, all four case hashes and 24 source entries.
  Full suite: 334 run, 333 passed, one skipped under Python 3.12 with pysam 0.24.0.
  Declared Python >=3.13 is still unverified. No rejected large matrix was rerun.
- **Increment 2 hypothesis:** baseline layouts may retain exact contig dovetails
  supported by at least two uniquely assigned raw molecules. Frozen settings:
  seed 199001; clean unique 20x and damage unique 10x; overlap >=30; exact only;
  raw support >=2; existing anchors/search thresholds. Truth never entered link
  discovery. Coordinates are 0-based and intervals half-open.
- **Clean result:** 416 reads produced 343 baseline contigs, N50 132. All 2,512
  exact physical dovetails were raw-length feasible (overlaps 30–192, median 66).
  Of 365 bridge-candidate molecules, 354 supported competing discovered edges.
  Eleven edges had one unique supporter; none had two. Zero joins/paths survived;
  projection remained 343 contigs and N50 132.
- **Damaged result:** 208 reads produced 163 baseline contigs, N50 145. There were
  400 exact, raw-length-feasible dovetails (30–181, median 64). Of 149 candidate
  bridge molecules, 133 were ambiguous; 11 edges had >=1 unique supporter and
  three had >=2 (maximum three). Three paths reduced output to 160 contigs and
  changed N50 only 145→146; longest projected contig was 213 bp.
- **Interpretation:** raw read length is not the immediate blocker because every
  exact proposal passed the coarse overlap+two-flank length check. The limiting
  evidence is unique bridge assignment/support. Candidate-local ambiguity is not
  a genome-wide uniqueness proof. Projection coordinate safety is intentionally
  not evaluated because the existing linker does not compose raw placements.
- **Disposition:** Increment 1 PASS; Increment 2 NO-GO. The clean 1.0x result cannot
  reach the frozen >=1.5x requirement, so the provenance adapter and development
  matrix are NOT AUTHORIZED. Do not lower support, enable near-exact matching,
  replace seeds or start another correction loop.
- **Evidence/source:** `results/post-layout-link-opportunity-dev-v1/`; gate SHA256
  `4b57c105586952dc90dd06e409cc3de979d3fbe70639b4b661ef22affcbefa2b`;
  validation SHA256 `a74412b424ef563e03bc0651305eb7cb677021b85d98a92c6b7bdc4733d4c92a`;
  opportunity script SHA256 `24ad7b803745dd8b8a8b1ff44ffe81be7b0a0403a858010ff9cab07ea08b6c93`.
  The restored assembly graph remains unchanged at
  `5bae98bb7dc5fb4fcedc3619907db82860c63be5d3e331c1f762395c094a828f`.
- **Next action:** matched, pinned CarpeDeam assessment followed by a concrete
  retain/integrate/narrow architecture decision. CarpeDeam is not installed in
  the checked WSL environment, and no conda/mamba executable was found. Do not
  invent or run actual-file commands until its CLI/version and exact input cohort
  are verified. User-run actual FASTQ needed now: no.

> **Reviewer correction, Astra, 2026-09-13:** follow the new active checkpoint in
> `sol_implementation_guide.md` and `astra_review_20260913.md`. Historical reports
> below remain intact. The 50k graph totals are 213 transitive removals and 3,768
> reciprocal edges (81.75%), not 33/3,981. The 178-bp proxy component is not a
> dataset-wide connectivity ceiling. Maximum-overlap selection is suspended.
> Gargammel DOI is 10.1093/bioinformatics/btw670; PyDamage is 10.7717/peerj.11845.

## Reviewer findings — Astra, 2026-09-13 (separate from implementation reports)

- **Evidence:** `results/astra-review-20260913/audit.py`, `audit.json`, source
  hashes and pre-edit document snapshots; detailed critique in the linked review.
- **Measured:** a redundant three-hop-reachable strict overlap makes four clean
  reads emit 120/100/100 bp; the chain-only control emits one 160-bp contig using
  the unchanged path rule. Clean seed199001 boundary-edge counts increase
  0/4/22/124 oriented edges at 2/5/10/20x. Neither result quantifies all ladder loss.
- **Verification:** 321 unit tests ran with pysam0.24.0: 320 passed, one external
  aligner integration skipped. Python3.12 is below declared >=3.13 support, so
  supported-runtime validation is still due. Standalone synthetic audit passed.
  Saved 50k FASTA equality was independently verified by streaming SHA256.
- **Interpretation:** recent work established operational progress and real
  fragmentation, but over-attributed the cause to path selection and extrapolated
  from the 1.156% mapping subset. The depth-aware comparison gate also needs work:
  its current case key lacks depth. Do not claim the new ladder passed that gate.
- **Disposition:** documentation/review only; assembler behavior unchanged.
  No real FASTQ run, staging or commit. Retain previous failed attempts as evidence.
- **Next action:** follow active Phase A, then one repair under Phase B gates.
  Implementation agents append their own findings below with model/date/source,
  measured results versus inference, frozen acceptance and explicit disposition.

## Implementation findings and historical reviewer entries

### Codex (GPT-5) — 2026-09-13 — Phase A graph-invariant trace

- **Hypothesis:** separate boundary containment, missed transitive witnesses and
  branch-preserving path emission before changing production graph policy.
- **Source/evidence:** `results/phase-a-invariant-trace-dev-v3/summary.json` and
  its source snapshot; clean development seed 199001 at target 20x. Truth is used
  only by diagnostics. Coordinates are 0-based, intervals half-open.
- **Tests:** focused layout-trace and paired-gate modules passed 16/16; full
  discovery passed 327 tests with one external-aligner integration skip. The paired
  gate now includes `target_coverage` in case identity and reports it per row;
  missing depth cases and duplicate same-depth cases are tested. It joined the
  saved ladder as 48 unique cases and 192 correlated mode/stage rows. Tests used
  WSL Python3.12 plus isolated pysam0.24.0; declared Python>=3.13 remains unchecked.
- **Measured baseline:** 416 reads/nodes, 3,036 physical graph edges, 767 after
  reduction, 743 branch-rejected, 24 reciprocal; 343 contigs, N50 132, longest 200.
- **Boundary ablation:** removed 62 physical boundary relations; contigs 346,
  N50 131, longest 200, and joined truth dovetails unchanged at 22. Rejected.
- **Redundancy audit:** all 690 strict residual branch physical edges had no
  sequence-equivalent witness within eight edges and 10,000 explored states; no
  search limit was hit. The ablation took zero actions and is rejected as a no-op.
- **Alternative-preserving diagnostic:** maximal non-branching paths produced
  732 contigs, N50 166 and longest 276, with 729 truth-labelled joins. It reused
  314 nodes across paths. The +25.8% N50 misses the frozen +50% gate, output exceeds
  input reads, and current fixed-slot provenance cannot represent it safely.
- **Finding:** 603/690 residual truth-labelled branch edges touch an ambiguously
  contained node. This association, together with prior strain failures from
  containment suppression, narrows the problem to alternative ownership rather
  than a missing longest-edge tie-break or broad transitive reducer.
- **Disposition:** diagnostic and gate changes retained; all three production
  topology candidates rejected. No actual FASTQ, new seed, staging or commit.
- **Next action:** specify and test a bounded coordinate-compatible ownership
  representation that does not suppress nodes or duplicate every branch path.
  Stop before a ladder if raw molecule ownership is ambiguous in tiny strain/repeat
  controls.

> **Astra review correction, 2026-09-09:** historical entries below are preserved
> as reports of their authors' work, not all endorsed conclusions. The final
> review entry corrects the claim that a no-op proves whole-input linkage is
> absent and records the clean-layout priority. Follow the current guide when
> older next-step wording conflicts.

## Reporting contract

Sol: append a dated entry after each implementation or experiment, even if it
fails or changes nothing. Do not edit the Astra review instructions in
[the implementation guide](sol_implementation_guide.md) to promote your findings.
Flag disagreements for review. Never delete a failed run or overwrite an earlier
result directory. Keep incomplete runs labeled incomplete.

Each entry must contain:

- Author/model, hypothesis, scope, baseline, branch and source hashes/snapshot.
- Files changed and the precise rule implemented; coordinate conventions.
- Exact commands, interpreter/dependency versions, test totals/passes/skips/errors.
- Result paths, seed roles, completed/expected case counts and matching-ID checks.
- Per-case before/after N50, correctly aligned span, mismatch numerator and
  denominator, mosaic/incompatible counts AND bases, nonambiguous coverage,
  unresolved bases, runtime and peak RSS where measured.
- Observed result versus interpretation; explicitly state unmeasured quantities.
- Gate pass/fail/neutral/incomplete, retained/reverted code, limitations and one
  next action. Equal summary metrics do not establish identical FASTAs.

## Historical findings transcribed by Astra — 2026-09-09

These entries are reviewer summaries of saved artifacts and this task's history,
not newly authored Sol reports. Attribution of all uncommitted changes to Sol is
not possible because there is no clean commit boundary.

| Implementation | Evidence and finding | Review status |
|---|---|---|
| Exact overlap graph | Clean unique synthetic controls demonstrated long, correct spans; damaged reads remain fragmented. | Useful permanent control; not damage-aware. |
| Damage-compatible graph edges | v5 development roots each contain 72 cases (36 damaged), as do their held-out counterparts. Terminal-directional acceptance can reconnect reads but uses positive profile values categorically. | Development evidence; not quantitative edge calibration. |
| Variable-length layout and containment | Raw-coordinate, containment and transitive-spelling repairs are documented in the guide. Reciprocal selection and transient linking had reported strain/recovery regressions. | Retain tested invariants; do not revive rejected heuristics without a new mechanism. |
| Repaired truth scorer | Distinguishes provenance ambiguity, candidate sequence agreement, injected damage and diagnostic mosaics; reports evaluated/unresolved denominators and repeat-masked coverage. | Necessary evaluation infrastructure, not an assembly improvement. |
| v5 consumed validation | Saved 95202-strains-3-0.4: exact N50 80, mosaic 135 bp, nonambiguous coverage 2,189; damage-aware N50 111, mosaic 254 bp, coverage 2,062. One mosaic each hides +119 unsafe bases and -127 coverage. | Failed. Seeds 95201–95204 are consumed, not fresh validation. |
| Junction-local phase guard | Task history reports neutral development spot checks; consumed failure remained 254 mosaic bases, with 416 uninformative checks. Code and partial outputs were deleted in that turn. | Rejected; numbers are historical reports, no longer reproducible from those deleted outputs. Initial claim of unique-panel fragmentation was retracted after comparison with the proper damaged baseline. |
| Phase-block column helper | Current source and four tests reviewed. Only tests call it; no graph adapter or assembly integration exists. Earlier turn reported 271 tests passing. | Preliminary linkage helper; no demonstrated N50/accuracy gain. |

Preserved v5 roots under `results/`:

- `controlled-damage-string-evaluator-v5-dev`
- `controlled-string-evaluator-v5-dev-control`
- `controlled-damage-string-evaluator-v5-heldout`
- `controlled-string-evaluator-v5-heldout-control`

## Astra review probes — 2026-09-09

Executed small synthetic probes against current `phase_blocks` with existing
Q38 observation helper. Four allele pairings, each supported by two molecules,
yield one connected block and all four pair counts. This is connectivity evidence,
not unique phasing. A supported two-pair block becomes two separate blocks after
adding a singleton third allele to one site. Both behaviors require explicit
semantics/regressions before integration. See the authoritative guide checkpoint.

## Sol entries — append below

No new Sol-authored entry has been recorded in this section yet.

### Codex (GPT-6, Astra-guided continuation) — 2026-09-09 — phase-link evidence repair

- **Hypothesis and scope:** Supported same-molecule allele pairs can be retained
  without treating every connected variable site as a resolved haplotype. This
  increment changes only the preliminary phase helper and its synthetic tests;
  it does not alter graph construction, layout, consensus or emitted FASTA.
- **Baseline and snapshot:** branch `fix/support-two-quality-gates`; baseline
  `damage_string_graph.py` SHA-256
  `603c6e092f78d21827c5653ee6ed350dc026f69d0a74bfa2934acd99600eccc9`.
  Intermediate helper-only `phase_blocks.py` SHA-256 before the component adapter
  `48fe4861e1a1dbe1bf2ded5d8da52730f93db445b6a4d13c67afa47a82112301`;
  updated test SHA-256
  `2d18128766c51d2e1f5cbebc26642ce7ab4aa7551f1b45629938e2334376d654`.
- **Implemented rule:** contig positions are 0-based. Marker alleles require at
  least two robust, deduplicated molecules. Lower-support robust alleles remain
  in an exclusion audit and cannot veto an otherwise complete link. All variable
  site pairs are examined. Molecule IDs are retained for eligible marker alleles
  and every supported allele pair. Only complete one-to-one pair mappings join
  a block. Branching, incomplete and cycle-inconsistent mappings remain explicit
  ambiguous evidence. A connected block is not reported as a complete haplotype.
  The consensus guard and the two-molecule heuristic were not changed.
- **Synthetic controls:** singleton third-allele contamination, all-four-pair
  ambiguity, nonadjacent links, a three-site chain with disjoint molecule cohorts,
  an inconsistent three-site cycle, deduplicated conflicting `None` observations,
  reverse orientation, terminal damage and ordinary singleton error. The reverse
  test supplies contig-oriented allele indices and reverse-complemented raw bases.
- **Verification:** WSL Python 3.13.15 with `pysam==0.24.1` in
  `/tmp/anvaya-audit-tests`. Focused command
  `PYTHONPATH=src /tmp/anvaya-audit-tests/bin/python -m unittest tests.test_phase_blocks -q -b`
  ran 11/11 tests successfully. Full command
  `PYTHONPATH=src /tmp/anvaya-audit-tests/bin/python -m unittest discover -s tests -q -b`
  ran 278 tests successfully with one skip. An earlier full-suite attempt using
  `/root/miniforge3/envs/anvaya/bin/python` ran 276 tests with 21 skips and one
  import error because that environment lacks `pysam`; no assertion failed.
  `git diff --check` reported only pre-existing line-ending warnings.
- **Observed result and gate:** the helper-semantics gate passes its listed
  regressions. Assembly accuracy, N50, aligned span, mismatch denominators,
  mosaic/incompatible counts and bases, runtime, RSS, and real-data behavior are
  unmeasured because the helper still has no production caller. No result
  directory or seed set was created. No actual FASTQ was run. The change is
  retained pending review; the next action is the bounded overlap-component
  adapter fixture required by the guide.

### Codex (GPT-6, Astra-guided continuation) — 2026-09-09 — bounded component adapter

- **Hypothesis and scope:** An already laid-out derived component contains enough
  truth-blind raw-placement data to project molecule-deduplicated, damage-aware
  phase evidence while abstaining on repeated physical placements. This does not
  yet establish that a phase rule can safely choose or reject graph edges.
- **Baseline and snapshot:** branch `fix/support-two-quality-gates`; unchanged v5
  `damage_string_graph.py` SHA-256
  `603c6e092f78d21827c5653ee6ed350dc026f69d0a74bfa2934acd99600eccc9`.
  Current uncommitted `phase_blocks.py` SHA-256
  `9a573c56a6ccb3479b5a1a2d5661dae8bf60d2b448d38d98422fdcbe1edb992d`;
  component test SHA-256
  `475c5a3bb93aca9315bfa81f0d3f7b9a2e18d46b8e33fb08f0061dba2431879b`.
- **Implemented rule:** `phase_component` consumes a derived record, immutable
  pool records, 0-based half-open `RawPlacement` intervals and a supplied damage
  profile. It accepts no reference, origin or simulated truth. A read with more
  than one physical `(offset, orientation)` placement makes its molecule
  ambiguous and contributes no observations. Consistent placements project raw
  bases and qualities in both orientations, deduplicate by molecule and exclude
  conflicts. Diagnostics count consistent/ambiguous placements and molecules,
  usable observations, distinct observed molecules, covered/uncovered sites,
  missing and low quality, ambiguous bases, outside-contig bases and duplicate
  conflicts.
- **Private-truth fixtures:** a 12 bp component with two linked marker positions
  supported by four molecules on both strands recovers the expected two allele
  pairs and exact molecule IDs; a component whose two variable sites have
  disjoint molecule cohorts remains unlinked; a repeated read placed at offsets
  0 and 4 is reported twice as ambiguous and contributes zero evidence. Tests
  also verify current sequence and immutable raw reads directly, plus missing,
  low-quality, ambiguous, outside-coordinate and duplicate-conflict accounting.
- **Verification:** focused command
  `PYTHONPATH=src /tmp/anvaya-audit-tests/bin/python -m unittest tests.test_phase_blocks tests.test_phase_component -q -b`
  ran 16/16 tests successfully. Full command
  `PYTHONPATH=src /tmp/anvaya-audit-tests/bin/python -m unittest discover -s tests -q -b`
  ran 283 tests successfully with one skip under WSL Python 3.13.15 and
  `pysam==0.24.1`. `git diff --check` reported only pre-existing line-ending
  warnings.
- **Observed result and gate:** the bounded projection fixtures pass. They use
  explicit synthetic derived records rather than output from the v5 graph, so
  graph-coordinate survival, pre-layout ambiguity and any assembly decision
  remain unmeasured. N50, aligned span, accuracy numerator/denominator,
  mosaic/incompatible counts and bases, runtime and RSS are unmeasured. No seeds,
  actual FASTQ or result directory were used. The next action is to drive a
  bounded fixture through `damage_string_graph.assemble`, verify its emitted raw
  placements against private truth, and only then define an integration rule.

### Codex (GPT-6, Astra-guided continuation) — 2026-09-09 — graph-coordinate gate

- **Hypothesis and scope:** v5 graph layouts should preserve raw-read coordinates
  and abstain when an exact child has multiple physical placements within one
  parent. This increment changes the string graph's handling of full-containment
  edges, adds graph-origin coordinate regressions, expands the variable-length
  manifest, and records the consumed 95202 diagnosis. It does not use phase
  evidence to accept, reject or spell a graph path.
- **Snapshot:** branch `fix/support-two-quality-gates`; current uncommitted hashes:
  `damage_string_graph.py`
  `200b421e4013dac3e3bb43563f36de9f9e39570032c3fbfd6b90933a6361a723`,
  `phase_blocks.py`
  `9a573c56a6ccb3479b5a1a2d5661dae8bf60d2b448d38d98422fdcbe1edb992d`,
  `test_string_layout.py`
  `ca3991ace18bd41881a2c4757635a175da39f54c653374c60d488cb21c9067b6`,
  `test_phase_component.py`
  `78e5682a144717f31a26868f30d154e9f0d53b814670f2277a98c25c2b70cf2e`,
  and `phase_component_diagnosis.py`
  `1b6761858390945e97a6f026bb534946388c7f3a325a90b31040f825ae5871d1`.
  The changed `variable_length_validation.py` manifest code has SHA-256
  `8b3a40ef2831594dc494006a010c6efa51928c3c0612c04abfc6a5fa9cb2110c`.
- **Private-truth graph fixtures:** a 120 bp unique chain assembled from reads at
  0-based offsets 0, 30 and 60, with the middle read reverse-complemented. The
  emitted sequence, orientations, offsets, raw reads and molecule IDs all match
  private expectations. A 40 bp non-palindromic periodic child occurs repeatedly
  in its 160 bp parent. Before the repair, the graph counted ambiguous containment
  but still emitted one merged path with one arbitrary child offset. The graph now
  removes only full-containment edges for children having multiple exact offsets
  or orientations in the same parent; both reads remain separate. Partial
  dovetails and children observed once in several parent reads remain eligible.
- **Rejected development attempts:** results are preserved as
  `variable-length-phase-coordinate-dev-v1` and `-v2`. V1 excluded every edge
  touching an ambiguously contained child and changed clean seed 95031 from the
  matched control's 153 contigs / N50 150 / longest 200 to 107 / 102 / 2358, and
  seed 95032 from 158 / 160 / 250 to 109 / 106 / 2359; five exact-versus-aware
  case gates failed. V2 limited exclusion to full-containment edges for every
  ambiguously contained child but still changed development outputs and failed
  four gates. Neither rule is retained.
- **Matched development result:** the retained same-parent-repeat rule is saved in
  `results/variable-length-phase-coordinate-dev-v3`. Its control, generated by
  temporarily removing only that rule from otherwise identical current source,
  is `results/variable-length-phase-coordinate-current-control`. Both contain the
  same 12 development case IDs from seeds 95031 and 95032 and exact plus
  damage-aware modes. Their complete `cases` payloads compare equal; only source
  manifests differ. The result is neutral on these cases, not a contiguity gain.
  The updated manifest hashes the experiment, scorer and all `src/anvaya/*.py`
  dependencies (23 files total). No older result was overwritten.
- **Commands and verification:** the development runs used
  `PYTHONPATH=src /tmp/anvaya-audit-tests/bin/python experiments/variable_length_validation.py --output <new-result-directory>`.
  The full unit command
  `PYTHONPATH=src /tmp/anvaya-audit-tests/bin/python -m unittest discover -s tests -q -b`
  ran 285 tests successfully with one skip under WSL Python 3.13.15 and
  `pysam==0.24.1`. `git diff --check` reported only pre-existing line-ending
  warnings.
- **Consumed-case diagnosis:** command
  `PYTHONPATH=src /tmp/anvaya-audit-tests/bin/python experiments/phase_component_diagnosis.py --output results/phase-component-diagnosis-95202-v1`
  completed and is explicitly labeled diagnosis, not validation. It reproduced
  one 254 bp diagnostic mosaic at `path_15`, with diagnostic positions 45, 95,
  145, 195 and 245 in 0-based contig coordinates. All 12 placements were
  consistent, all 254 sites covered, and 12 distinct molecules observed, but the
  final path had zero robust variable sites and therefore no phase links or
  blocks. This directly shows that layout selection removed the competing allele
  evidence before the post-layout adapter ran.
- **Gate and next action:** graph-origin coordinate survival passes, and the
  repeated-same-parent coordinate bug is repaired without changing matched
  development case payloads. The intended phase mechanism remains incomplete:
  post-layout integration would be a no-op on the known failure. The next bounded
  implementation must retain alternative graph nodes and solve their coordinate
  constraints before reciprocal path selection, reporting inconsistent components
  rather than choosing a spanning-tree offset. No new validation seeds or actual
  FASTQ were used; runtime and peak RSS were not measured.

### Codex (GPT-6, Astra-guided continuation) — 2026-09-09 — pre-layout phase audit

- **Failure addressed:** the consumed 95202 diagnosis showed that reciprocal
  path selection had already discarded competing alleles before the existing
  post-layout adapter ran. The retained implementation therefore observes the
  complete accepted overlap graph before containment retirement, transitive
  reduction, reciprocal-edge selection and path spelling. It remains
  output-neutral and does not yet split or choose a path.
- **Implemented rule:** `phase_graph.py` solves every accepted edge constraint as
  `target_offset = source_offset + shift` in a 0-based coordinate system. It
  reports each bidirected physical component once. Inconsistent cycles or both
  orientations of one physical node make a component unresolved; no arbitrary
  spanning-tree coordinate is treated as unique. Unique exact containments add
  explicit forward and reverse coordinate constraints. All consistent
  alternative nodes project their immutable raw reads, qualities, orientations
  and molecule IDs into `phase_placements`. Palindromic reads remain ambiguous.
  Graph connectivity itself never links markers; only the same molecule observed
  at both positions can do so.
- **Synthetic controls:** a parallel linked branch on both strands retains two
  supported allele pairs; disjoint marker cohorts joined by a conserved bridge
  remain unlinked; terminal C-to-T damage and one ordinary substitution error do
  not create robust phase markers. Contradictory coordinate cycles, physical-node
  orientation conflicts, bidirected mirror deduplication, palindromic ambiguity,
  and reverse unique containment are covered. All positions are 0-based and raw
  intervals are half-open.
- **Development observation:** preserved final-source result
  `results/variable-length-prelayout-phase-audit-dev-v2` uses the existing
  development seeds 95031 and 95032. Across 12 cases, the exact mode reported
  191 components (189 coordinate-consistent), 40 variable sites, 9 resolved
  links, 18 ambiguous links and 266 unlinked site pairs. Damage-aware mode
  reported 89 components (87 consistent), 38 variable sites, 9 resolved links,
  17 ambiguous links and 316 unlinked pairs. All damage-aware unique and clean
  cases had zero variable sites; strain cases retained 38 sites and all 9
  resolved links. Each repeat case had one coordinate-unresolved component.
  These are opportunity counts, not evidence that a path rule improves assembly.
  The earlier `-v1` result is also preserved; `-v2` followed a mirror-component
  deduplication tightening and has the same aggregate counts.
- **Output-neutrality and verification:** every layout, consensus score and
  consensus diagnostic in all 24 seed/scenario/mode comparisons exactly matched
  `variable-length-phase-coordinate-dev-v3`. Focused phase/adapter/layout tests
  ran 43/43 successfully. The full WSL suite using Python 3.13.15 and
  `pysam==0.24.1` ran 296 tests successfully with one external-tool integration
  skipped. `git diff --check` reported only pre-existing line-ending warnings.
- **Hashes recorded by the final-source v2 result:** `damage_string_graph.py`
  `e2dbb70dd604f7f6aef81c1130a0576a39ea3ae5f2108daf7114297971e8cb45`,
  `phase_blocks.py`
  `fc59f7f080fe338cfaa86fb11c718f28b2026dbf8aa07429cee63325bbdd5d14`,
  `phase_graph.py`
  `cac5e82078ca31970f2f5bd8dc44994b22e25933e59a37a05d1f6c06ce666d93`,
  `test_phase_graph.py`
  `e1aee5158f083d7f337451dd9a01c8978c554b86c3cb7d9777eac14fbb2ebdc9`,
  and `test_phase_component.py`
  `beeb2ee19a47651e79d89bbf770efaa5ba6bbc2332c8a651babccfad9416dc87`.
- **Gate and next action:** available pre-layout linkage is now demonstrated on
  bounded fixtures and measured on development inputs. The next single action is
  a conservative path-splitting fixture: an unlinked conserved bridge must split,
  a same-molecule linked bridge must remain joined, and coordinate-unresolved
  components must abstain. Only after that rule passes the frozen development
  safety/coverage/N50 gate is a real FASTQ smoke test informative. No actual
  FASTQ, new validation seed, staging or commit was used.

### Codex (GPT-6, Astra-guided continuation) — 2026-09-09 — rejected unlinked-block path split

- **Tested hypothesis:** if a selected reciprocal path has robust marker evidence
  exclusive to each side of a junction and those markers belong to separate phase
  blocks, remove one physical path edge near the marker midpoint. The exact mode
  remained untouched. The implementation was deterministic across forward and
  reverse orientations, never added an edge, and abstained on coordinate conflicts.
- **Bounded result:** synthetic fixtures proved the intended mechanics: an
  unlinked conserved bridge removed exactly one physical edge, a bridge supported
  by the same molecules remained joined, and a coordinate-unresolved component
  abstained. An initial focused run caught a disabled-callback bug and a
  strand-dependent cut tie; both were corrected before evaluation. The corrected
  focused suite ran 46 tests successfully, followed by 299 full-suite tests with
  one external-tool integration skipped.
- **Frozen development result:** preserved directory
  `results/variable-length-phase-split-dev-v1`, using only development seeds 95031
  and 95032. The rule split one damage-aware strain path in each seed and changed
  no clean, unique, profile-mismatch or repeat case. It reduced neither diagnostic
  mosaic bases nor incompatible-candidate bases and added zero nonambiguous unique
  reference coverage. Instead, splitting duplicated 35 overlap bases for 95031
  and 83 for 95032. Seed 95032 N50 fell from 157 to 155; 95031 N50 was unchanged.
  Evaluated and resolved bases rose only by the duplicated overlap lengths, with
  no mismatch reduction. All exact-mode layout, consensus and consensus-diagnostic
  payloads matched the prior control.
- **Conclusion:** the development gate failed. Separate phase blocks encode lack
  of demonstrated linkage, not positive evidence that the selected alleles are
  incompatible. Splitting every such path can fragment a correct strain path and
  inflate summed output length through repeated overlap sequence. No threshold
  sweep is justified by this result.
- **Disposition and next rule:** all output-changing code and its temporary tests
  were reverted. The pre-layout observation oracle and its audit diagnostics are
  retained. After the revert, focused tests ran 43/43 and the full suite ran 296
  tests successfully with one skip. The next bounded candidate must require a
  positive contradiction between selected path alleles and supported
  same-molecule phase pairs; unlinked and coordinate-unresolved cases must abstain.
  Retained hashes immediately after this rejection were `damage_string_graph.py`
  `e2dbb70dd604f7f6aef81c1130a0576a39ea3ae5f2108daf7114297971e8cb45`,
  `phase_blocks.py`
  `fc59f7f080fe338cfaa86fb11c718f28b2026dbf8aa07429cee63325bbdd5d14`,
  `phase_graph.py`
  `de327f955206733cab84f2dcac0deaa9a9ce5b2ab6d0df3e6b21862a879cf866`,
  and `test_phase_graph.py`
  `e1aee5158f083d7f337451dd9a01c8978c554b86c3cb7d9777eac14fbb2ebdc9`.
  No actual FASTQ, validation seed, staging or commit was used.

### Codex (GPT-6, Astra-guided continuation) — 2026-09-09 — neutral positive-contradiction veto

- **Hypothesis:** the unlinked-block failure might be avoided by acting only on
  positive evidence: split a selected path when its robust marker allele pair is
  absent from an already-supported same-molecule phase link. Unlinked sites,
  coordinate conflicts, ambiguous selected bases and markers changed by an edge
  correction abstained. Exact mode remained the unchanged control.
- **Bounded verification:** fixtures covered a supported incompatible pair, a
  supported compatible pair, unlinked markers, a corrected marker, a
  coordinate-conflicting component and forward/reverse edge identity. The focused
  phase suite passed 15 tests; the full suite passed 301 tests with one external
  integration skipped before development evaluation.
- **Frozen development result:** preserved directory
  `results/variable-length-positive-phase-dev-v1`, using development seeds 95031
  and 95032. The rule found zero contradictory selected paths and split zero
  edges in every exact and damage-aware case. Consequently every N50, longest
  contig, mosaic/incompatible count and base total, nonambiguous unique coverage,
  evaluated base denominator and mismatch count was unchanged. The exact control
  payloads also matched. This is a verified no-op, not an improvement.
- **Interpretation:** the nine resolved links observed by the pre-layout audit are
  compatible with all selected development paths. The problematic conserved
  bridges occur between unlinked blocks, where current molecules contain no
  truth-blind evidence for choosing or rejecting a strain combination. Lowering
  support would manufacture linkage rather than recover it.
- **Disposition:** the positive veto, its generic graph callback and temporary
  fixtures were reverted. The pre-layout phase audit remains. After reversion,
  43 focused tests and 296 full-suite tests passed, with one skip. Current hashes
  are `damage_string_graph.py`
  `e2dbb70dd604f7f6aef81c1130a0576a39ea3ae5f2108daf7114297971e8cb45`,
  `phase_blocks.py`
  `fc59f7f080fe338cfaa86fb11c718f28b2026dbf8aa07429cee63325bbdd5d14`,
  `phase_graph.py`
  `da9a1fc09b1f32ba5e9191c17f007036048bb74551b4eae2a3d55cd6dc1f9981`,
  `overlap_graph.py`
  `521130aa735c6fe7d82b8e893ac6d659fbfd9ee582594ed47944d6b544206be7`,
  and `test_phase_graph.py`
  `17853e28be868be2a74aaa8dcecf5e3088f73ff9fe48bfb924c88a2b41772f80`.
  The phase-splitting route is closed for the available evidence. The next work
  targets truth-known overlap candidate recall and safe extension, which can
  improve contiguity rather than repeatedly correcting an unidentifiable join.
  No actual FASTQ, validation seed, staging or commit was used.

### Astra — 2026-09-09 — source/results review and revised development priority

- **Scope:** documentation review and bounded synthetic diagnostic only; no
  assembly policy change. No separate active Sol task was identified. Existing
  entries identify their authors; the dirty tree cannot be attributed wholesale
  to Sol. The prior guide and findings are preserved in
  `results/astra-guide-review-20260909/`.
- **Independent verification:** the current full suite ran 296 tests in 3.658 s:
  295 passed, one external remapping integration skipped. Command:
  `PYTHONPATH=src /tmp/anvaya-audit-tests/bin/python -m unittest discover -s tests -q -b`
  from `/mnt/e/Projects/anvaya`. WSL required tool escalation after its sandboxed
  call returned access denied. No actual FASTQ was executed.
- **Saved-result cross-check:** `review.py` joins three result sets by unique
  `(seed, scenario)` keys and asserts identical 12-case membership. The positive
  phase run has identical complete layout, consensus and consensus-diagnostic
  payloads to pre-layout audit v2 for both modes. The unlinked split changes those
  three payloads only for damage-aware strain cases 95031 and 95032. Equality of
  these payloads does not establish FASTA byte equality.
- **Reproducibility gap:** each of the three recent result directories contains
  only `summary.json`. Current manifests cover driver, scorer and source modules,
  correcting the old guide's four-file claim. But hashes do not recover removed
  implementations or temporary tests. Preserve recoverable snapshots for future
  rejected experiments; do not claim past source was recovered by this review.
- **Interpretation correction:** the final-path 95202 diagnosis is not whole-input
  analysis. The positive veto was neutral on 95031/95032 and its summary does not
  expose every abstention reason. It therefore does not establish that all nine
  links were actionable checks, that the consumed 95202 failure was tested by this
  rule, or that the entire input lacks usable phase evidence. Earlier statements
  that the phase route was conclusively closed overstate the evidence. Defer it
  for lack of demonstrated benefit, not proven impossibility.
- **Clean containment probe:** existing fixture seeds, no damage and no sequencing
  errors; exact substring placements enumerated in both orientations. Coordinates
  are 0-based, half-open. Each seed has 208 canonical sequence nodes. Seed 95031
  has 35 children with one containment placement, 105 with multiple placements,
  and 68 with none; 95032 has 27, 108 and 73 respectively. Saved clean layouts
  have 153/158 contigs, N50 150/160, longest 200/250, 6/4 pre-layout components
  and 127/133 ambiguous ends. These are measured associations, not a proof that
  removing contained nodes improves assembly. Multiple parents may represent one
  locus; repeats/strains prevent arbitrary parent selection.
- **Guide corrections:** removed contradictory active instructions to split
  unlinked blocks; marked helper/coordinate repairs completed; updated scorer,
  adapter ownership, manifest coverage, cap semantics and test counts; separated
  operational actual-input smoke testing from scientific performance validation.
  Every old guide section has an explicit reconciliation entry.
- **Literature:** checked primary records and official repository documentation for
  string graphs, CarpeDeam, POLYTE and contained-read preservation/RAFT. Together
  these support auditing layout redundancy while retaining haplotype evidence;
  they do not justify a guaranteed N50 gain or blanket containment removal. Links
  and transfer limitations are in the guide. Some full-text endpoints were blocked;
  this was not a line-by-line audit of external implementations.
- **Next action:** one bounded full-graph loss trace on existing clean variable-length
  controls, with a smallest forward/reverse reproducer. Rank actual lost correct
  span, including containment and reciprocal branching, before selecting a minimal
  change. Then test that change against exact/damage-aware controls and repeat/
  strain negatives. No new phase splitter, correction loop or validation seeds now.
- **Artifacts:** `review.py`, `review.json`, previous documentation and retained
  source snapshots under `results/astra-guide-review-20260909/`. The probe records
  source hashes. No staging, commit, source-policy edits or actual-data run.

### Codex (GPT-5, Astra-guided continuation) — 2026-09-09 — clean loss trace and rejected coordinate containment

- **Increment A:** added `experiments/layout_loss_trace.py` and its focused tests.
  It exhaustively classifies proper clean dovetails through candidate discovery,
  containment, two-hop reduction and reciprocal selection on development seeds
  95031/95032. Coordinates and truth intervals are 0-based, half-open. Truth is
  evaluator-only. Result `results/layout-loss-trace-dev-v3` includes source
  snapshots, hashes, environment/resource metadata and a four-read forward/reverse
  reproducer.
- **Measured loss:** every selected edge was truth-compatible. Candidate discovery
  missed 40/1072 and 26/1104 proper dovetails. Reduction left 247/268 proper edges
  branching; 189/208 touched multiply contained nodes. Of 140/135 contained nodes,
  105/108 had multiple exact parent representations of one truth locus. A private
  truth-only ablation reduced 153/158 contigs to 2/1 and raised N50 150/160 to
  2358/2359. This established causal clean topology loss, not a production rule.
- **Candidate hypothesis:** candidate parents that imply one coordinate in a
  consistent accepted graph might allow the contained node to retire safely while
  projecting its raw molecule through all compatible parents. Forward/reverse,
  repeated-coordinate and disconnected/different-coordinate fixtures were added.
  Focused tests passed; the full pre-matrix suite ran 301 tests with one external
  remapping integration skipped.
- **Frozen development outcome:** preserved under
  `results/variable-length-coordinate-containment-dev-v1`, including source/test
  snapshot, gate, matrix, stdout/stderr and `/usr/bin/time` output. Clean exact and
  damage-aware contigs changed 153→2 and 158→1, with N50 150→2358 and 160→2359.
  Repeat cases were unchanged. The run took 4.94 s and 43,068 KiB maximum RSS.
- **Gate failure:** seed 95031 strain unresolved bases increased 0→650 in exact
  and damage-aware modes. Seed 95032 damage-aware increased 237→485. These were
  cross-reference transforms. Evaluated mismatch rates also worsened in 16/24
  comparisons, partly because removing duplicated correct output shrank the
  denominator. The independent strain-provenance failure is sufficient to reject.
- **Disposition:** output-changing source and temporary candidate tests were
  reverted; increment-A diagnostic and results remain. Measurement: coordinate
  equivalence resolves topology. Inference: it does not establish unique strain
  ownership. The next candidate must separate topology suppression from evidence
  ownership, attaching raw evidence only after unique output ownership or retaining
  it explicitly unresolved and outside allele counts. No validation seeds, actual
  FASTQ, staging or commit were used.

### Codex (GPT-6, Astra-guided continuation) — 2026-09-10 — deduplicated evaluator and rejected exact-edge owner

- **Evaluator hypothesis:** emitted-base mismatch rate rewards duplicate correct
  substrings and cannot by itself judge containment removal. `score` now projects
  uniquely resolved outputs onto `(reference, 0-based position)`, complements
  reverse-strand observations, evaluates concordant output once and reports
  discordant outputs separately. Truth remains evaluator-only. Duplicate-correct,
  duplicate-error, conflict and reverse-coordinate regressions pass.
- **Rescore finding:** across 48 mode/stage comparisons for the unique-owner
  candidate, the legacy rate worsened 28 times while the nonduplicated mismatch
  rate worsened zero times. Six rows, representing two distinct damage-aware
  inputs, nevertheless gained one discordant coordinate and lost one evaluated
  coordinate. This keeps the candidate rejected without relying on the biased rate.
  Artifacts are under `results/variable-length-unique-owner-containment-rescore-v1`.
- **Localized evidence:** the new sites are reference positions 2143 (seed 95031
  `damage_error_unique`) and 2101 (seed 95032 `damage_unique` and its same-layout
  profile-mismatch twin), all 0-based. Each has one A-bearing output against several
  G-bearing outputs and persists through consensus. The A-bearing candidate
  mismatches occur at 0-based contig positions 3 and 106.
- **Exact-edge hypothesis:** coordinate equivalence inferred through a
  damage-corrected graph edge might create those conflicts. The frozen candidate
  therefore used correction-free accepted edges for coordinate inference while
  retaining unique post-projection ownership. Focused tests passed; the full suite
  ran 305 tests with one external remapping integration skipped.
- **Frozen outcome:** `results/variable-length-exact-coordinate-containment-dev-v1`
  retains source/tests, gate, matrix and resources (5.17 s wall, 65,644 KiB maximum
  RSS). Clean contigs still changed 153/158→2/1 and N50 150/160→2358/2359. The same
  two discordances remained. The candidate failed six of 48 mode/stage comparisons.
- **Root-cause trace:** baseline reads 21 and 16 occur in correction-aware two-read
  paths of 119/111 bp with zero mismatches. Suppression instead emits 97/109 bp
  singletons, loses partners 134/146 and retains the terminal damaged base. The raw
  base is high-quality T at 5-prime distance 3/2; reverse orientation yields A at
  reference positions 2143/2101. This is a path-spelling regression, not merely a
  denominator or transform-category change. Evidence is in the rescore directory.
- **Disposition and next mechanism:** exact-edge restriction was rejected and
  reverted, falsifying its causal hypothesis. The next truth-blind candidate should
  provisionally suppress containment but restore it when a baseline-selected
  correction-bearing adjacency disappears or changes relative coordinate. No new
  seeds, actual FASTQ, staging or commit were used.

### Codex (GPT-6) — 2026-09-10 — correction-adjacency candidate and failed validation

- **Reference boundary:** Myers and Simpson–Durbin support removing redundancy only
  when the represented string-graph walk is preserved. Kamath et al. demonstrate
  that contained-read removal can break a necessary haplotype path. Their long-read
  thresholds do not transfer here; the applicable invariant is path preservation.
- **Candidate:** compared unsuppressed and provisionally suppressed layouts, retained
  every selected correction-bearing adjacency at the same bidirected orientation
  and relative 0-based offset, and restored suppressions monotonically by affected
  physical component. Added a separate error-only development panel.
- **Development:** 305 tests ran with one external integration skipped. All 72
  normal/error-only comparisons passed. Clean N50 changed 150/160→2358/2359.
  Frozen evidence is under
  `results/variable-length-correction-adjacency-containment-dev-v1`.
- **Frozen validation:** previously unused seeds 96011/96012 were predeclared and
  consumed once. Seventy of 72 comparisons passed. Seed 96012 damage-aware
  `damage_error_strains_80_20` gained one discordant coordinate at reference 1,
  position 1539 (0-based), in layout and consensus. N50 improved 150→155 and
  nonambiguous recovery 4178→4341, but the per-case ambiguity gate rejects it.
- **Disposition:** candidate preserved in results and reverted without tuning on
  validation. Correction-bearing adjacency preservation fixes the known development
  failures but does not protect every strain alternative. The next diagnostic must
  classify position 1539 as an exact-path-spelling or provenance-resolution change
  before proposing another rule. No actual FASTQ, staging or commit was used.

### Astra reviewer — 2026-09-10 — independent audit and revised next phases

- **Scope:** reviewed current implementation and frozen containment attempts, paired
  result summaries, comparison scripts and guide. No assembler policy changed.
  The dirty tree contains several authors' work; this is not blanket attribution
  of all changes to Sol. Earlier entries remain historical evidence.
- **Reproduced defect:** the rejected correction-adjacency candidate canonicalizes
  reverse endpoints without transforming shift. For lengths 100/80 and overlap 60,
  forward/reverse shifts are 40/20; its signatures disagree for the same physical
  edge. `results/astra-review-20260910/audit.py` ran and saved `audit.json`.
  Historical rollback impact remains to be measured.
- **Inspection finding:** coordinate co-occurrence of accepted-edge endpoints is
  weaker than preservation of selected traversal and spelled correction alleles.
  A fixed-point stop also lacks an explicit final preservation assertion. Thus the
  preceding candidate's claim to retain every selected corrected adjacency is not
  established. An end-to-end regression is required before reusing the guard.
- **Result interpretation:** development 72/72 and validation 70/72 are correlated
  mode/stage rows over 18 cases per panel, not independent replicates. The failure
  is one seed/scenario at two stages. At reference index 1 position 1539 (0-based),
  baseline already has three uniquely projected covering contigs; candidate has
  four. Stable raw-provenance matching must determine whether the extra conflicting
  output changed spelling or merely became resolved. The historical rejection
  remains; the cause is not yet proven to be a new false join.
- **Evaluation gaps:** duplicate/equal-key checks passed independently on all four
  saved panels, but the comparison script does not enforce them. It also lacks
  a positive-benefit assertion and a failing exit status. Coordinate disagreement
  can hide mismatches from the denominator. Coverage is not exact correct span;
  the variable-length driver does not pass repeat masks. Explicit validation seeds
  are mislabeled development in summary metadata, and default held-out seeds have
  already been consumed.
- **Verification:** system Python ran 300 tests: 278 passed, 21 skipped and one
  import error because `pysam` is missing. Prior temporary environment is absent.
  This is incomplete suite verification, not a claimed code regression. No actual
  FASTQ or new assembly matrix was run.
- **Disposition:** updated reviewer guide with Phase 1 paired diagnosis/invariant
  tests, Phase 2 conditional single-mechanism changes and Phase 3 separate validation
  and operational tracks. New findings must use the author-labeled template in
  section 14. No automatic all-adjacency guard or additional polishing loop.
- **Next action:** complete guide section 11.1 for the consumed failed case, then
  use the result to select one implementation branch. No new confirmation seeds.

### Codex (GPT-5) — 2026-09-10 — Phase 1 evaluator and provenance diagnosis

- **Branch/source:** `fix/support-two-quality-gates`; retained graph source SHA256
  `f6101c667ef9f39cae4fac5131a990816f958e3736e4b1b9e48569677ed88a9e`;
  frozen candidate graph source SHA256
  `9ce3d6cca2f3fcf2890fb46e7f21e90649564444c4ebe372d4602c7c7a040264`.
  Imported paths and hashes are recorded in the trace and matrix summaries.
- **Hypothesis/single mechanism:** first determine whether validation position 1539
  is a new spelling error or an old damaged observation newly resolved by containment.
  No containment policy was enabled in the retained assembler. Evaluator changes
  expose paired coordinate states and supply simulator-known damage events/repeat
  masks; validation metadata now prevents consumed seeds being labeled fresh.
- **Predeclared criteria:** zero new incorrect/discordant/lost previously-correct
  coordinates; positive benefit in damage-aware mode on a damaged development input;
  error-only is a safety control and may be a no-op. Cases and correlated mode/stage
  rows are counted separately. Coordinates are 0-based, intervals half-open.
- **Tests:** focused tests: 54 run, 54 passed, 0 skipped/failed. Broad suite:
  314 run, 292 passed, 21 skipped, one import error because system Python lacks
  `pysam`; this is incomplete environment verification. No actual FASTQ was run.
- **Inputs and joins:** development seeds 95031/95032; consumed diagnosis seeds
  96011/96012; main and separate error-only panels. Case keys and per-coordinate
  keys are required unique/equal. Baseline/candidate trace input SHA256 values match.
- **Measured finding at 1539:** corrected scoring gives discordant A/G in both
  baseline and candidate at layout and consensus. The A is the Q35 simulated C→T
  observation from `molecule_202`, raw position 2. Baseline path 8 contains molecules
  111/202 and has cross-reference provenance; the known event resolves it to the
  reference-1 candidate, which already spells A. Candidate post-projection containment assigns
  molecule 111 to parents 58/164 and emits molecule 202 as a unique singleton.
  This changes membership and resolution type, but neither spelling nor coordinate
  state. The old evaluator made it appear newly resolved because it omitted events.
  Evidence: `results/containment-failure-diagnosis-phase1-v4`.
- **Corrected matrices:** development main has 12 cases/48 rows and zero safety
  failures, but fails the global gate because no damaged case improves in
  damage-aware mode. Development error-only passes 6 cases/24 rows. Consumed main
  passes 12 cases/48 rows and gains one correct coordinate plus N50 155→160 for
  seed 96011's damage-aware strain case; consumed error-only passes 6 cases/24 rows.
  These are diagnostic, not new confirmation. Evidence is in the two
  `correction-adjacency-rescore-*-phase1-v1` directories; use the v2 comparisons.
- **Actions:** the frozen candidate had 7 post-projection-owned containments in the
  traced case; molecule 111 was one. Its archived reverse-edge signature bug remains.
  The current source adds only a tested physical-edge signature: unequal 100/80 bp
  reverse shifts 40/20 and mirrored correction coordinates identify one edge.
- **Disposition:** Phase 1 diagnosis and prospective evaluator contract retained.
  Frozen containment candidate remains rejected because it provides no damage-aware
  development benefit and its preservation semantics are incomplete. The former
  validation discordance claim is superseded by the corrected paired diagnosis.
- **Next action:** follow section 12's localized-restoration branch. Build a minimal
  guard around actual selected traversal/spelling, use the corrected physical edge
  identity, and prove one disconnected safe containment can remain active. Do not
  restore the archived component-wide guard unchanged or consume new seeds.
- **User-run test needed:** no. This phase uses small synthetic fixtures only.

### Codex (GPT-6) — 2026-09-10 — Phase 2 selected-edge localized containment

- **Branch/source:** `fix/support-two-quality-gates`. Candidate evidence is frozen
  under `results/local-selected-edge-containment-dev-v3` and the separate error-only
  control under `results/local-selected-edge-containment-error-dev-v3`, including
  recoverable source and tests. Core SHA256 values are
  `damage_string_graph.py=39b6b0e577d4ded49851a1ce4fb73405b8680274d6e07fd1ea8ed630e8eb7e61`
  and `overlap_graph.py=8f33d1164747244d8939e62ee6323dd872fd750903c68fd059b0db1a9f3a14b4`.
- **Hypothesis and single changed mechanism:** allow a multi-parent exact child to
  retire only when every accepted parent implies the same coordinate and the child
  has one post-projection owner. Preserve correction-bearing physical edges from
  the actual reduced traversal. Restore a directly responsible suppressed endpoint;
  if direct responsibility cannot be isolated, restore suppressions only in the
  affected physical component. The final invariant fails closed to the unique-exact
  baseline. No identity threshold, correction model, phase splitter, or DBG changed.
- **Literature basis:** Myers's string graph formalizes overlap-graph reduction
  (doi:10.1093/bioinformatics/bti1114). ContainX shows that blanket contained-read
  deletion is not coverage preserving and motivates retaining only necessary
  contained reads (doi:10.1093/bioinformatics/btad124). Kamath et al. demonstrate
  assembly gaps caused by contained-read deletion (doi:10.1101/gr.279311.124).
  CarpeDeam treats terminal damage as position dependent and separates raw-fragment
  extension from corrected-contig merging (doi:10.1186/s13059-025-03839-5). These
  sources support preserving demonstrated path evidence; their long-read heuristics
  and thresholds were not transferred to these short aDNA fragments.
- **Predeclared criteria and seed roles:** zero new incorrect, discordant, or lost
  previously correct 0-based coordinates; at least one damaged development input
  improves in damage-aware mode; error-only is a safety control. Seeds 95031/95032
  remain development inputs. Mode/stage rows are correlated and are not counted as
  independent replicates.
- **Tests:** focused string-layout 30/30; validation tests 8/8, including a fixture
  where one correction dependency is restored while a disconnected safe containment
  remains suppressed. Full discovery ran 320 tests: 298 passed, 21 skipped, and one
  import error because this WSL Python lacks optional `pysam`. There were no code
  assertion failures. Actual FASTQ was not run.
- **Inputs and joins:** main panel 12 independent cases/48 correlated mode-stage
  rows; error-only 6 cases/24 rows. `paired_coordinate_gate.py` verified equal unique
  case and coordinate keys, matched modes/stages, complete states, and nonzero
  evaluated denominators against the corrected frozen Phase 1 baseline.
- **Results:** both comparisons pass with zero failed rows and zero global failures.
  Seed 95031 `damage_unique` gains three correct coordinates in damage-aware layout
  and consensus with no new error/discordance; `profile_mismatch_unique` repeats the
  same gain on the same simulated reads under a scaled profile, so this is one
  independent positive input, not four replicates. Clean unique N50 changes
  150→2358 for seed 95031 and 160→2359 for seed 95032. Error-only coordinate states
  are unchanged; unique N50 changes 153→156 and 155→166, while strain N50 changes
  153→157 and 155→163.
- **Actions:** across the main matrix's exact and damage-aware modes, 1395
  coordinate-equivalent containments were proposed, 510 remained suppressed, 30
  were restored, and 48 protected physical edges plus 1210 selected traversal edges
  were audited. Twenty-five of the restorations required affected-component scope;
  five were direct endpoint restorations. Whole-graph fallback count was zero.
- **Resources:** timed main development run used 8.34 s wall time and 393132 KiB
  peak RSS; error-only used 4.50 s and 223164 KiB on Python 3.12.3 / WSL2 Linux
  6.18.33.2. These small synthetic measurements do not establish 100k-read scaling.
- **Finding (measured):** the selected-edge guard repairs the clean containment
  blockage and yields a small positive damaged-data coordinate gain without a paired
  safety failure. Unlike the rejected guard, it does not become a main-panel no-op
  and it never falls back across the whole graph.
- **Interpretation:** this is sufficient to retain Phase 2 as a development
  candidate. The biological evidence is still narrow: one independent damaged read
  set improves, the simulator has no indels/PCR duplicates/real abundance spectrum,
  and some cases need component-local restoration rather than direct attribution.
- **Disposition:** retained as the Phase 2 candidate, pending fresh independent
  validation and bounded operational profiling. No claim of general aDNA assembly
  superiority is made.
- **Next action:** freeze genuinely unused validation seeds and run the guide's
  Phase 3 biological panel once, then profile bounded synthetic sizes before issuing
  an actual FASTQ smoke command. Do not tune this mechanism on those validation
  seeds after inspection.
- **User-run test needed:** no at this point. Actual FASTQ testing starts only after
  independent validation and resource profiling pass.

### Codex (GPT-6) — 2026-09-12 — Phase 3 validation rejection and operational profile

- **Branch/source:** `fix/support-two-quality-gates`. The frozen Phase 2 candidate
  uses `damage_string_graph.py=39b6b0e577d4ded49851a1ce4fb73405b8680274d6e07fd1ea8ed630e8eb7e61`
  and `overlap_graph.py=8f33d1164747244d8939e62ee6323dd872fd750903c68fd059b0db1a9f3a14b4`.
  Candidate source, tests and evaluator are under
  `results/phase3-candidate-validation-v1/source`. The active assembler was restored
  to the Phase 1 safe baseline after rejection; its profiled graph hashes are
  `damage_string_graph.py=d1297ff2844b003e2e57b237f7bb059abbe8dcddf241c295ed6a2295956c213e`
  and `overlap_graph.py=521130aa735c6fe7d82b8e893ac6d659fbfd9ee582594ed47944d6b544206be7`.
- **Hypothesis and single changed mechanism:** independently validate the Phase 2
  localized containment rule once, then profile the restored safe CLI without
  changing assembly policy. The only retained code additions are evaluator import
  provenance and optional stage timings; timings are excluded from deterministic
  graph decisions.
- **Predeclared criteria and seed roles:** seeds 197101–197104 were verified absent,
  declared `validation`, run once and then added to the consumed-seed guard. The
  biological gate permits zero new incorrect/discordant/lost previously-correct
  coordinates and requires one safe independent damaged benefit. Coordinates are
  0-based; intervals are half-open. The operational 1000-read stop gate was fixed
  at 120 seconds wall and 2,097,152 KiB RSS before profiling.
- **Tests:** focused string-layout and CLI tests ran 38/38. After restoring the
  baseline, string-layout ran 27/27, layout-loss trace 3/3 and validation 6/6. The
  final broad run discovered 315 tests: 293 passed, 21 skipped and one import error
  because optional `pysam` is absent; no assertion failed. Python compilation of
  the changed evaluator, fixture generator, graph and CLI passed. No full actual
  FASTQ dataset was run.
- **Inputs and joins:** the main validation has 24 independent cases and 96
  correlated exact/damage-aware layout/consensus rows. The error-only control has
  12 cases and 48 correlated rows. The paired gate enforced unique/equal case and
  coordinate keys. `results/phase3-independent-validation-manifest-v1.json` fixes
  source/evaluator hashes and limitations. The validator now records the actual
  imported package paths and hashes instead of incorrectly hashing the current tree
  during frozen-baseline runs.
- **Results:** main validation failed 10 rows and the global positive-benefit gate;
  error-only failed 8 rows. All row failures occur in uneven 80/20 strain cases.
  Seed 197101 exact layout/consensus gained 67 correct coordinates but made reference
  1 positions 1518 and 1525 newly discordant. Seed 197103 lost 17 previously correct
  reference 0 positions 1703–1719 in exact and damage-aware modes; damage-aware mode
  gained two other coordinates. Seed 197104 lost 74 previously correct reference 1
  positions 1693–1766 and worsened unresolved/mosaic evidence; its damage-aware mode
  gained two other coordinates. The error-only panel reproduced 18-coordinate and
  74-coordinate strain losses for seeds 197103 and 197104. Biological evidence is
  frozen under the four `results/phase3-*-validation-v1` directories.
- **Actions:** every validation gain was unsafe, so the candidate was rejected and
  not tuned on consumed seeds. Exact-mode failures establish that coordinate
  containment under strain ambiguity, rather than damage scoring alone, caused the
  regression. The Phase 2 implementation remains recoverable only in frozen results.
- **Operational result:** deterministic fixtures at 250/500/1000 reads contained
  30,755/56,041/113,791 bases. External wall/RSS were 0.72 s/23,972 KiB,
  0.92 s/27,068 KiB and 1.92 s/38,972 KiB. At 1000 reads, candidate discovery took
  0.654 s, containment 0.322 s, phase audit 0.183 s and consensus 0.092 s. A repeat
  produced byte-identical FASTA, placements and decision reports. Evidence is in
  `results/phase3-operational-profile-v1`.
- **Finding (measured):** the selected-edge containment mechanism is not strain-safe
  under independent validation. The restored baseline CLI is deterministic and
  comfortably within the declared resource gate at its current 1000-read cap.
- **Interpretation:** the prototype is ready for a bounded operational smoke test,
  but no containment-based contiguity improvement is retained. Literature-supported
  realism such as empirical fragment/error profiles, indels, duplicates,
  contamination and uneven community abundance remains outside the current narrow
  simulator. A full-data or superiority claim is unsupported.
- **Disposition:** Phase 2 candidate rejected; evaluator provenance fix, consumed
  validation seeds, stage timings and bounded fixture generator retained. Active
  assembly behavior is the prior safe baseline.
- **Actual-input smoke:** the user ran exactly 1,000 records from the existing 100k
  EMN001 FASTQ with the bundled `dhigh` profile. It completed in 1.35 s wall and
  26,504 KiB RSS. From 59,370 input bases it emitted 994 contigs/59,146 bases, N50
  64, N90 38 and longest 133. Four linear paths used eight reads; two reads were
  exact containments. Layout and consensus FASTA are byte-identical. All 59,144
  covered sites lacked the required molecule support, so zero bases changed. This
  passes the operational gate but is too sparse to test assembly or damage benefit.
- **Consumed-seed diagnosis:** at seed 197104, reference 1 position 1693 (0-based),
  baseline preserves `molecule_199` as a correct 143-bp `[1693,1836)` contig. The
  candidate joins it to `molecule_26` through an exact 46-bp overlap, producing a
  254-bp cross-reference mosaic with four unexplained mismatches against either
  strain. This makes positions 1693–1766 unevaluable. Suppressing multiply contained
  `molecule_0` and `molecule_111` before reduction removes their selected exact
  63-bp path, changes graph degrees and makes the unrelated 26→199 edge reciprocal.
  The correction-edge-only guard cannot see this nonlocal exact-path dependency.
  Frozen evidence is `results/phase3-strain-failure-diagnosis-v2`.
- **Next action:** close the pre-reduction containment branch. Replace its quadratic
  exact-containment enumeration so a larger bounded sample is measurable, but keep
  ambiguous contained reads in reduction. Only then develop one molecule/strain-
  evidence path-selection rule; do not preserve every prior exact path as another
  no-op guard or tune on consumed validation seeds.
- **User-run test needed:** no until bounded scaling tests pass. The completed smoke
  is recorded under `results/damage-string-graph-emn001-smoke-1000-v1`.

### Codex (GPT-6) — 2026-09-12 — indexed exact-containment scaling

- **Branch/source:** `fix/support-two-quality-gates`. Retained runtime source hashes
  are `damage_string_graph.py=3a2799220629b8b1fbfdbd52793cb1c3f50b1d0587762728e45de7fff4cd0553`,
  `overlap_graph.py=521130aa735c6fe7d82b8e893ac6d659fbfd9ee582594ed47944d6b544206be7`
  and `cli.py=72c527de361a0b4a51e4e41508fee4f237266ae8353c8f3573dd210e96d3b662`.
  Frozen source and results are under `results/indexed-containment-scale-v1`.
- **Hypothesis and single changed mechanism:** replace the quadratic exhaustive
  child-by-parent exact-containment scan with indexed prefix-seed discovery followed
  by full oriented-child verification. Containment policy, ambiguity handling,
  transitive reduction and path selection are unchanged. Offsets are 0-based and
  intervals half-open.
- **Literature basis:** Myers separates contained-read handling and transitive
  reduction in string-graph construction (doi:10.1093/bioinformatics/bti1114).
  Seed-and-extend indices use exact anchors to generate candidates followed by full
  verification (Minimap2, doi:10.1093/bioinformatics/bty191). Kamath et al. show
  that contained-read deletion can change assembly connectivity
  (doi:10.1101/gr.279311.124). These support optimizing candidate discovery without
  changing the rejected deletion policy.
- **Predeclared criteria and seed roles:** operational seeds 198201/198211/198221
  were absent before use. Indexed placements must equal an exhaustive oracle; four
  frozen 1,000-read outputs must be byte-identical. Sequential 2k→5k→10k profiling
  stops if wall time exceeds 120 s or RSS exceeds 2,097,152 KiB at 10k. This is an
  operational gate and supplies no biological acceptance evidence.
- **Tests:** focused containment/layout/CLI tests ran 42/42. Python compilation
  passed. Full discovery ran 316 tests: 294 passed, 21 skipped and one import error
  because optional `pysam` is absent; there were no assertion failures.
- **Inputs and equivalence:** the oracle fixture includes random reads, forward and
  reverse containment, reads below the 31-base maximum seed length and repeated
  placements. The indexed mapping has the same child keys and exact
  `(parent, reverse, offset_0based)` sets. A first implementation preserved data but
  reordered placement-report rows; sorting child keys fixed it before profiling.
  Layout FASTA, consensus FASTA, placements TSV and decisions TSV then matched the
  frozen 1,000-read outputs byte-for-byte. Contained/ambiguous counts remained
  92/293.
- **Results:** frozen-profile containment time changed 0.322→0.013 s (24.1-fold).
  The 2k profile used 2.77 s and 44,448 KiB; 5k conservatively used 6.54 s and
  84,260 KiB; 10k used 15.13 s and 132,160 KiB. At 10k, candidate edges used
  8.594 s, phase audit 2.171 s and containment 0.245 s. N50 stayed 151/150/150 at
  2k/5k/10k because the generator holds approximate coverage constant. A WSL timing
  inconsistency at 5k made application time exceed external wall time; the larger
  value is reported.
- **Actions:** every candidate placement is full-string verified. No occurrence cap
  is applied because dropping repetitive anchors would violate complete exact
  containment discovery. The rejected multi-parent containment behavior was not
  restored.
- **Finding (measured):** quadratic containment enumeration is no longer the scale
  blocker through 10k reads. Candidate-edge discovery is now dominant. Resource and
  deterministic-equivalence gates pass.
- **Interpretation:** an exactly 10k-record actual-input smoke is justified. Scaling
  beyond it requires profiling candidate-edge discovery; these synthetic timings
  cannot support a direct jump to 100k or a biological-improvement claim.
- **Disposition:** retained as a behavior-preserving performance change.
- **Next action:** user runs one exact 10k EMN001 smoke. If it passes resources but
  remains read-through, optimize only candidate-edge discovery before increasing
  coverage; path-policy work remains separately gated against strain mosaics.
- **User-run test needed:** yes; exact command is supplied after final verification.

### Codex (GPT-6) — 2026-09-12 — EMN001 10,000-read operational smoke

- **Input and scope:** exactly 10,000 records and 572,419 bases from the existing
  EMN001 subset, using the supplied CarpeDeam `dhigh` profile. This is an
  operational and overlap-density observation, not reference-backed biological
  validation. Coordinates in reports are 0-based and intervals are half-open.
- **Resources:** 5.57 seconds external wall time and 83,924 KiB peak RSS. The
  application graph stage completed in 5.01 seconds. Candidate-edge discovery
  used 2.241 seconds, phase audit 1.317, indexing 0.489, consensus 0.428 and
  indexed containment 0.077 seconds.
- **Assembly result:** read N50 was 60 bp and longest read 125 bp. The assembler
  emitted 9,343 contigs, 543,772 bases, N50 61 bp and longest contig 211 bp.
  It selected 351 linear paths containing 727 reads, retained 281 exact
  containments and removed eight transitive edges. The output-to-input count is
  93.43%, so the result remains mostly read-through.
- **Damage-aware evidence:** candidate classification recorded 1,377 exact, 24
  damage-compatible, 641 insufficient-molecule, 63 non-directional/internal and
  two low-quality candidates. Seven selected-overlap bases were corrected. This
  shows the damage-aware edge path is active, but does not establish accuracy.
- **Consensus:** layout and consensus FASTAs are byte-identical. Of 543,751
  covered sites, 542,246 had fewer than three independent molecules. The other
  1,505 already matched the model-best allele. There were no posterior, conflict,
  placement or linkage rejections and zero changed bases.
- **Decision:** the operational gate passes. Candidate discovery is the largest
  stage but is only 2.241 seconds, so optimizing it now would not improve the
  biological result. Do not tune the consensus thresholds: physical overlap
  density is the measured constraint. Run a bounded 25,000-record gate next,
  stopping above 60 seconds or 524,288 KiB. If output/read count remains above
  90% or assembly N50 is less than 1.10 times read N50, audit candidate recall and
  rejection categories on that frozen subset before raising coverage again.
- **Evidence:** `results/damage-string-graph-emn001-smoke-10000-v1` contains the
  input subset, outputs, reports, diagnostics, stdout/stderr, timing and this
  interpretation.

### Codex (GPT-6) — 2026-09-13 — provisional EMN001 25,000-read smoke

- **Measured result:** 25,000 reads and 1,425,995 input bases produced 22,702
  contigs and 1,328,030 bases. Read N50 was 60 bp; output N50 was 62 bp and the
  longest contig was 211 bp. External wall time was 16.56 seconds, application
  time 18.72 seconds and peak RSS 195,688 KiB.
- **Gate:** resources pass. Output/read is 90.81% and N50 gain is 3.3%, so both
  response gates fail and 50k is not authorized.
- **Validity issue:** `diagnostics.json` lacks the `candidate_search` and
  `raw_candidate_search` objects emitted by the current workspace assembler.
  The current `damage_string_graph.py` also has a later modification time than
  the result. The run therefore predates the audit instrumentation or imported a
  different package and cannot locate the first candidate-loss stage.
- **Verification:** current workspace `tests.test_string_layout` and
  `tests.test_cli` pass 39/39 under `PYTHONPATH=src`.
- **Next action:** rerun the frozen `reads-25000.fq.gz` under `PYTHONPATH=src`,
  recording `anvaya.__file__`, `damage_string_graph.__file__`, and the source
  SHA256. Require both candidate diagnostic objects before analyzing rejection
  stages. Do not change assembly policy or run 50k first.

### Codex (GPT-6) — 2026-09-13 — verified 25k candidate-stage audit

- **Source and reproducibility:** `PYTHONPATH=src` imported
  `/mnt/e/Projects/anvaya/src/anvaya/damage_string_graph.py`, SHA256
  `5bae98bb7dc5fb4fcedc3619907db82860c63be5d3e331c1f762395c094a828f`.
  The result reproduced 22,702 contigs. External wall time was 19.26 seconds,
  application time 18.78 seconds and peak RSS 195,328 KiB.
- **Candidate funnel:** 399,982 target anchors yielded 17,716 distinct tested
  offsets. Of these, 9,421 were selected, 4,209 were shorter than 30 bp, 1,744
  failed DNA identity and 2,342 failed R/Y identity. No anchor was occurrence-
  capped and no molecule had tied best placements. The 199,992 unavailable hits
  are the expected eight excluded self-anchor hits for 24,999 canonical nodes,
  not missing candidate evidence.
- **Damage classification:** 5,096 selected alignments were exact, 141 were
  damage-compatible, 3,145 lacked three independent molecules, 1,031 were
  non-directional/internal and eight lacked usable quality. Forty-one selected
  overlap bases were corrected.
- **Graph funnel:** 1,420 physical edges entered reduction, 1,387 remained after
  removal of 33 transitive edges, and 1,302 remained reciprocal after 85 branch
  rejections. Graph survival is 91.7%, so reduction and branching are not the
  dominant cause of read-through output.
- **Interpretation boundary:** aggregate rejection counts do not label which
  rejected offsets are true same-origin dovetails. The prior mapping proxy found
  565/610 expected directed overlaps selected (92.6%), but its
  `overlap_assembly.py` SHA256 differs from current source. It is suggestive only.
- **Next action:** rerun `experiments/overlap_opportunity.py` on the existing
  query-name BAM under current `PYTHONPATH=src`. If current proxy recall remains
  high, proceed to a bounded 50k coverage gate. If one rejection stage dominates
  true expected overlaps, develop one stage-specific candidate under synthetic
  strain/error safety gates first. Do not tune consensus or graph reduction.

### Codex (GPT-6) — 2026-09-13 — current-source overlap-opportunity audit

- **Source:** `overlap_opportunity.py=a2473c573dbe62b38d3c5d96e001741c918307c10a6f884470f54a4c89649390`,
  `overlap_assembly.py=580fab4033af1a0014b850e473132299dd456647b32bd473f6a54faec3b0ec35`,
  and `overlap_index.py=738fc5735e650edd4b422c3cffce53ba4c4541bf00fff28a8e61112d26d5e076`.
- **Result:** 565 of 610 reference-expected directed dovetails were selected
  (92.62%); 23 had no anchor vote and 22 failed R/Y identity. There were no
  DNA-identity or occurrence-cap losses among expected proxy overlaps. Runtime
  was 1.47 seconds and peak RSS 49,112 KiB.
- **Evidence limit:** only 1,264 of 97,964 query groups (1.29%) were eligible.
  Alternative/split mappings excluded 76,763, low/unknown MAPQ excluded 12,456,
  and 7,467 were unmapped. The result is strongly selected toward easy unique
  reads, uses mapping as a proxy rather than origin truth, and does not test graph
  ownership or final layout.
- **Decision:** candidate recall is adequate on the measurable subset and exactly
  reproduces the earlier 565/23/22 stage counts under current code. Do not relax
  anchors, DNA identity, R/Y identity, molecule support or graph reduction. The
  next experiment is an unchanged 50k coverage gate, predeclared at 120 seconds,
  1,048,576 KiB, output/input at most 90%, and assembly/read N50 at least 1.10.
  Failure stops escalation to 100k; success triggers reference-backed accuracy
  evaluation before any algorithm change.
- **Evidence:** `results/overlap-opportunity-100k-current-v2` and
  `results/damage-string-graph-emn001-smoke-50000-manifest-v1.json`.

### Codex (GPT-6) — 2026-09-13 — 50k response-gate failure and cohort audit

- **Source/resources:** the current `damage_string_graph.py` SHA256 is
  `5bae98bb7dc5fb4fcedc3619907db82860c63be5d3e331c1f762395c094a828f`.
  External wall time was 62.65 seconds, application time 67.30 seconds and peak
  RSS 326,504 KiB. Both 120-second/1,048,576-KiB limits pass.
- **Assembly:** 50,000 reads and 2,919,124 bases, read N50 62 bp and maximum 125
  bp produced 43,478 contigs and 2,636,221 bases, N50 65 bp and maximum 287 bp.
  Output/read is 86.96%, passing its gate. Assembly/read N50 is 1.048, failing
  the required 1.10 response. The predeclared decision stops before 100k.
- **Graph/consensus:** 4,609 physical edges entered reduction and 3,981 reciprocal
  edges remained after 33 transitive removals and 628 branch rejections (86.4%
  survival). The assembler accepted 1,674 damage-compatible candidates and
  corrected 373 overlap bases. Consensus changed zero bases; layout and consensus
  FASTAs are byte-identical.
- **Implemented audit support:** `experiments/overlap_opportunity.py` now accepts
  `--read-names-fastq` and streams the query-name BAM while retaining only a set
  of cohort names. It validates FASTQ structure, lengths and unique names; missing
  alignment names fail instead of shrinking the denominator silently. Summary
  fields distinguish source, requested and matched query groups. No assembler
  behavior changed.
- **Verification:** Python compilation passes. The focused test module could not
  execute in the agent's noninteractive WSL Python because optional `pysam` is
  absent; no assertion result is claimed. The user's data environment already ran
  the same pysam-based audit and will exercise the new exact-cohort path.
- **Next action:** user runs the exact 50k cohort audit. If eligible overlap
  components are intrinsically short despite high candidate recall, the dataset
  prefix lacks enough physical connectivity and another graph heuristic cannot
  solve its N50. If expected dovetails are lost at one stage, develop exactly one
  measured candidate and validate it against uneven-strain negatives.

### Codex (GPT-6) — 2026-09-13 — exact-cohort limit and coverage-ladder diagnosis

- **50k cohort audit:** requested and matched query groups are both 50,000. Only
  578 were eligible; they formed 535 >=30-bp reference-overlap components and the
  longest optimistic component was 178 bp. Candidate recall was 65/70 (92.86%);
  all five misses lacked an anchor vote. This proxy is selection-biased but shows
  that the uniquely measurable cohort cannot support long contigs.
- **Literature basis:** Gargammel models fragmentation, post-mortem damage and
  sequencing errors for truth-known aDNA simulations (Renaud et al. 2017,
  doi:10.1093/bioinformatics/btx013). PyDamage evaluated damage over explicit
  coverage bins (Borry et al. 2021, doi:10.1093/nargab/lqab096). metaSPAdes treats
  local coverage and strain/repeat ambiguity explicitly rather than applying one
  global coverage rule (Nurk et al. 2017, doi:10.1101/gr.213959.116).
- **Implemented evaluator change:** `variable_length_validation.py` accepts
  positive `--depths`; the existing 10x default is unchanged. In strain fixtures,
  total target coverage is split 80:20. The output records target coverage for
  every case. No assembler behavior changed. Focused tests pass 46/46.
- **Development ladder:** seeds 199001/199002, target depths 2/5/10/20, six
  scenarios and exact/damage-aware modes produced 48 cases. These seeds are now
  development evidence. Runtime was 30.99 seconds; peak RSS 1,467,144 KiB reflects
  retained per-contig evaluator evidence.
- **Root-cause result:** clean exact and damage-aware outputs are identical and
  error-free, yet mean N50 declines 237→180→149→133.5 bp as depth increases.
  Mean contigs rise 16.5→53→149.5→343. At 20x, 3,024.5 physical edges yield only
  20 reciprocal edges. The branchless reciprocal rule, not damage scoring or
  physical coverage, is the first demonstrated contiguity defect.
- **Damage result:** at 2x damage-only coverage, damage-aware mean N50 improves
  191→226.5 bp and discordant reference bases fall 11→8. Higher-depth benefit is
  suppressed by path selection. Damage-aware 5x strain unresolved bases are worse
  than exact (345.5 versus 308), so strain safety cannot be inferred from unique
  cases.
- **Next candidate/gate:** after unchanged transitive reduction, select each
  orientation's unique maximum-overlap outgoing and incoming edge and retain only
  reciprocal choices; tied best edges abstain. Require >=50% clean-20x N50 gain,
  zero clean error, and no increase in discordant, mosaic, incompatible or
  unresolved evidence across every repeat/strain development row before fresh
  validation. Do not tune thresholds after a failure.

### Codex (GPT-6) — 2026-09-13 — Phase B coordinate-owner candidate rejected

- **Author/date/source:** Codex GPT-6, 2026-09-13. No branch change. Baseline
  `damage_string_graph.py=5bae98bb7dc5fb4fcedc3619907db82860c63be5d3e331c1f762395c094a828f`;
  rejected candidate `a496dce3dfd9fd1d00897f916aedf691190fc407dbff146c1958fc27282de83d`.
  Frozen imported sources, commands and traces are under
  `results/phase-b-coordinate-ownership-dev-v1/`.
- **Hypothesis and one mechanism:** retire a multiply contained node only when
  all parent representations imply one 0-based graph coordinate, component
  spelling is correction-compatible, supported phase variation is absent, and
  provisional projection gives the node one output coordinate. Attach its raw
  placement once to that output.
- **Predeclared gate/roles:** development seeds 199001/199002, depths 2/5/10/20,
  six main scenarios plus separate error-only unique/repeat/80:20-strain controls;
  exact/damage-aware modes and layout/consensus stages. Require at least 50% clean
  20x N50 gain per seed, a damaged benefit, and zero per-case coordinate, error,
  discordance, unresolved, repeat or strain regression.
- **Tests/environment:** focused ownership and Phase A regressions passed 41/41.
  Full suite passed 333/333 under Python 3.12.3 and `pysam 0.24.0`, with one
  external-aligner test skipped. Python >=3.13 remains unverified.
- **Main-panel result:** the paired gate passed all 48 independent cases and 192
  correlated mode/stage rows. Clean 20x N50 changed 132→2398 and 135→2382;
  each candidate output had one contig and retained 2398/2382 uniquely recovered
  correct bases. Candidate runtime was 33.22 seconds and peak RSS 1,441,672 KiB,
  dominated by retained evaluator evidence. At least one damage-aware low-depth
  case improved, while repeat/strain rows in this panel had no gate failure.
- **Error-only result:** rejected. One independent case, seed 199001 2x 80:20
  strains, lost reference-1 interval `[2069, 2120)` (0-based, half-open) in exact
  and damage-aware layout; consensus repeated the layout loss. Evaluated reference
  bases fell 2616→2565 and the mismatch-rate denominator regressed. Baseline
  `node_33` was a correct 51-bp minor-strain output. The candidate assigned it to
  a 360-bp cross-reference output with five unexplained mismatches.
- **Actions in the failed case:** three multiply contained nodes were classified
  coordinate-compatible and all three received one provisional output owner.
  Eleven nodes were retired versus eight at baseline; contigs fell 21→17. The
  failing child had two exact reference-0 parent placements even though its raw
  origin was reference 1. The component had zero robust phase variables, so the
  required ownership evidence was absent rather than contradictory.
- **Finding:** one graph/output coordinate is insufficient evidence ownership for
  an exactly shared short fragment. The discriminating strain allele can lie
  outside the child, and no threshold over the child's sequence can recover that
  missing linkage.
- **Disposition:** rejected, source/results preserved, production baseline restored.
  This is Phase B mechanism failure 1. Do not tune or restore it.
- **Next action:** specify a minimal primary-versus-unresolved containment evidence
  representation before a second graph change. Preserve one nonduplicated raw
  owner without assigning an unlinked child to a primary contig. If the existing
  clean-benefit and coordinate-preservation gates cannot both be met truth-blind,
  stop this containment line.
- **User-run test needed:** no. This phase used only bounded synthetic fixtures.

### Codex (GPT-6) — 2026-09-14 — Phase B unresolved-owner candidate retained

- **Source/mechanism:** `ProgressiveSequencePool` gained an additive unresolved
  containment tuple while `active_derived` retained its primary-layout meaning.
  Coordinate-compatible children are removed from topology only after one
  provisional output coordinate is established. Their raw placements are stored
  once in unresolved records and are never attached to primary contigs. Consensus
  excludes them; the evaluator uses primary plus unresolved evidence for coordinate
  safety while computing N50 from primary layouts only.
- **Frozen criteria/roles:** development seeds 199001/199002 at 2x/5x/10x/20x;
  six main and three separate error-only scenarios; exact/damage-aware layout and
  consensus. The unchanged Phase B benefit and per-case zero-regression gates apply.
- **Tests:** mandatory clean, reverse, repeat, nested, damaged, CLI-sidecar and
  51-bp strain-regression controls passed. The full suite passed 332 tests with one
  external-aligner skip under Python 3.12 and `pysam 0.24.0`. Python >=3.13 remains
  unverified.
- **Results:** main gate passed 48/48 independent cases and 192/192 correlated
  mode/stage rows. Error-only passed 24/24 cases and 96/96 rows. Clean 20x primary
  N50 changed 132→2398 and 135→2382 in both modes; evaluated coordinate bases stayed
  2398 and 2382. The outputs retain 272 and 279 unresolved records separately.
  Low-depth damage-aware cases supplied the mandatory damaged-data benefit.
- **Actions:** across main mode/cases, 2,170 coordinate-compatible candidates were
  proposed, 1,750 had a single provisional owner and became unresolved evidence,
  and 420 abstained. The former failing error-only case retains all 51 concordant
  positions in `[2069, 2120)`, with three unresolved records totaling 153 bases.
- **Resources:** baseline/candidate main runs used 32.28/33.66 seconds and
  1,468,644/1,463,332 KiB peak RSS. Error-only used 16.39/16.42 seconds and
  780,456/780,560 KiB. These include the memory-heavy evaluator evidence.
- **Disposition:** retained for Phase C independent validation. Evidence and exact
  source snapshots are under `results/phase-b-unresolved-containment-dev-v1/`.
- **Next action:** predeclare genuinely uninspected validation seeds, freeze the
  exact candidate/evaluator, and run the same main and error-only panels once.
  Failure rejects without replacement seeds or tuning.

### Codex (GPT-6) — 2026-09-14 — Phase C unresolved-owner validation rejected

- **Frozen protocol:** seeds `20261401`/`20261402` were preregistered before
  execution and run once at 2x/5x/10x/20x across the unchanged six-scenario main
  and three-scenario error-only panels. Baseline, candidate, evaluator, commands,
  stdout/stderr, gates and source hashes are preserved under
  `results/phase-c-unresolved-containment-validation-v1/`. Coordinates are 0-based
  and intervals half-open. `realized-coverage.json` records sampled bases and
  realized per-reference depth; 80:20 remains target allocation.
- **Result:** error-only passed 24/24 cases and 96/96 mode/stage rows. Main failed
  four of 192 rows: seed `20261401`, damage-aware 20x `damage_unique` and
  `profile_mismatch_unique`, at layout and consensus. Reference-0 positions 1876
  and 2333 changed from concordant-correct C to C/T discordance. Evaluated coverage
  decreased, discordance and mismatch rate increased, and N50 changed 137→136.
- **Benefit that did not clear safety:** clean 20x primary N50 reproduced at
  135→2381 and 137→2400, with 272/276 unresolved records and unchanged evaluated
  spans. Per-coordinate safety is a hard gate, so averages and clean benefit do
  not rescue the candidate.
- **Causal observation:** the new T observations at both failed coordinates came
  from primary records; unresolved records contributed C. Separately retaining a
  removed child therefore does not isolate the primary survival set from topology
  suppression. This is evidence against the mechanism, not a threshold problem.
- **Resources:** baseline/candidate main runs used 29.87/33.06 seconds and
  1,465,756/1,462,220 KiB peak RSS; error-only used 13.60/15.38 seconds and
  769,864/770,036 KiB. All stderr files were empty.
- **Disposition:** rejected without tuning or replacement seeds. Seeds are marked
  consumed, candidate-specific production code/tests were removed, and the Phase A
  graph source was restored at SHA256 `5bae98bb7dc5fb4fcedc3619907db82860c63be5d3e331c1f762395c094a828f`.
- **Next action:** stop ambiguous-containment suppression. Reduce the primary
  survival-set failure to a tiny seed-independent fixture, then evaluate a local
  path-emission mechanism that leaves containment membership unchanged and
  requires direct raw-molecule support for emitted joins.

### Codex (GPT-6) — 2026-09-15 — four-read survival-set counterexample

- **Hypothesis:** attempt 2 exposes a terminally damaged primary read because
  topology suppression removes its clean damage-correction partner, independently
  of unresolved-evidence ownership.
- **Reduction:** deterministic delta debugging reduced the frozen seed-20261401
  failure from 416 to four literal reads with 0-based half-open origins `[0,184)`,
  `[21,159)`, `[15,194)` and `[29,165)`. The final read has one terminal C→T
  event at raw position 1. Synthetic truth was used only by the reducer/audit.
- **Measured baseline:** two physical edges, two reciprocal paths, one corrected
  overlap base and primary placement groups `{0,2}` and `{1,3}`. Both outputs
  carry reference C at the affected coordinate.
- **Measured rejected candidate:** the 138-bp clean read becomes one
  coordinate-compatible uniquely owned containment. Its incident correction edge
  disappears; the graph retains one physical edge/one path and zero corrected
  overlap bases. The 136-bp damaged read survives alone with T as primary output.
- **Implementation:** added a literal, seed-independent regression test asserting
  the Phase A sequence, placement groups, ambiguity classification, edge counts
  and path counts. Production assembly code remains unchanged.
- **Evidence/tests:** reducer, full literal fixture and frozen traces are under
  `results/phase-c-unresolved-containment-validation-v1/`. The focused fixture
  test passed 1/1. No actual FASTQ was run.
- **Disposition:** containment suppression is closed. This is the second focused
  Phase B mechanism failure, activating the guide's stop-and-review rule.
- **Reviewer question:** open a new diagnostic-only phase for local path emission
  that leaves containment membership and primary-output eligibility unchanged,
  or retain Phase A and move to a different assembly architecture.
- **User-run test needed:** no. Actual FASTQ testing remains after Phase C.

### Codex (GPT-6) — 2026-09-16 — audited CarpeDeam comparator boundary

- **Decision:** the exact post-layout linker is closed and Increment 3 remains
  prohibited. A separately installed CarpeDeam safe-mode run is retained only as
  a benchmark comparator. It is not Anvaya's backbone. The internal string graph
  remains the prototype baseline; no additional graph or consensus heuristic was
  added.
- **Implementation:** `carpedeam_backend.py` stream-validates FASTA/FASTQ, validates
  the damage-profile pair, verifies the documented CLI, executes an argument list,
  and records executable/input/profile/output hashes, exact arguments, logs,
  elapsed time and FASTA statistics. It forces unsafe=false, min length 31, five
  raw-read iterations, ten total iterations, merge identity 0.99 and safe coverage
  five. Existing paths are refused. Missing, empty or malformed successful output
  is recorded as `invalid_output` before failure is raised.
- **Boundary:** no GPL source was copied or vendored. The wrapper invokes an
  external executable and records its binary hash. CarpeDeam FASTA exposes no
  raw-read placements, so no Anvaya provenance, coordinate-safety, strain-retention
  or post-assembly correction claim is authorized.
- **Source identity:** backend
  `a6563597503c1237bc71808f69b2d36d00bb7f77d9df669cbe6efc04ed0a28f2`;
  CLI `0f77d4c83fe73e64891b03070197324976a8295f78b2e848af45439d613432d2`;
  package exports
  `1d2fd7cc0b46636b5c580c0a1542500592ea1fcdbca1480acd35360246bf5c91`.
- **Tests/environment:** 340 tests run under Python 3.13.15 and pysam 0.24.0;
  339 passed and one external-aligner test skipped. The wrapper tests cover input
  scanning, pinned arguments, short-read rejection, backend failure and malformed
  successful output. Evidence is in
  `results/carpedeam-backbone-integration-dev-v1/unit-tests.txt`.
- **Real-binary synthetic smoke:** the installed executable at
  `/root/miniforge3/envs/carpedeam/bin/carpedeam` exposed every required option and
  has SHA256
  `3c8bf53cd7d0c429451b5dcf4ba0d08e9311960765121a755ca729624cba45e8`.
  A deterministic 300-read, 30,000-base synthetic FASTQ completed through the
  wrapper in 0.53 seconds and emitted one 1,196-bp contig. Diagnostics and logs
  are preserved with the unit-test evidence. This verifies integration, not
  biological accuracy.
- **Matched input:** the frozen 50k FASTQ SHA256 is
  `d501d6635c656d495868fc43678dec1cadf90348f6ff5eb56d0932e3b5f1b18b`.
  Its profile hashes are 5p
  `b094d9ffa0c54cd219c3ca9b0fa4580ef5a8b19a9205531d28d672871d9570d4`
  and 3p
  `5b50fd8fc7b1710fa50510385381436c5701094882e28a906e6311256ef9cd99`.
  The frozen internal result has 50,000 reads, 43,478 contigs and elapsed time
  67.30 seconds; its earlier N50 response was 62→65 and failed the 1.10 gate.
- **Current disposition:** implementation and installed executable are ready for
  one user-run actual-data checkpoint. No actual FASTQ was run and no positive
  result is claimed. After execution, apply operational gates first, then an
  identically configured reference evaluation.

### User-run result — 2026-09-16 — matched 50k comparator response

- **Identity/status:** the frozen 50k input and profile hashes matched the
  preregistered values. CarpeDeam safe mode completed with exit status zero using
  binary SHA256
  `3c8bf53cd7d0c429451b5dcf4ba0d08e9311960765121a755ca729624cba45e8`.
- **Operational gates:** 34.75 seconds wall time and 75,044 KiB peak RSS pass the
  120-second/1-GiB limits. Output contained 4,962 contigs from 50,000 reads, ratio
  0.0992, and N50 was 112 versus input-read N50 62, ratio 1.806. All preregistered
  response gates passed.
- **Required caution:** output totals only 510,650 bases, versus 2,636,221 bases
  in the internal 50k assembly. CarpeDeam's documented default emits extended
  contigs only. N50 and contig-count improvements therefore include an output
  selection/denominator difference and are not evidence of better recovery or
  accuracy by themselves.
- **Disposition:** proceed once to the matched reference-backed comparison using
  the same `refs/all.fa`, minimum contig 31, identity 90% and two threads for both
  outputs. Keep CarpeDeam external; do not integrate it as Anvaya's assembly
  engine. Use aligned length, genome fraction, unaligned length, NA50, mismatches,
  indels and misassemblies to decide which internal stage remains deficient.

### User-run result — 2026-09-16 — matched 50k reference evaluation

- **Evaluator:** MetaQUAST used the same 474,607,190-bp `refs/all.fa`, minimum
  contig 31, minimum identity 90% and two threads for both assemblies. It completed
  successfully in 21:14.48 with 3,954,272 KiB peak RSS. Anvaya `before.fasta` and
  `consensus.fasta` are byte-identical, SHA256
  `5e83212b909b7d52db9451a1445ed44f3e7552c82b106491a5dd415d6991ae1e`.
- **Anvaya versus CarpeDeam:** total aligned length 1,155,232 versus 455,125;
  genome fraction 0.396% versus 0.190%; mismatches/100 kbp 749.46 versus 1,066.74;
  indels/100 kbp 1.21 versus 3.74; duplication ratio 1.080 versus 1.016. Both had
  zero reported misassemblies. CarpeDeam led raw N50 112 versus 65, auN 116.5
  versus 71.4 and longest 346 versus 287.
- **Denominator finding:** Anvaya emitted 2,636,221 bases; CarpeDeam emitted
  510,650. Anvaya's 1,471,309 unaligned bases include 1,302,257 bases in 28,536
  contigs shorter than 65 bp; only 169,052 unaligned bases occur in 2,368 contigs
  at least 65 bp. CarpeDeam had 51,889 unaligned bases, including 37,007 below
  65 bp and 14,882 at least 65 bp. Zero misassemblies among these short outputs is
  weak evidence.
- **Like-for-like preread:** selecting Anvaya contigs with at least two distinct
  raw molecule IDs gives 5,125 records, 478,252 bases, N50 101 and longest 287.
  At three molecules it gives 1,067 records, 124,354 bases and N50 125. The
  two-molecule set is the appropriate next comparison with CarpeDeam's 4,962
  extended-only records; no sequences have been changed.
- **Interpretation:** Anvaya is already stronger on recovery and base accuracy in
  the full-output comparison. Its reported N50 is dominated by preserving
  singleton/unresolved fragments, while a smaller genuine contiguity gap remains
  within assembled records. CarpeDeam remains an external comparator, not a
  dependency or source of implementation code.
- **Next action:** perform one immutable output split and matched evaluation before
  changing assembly logic. If it passes, formalize assembled versus unassembled
  outputs. Subsequent algorithm work is restricted to evidence-backed branch
  resolution in `_project_master_edges`; repeated containment, linking and
  correction mechanisms remain closed.
- **Diagnostic implementation:** `experiments/split_assembled_output.py` streams
  FASTA records and grouped placements, counts distinct molecule IDs, writes
  disjoint assembled/unresolved FASTAs and an audit manifest, and never reads or
  transforms placement coordinates. Two focused synthetic tests pass.

### Codex (GPT-6) — 2026-09-16 — assembled-output contract

- **Matched result:** Anvaya assembled-only versus CarpeDeam had aligned length
  408,753 versus 455,125 (89.81%), genome fraction 0.169% versus 0.190%, N50
  101 versus 112, NA50 100 versus 111 and longest 287 versus 346. Anvaya retained
  lower mismatches/100 kbp, 815.65 versus 1,066.74, and lower indels/100 kbp,
  2.45 versus 3.74. Both reported zero misassemblies. The bounded normalization
  gate passes: recovery/contiguity are competitive and accuracy remains better.
- **Partition accounting:** all 43,478 records and 2,636,221 bases occur exactly
  once across 5,125 assembled records/478,252 bases and 38,353 unresolved
  records/2,157,969 bases. Coordinates in the placements input are 0-based and
  half-open but were not read or transformed. The assembled/unresolved SHA256s
  are `ff7ce7460377303d525b68a3e1b0613be2a1e6ec501e53a8ab2f0f1fb7aa82b2`
  and `08fb8237aac9487303fb5697d2f619ead398a0a714bc5fa8059cbde8e8addf63`.
- **CLI contract:** optional paired `--output-assembled` and
  `--output-unresolved` now partition consensus records by at least two distinct
  uniquely placed molecules. Existing all-record outputs remain unchanged.
  Original `unitig_N` identifiers are stable across full and partition FASTAs.
- **Verification:** seven focused tests passed after one caught-and-fixed loop
  parenthesis error. Full suite: 343 run, 342 passed, one external-aligner test
  skipped under Python 3.13.15/pysam 0.24.0. Evidence is
  `results/carpedeam-backbone-integration-dev-v1/unit-tests-after-output-contract.txt`.
- **Next action:** user runs one 50k CLI contract check. Exact hash agreement with
  the validated split closes output normalization; no repeat MetaQUAST is needed.
  Only then begin a bounded raw-thread branch-resolution feasibility audit.

### User-run result — 2026-09-16 — output contract verified

- **Determinism:** all four comparisons were byte-identical. Before and consensus
  retained SHA256
  `5e83212b909b7d52db9451a1445ed44f3e7552c82b106491a5dd415d6991ae1e`;
  assembled retained
  `ff7ce7460377303d525b68a3e1b0613be2a1e6ec501e53a8ab2f0f1fb7aa82b2`;
  unresolved retained
  `08fb8237aac9487303fb5697d2f619ead398a0a714bc5fa8059cbde8e8addf63`.
- **Accounting:** 50,000 input reads produced 43,478 all-record contigs and
  2,636,221 bases, partitioned into 5,125 assembled contigs/478,252 bases and
  38,353 unresolved records/2,157,969 bases. Consensus again changed zero bases.
- **Resources:** external wall time 62.07 seconds and peak RSS 331,148 KiB, within
  the frozen 120-second/1-GiB limits. Application elapsed time was 62.84 seconds.
- **Disposition:** output normalization is complete and requires no repeat
  MetaQUAST. Next work is a diagnostic-only independent raw-thread audit of the
  628 physical edges rejected at 329 ambiguous ends. No production branch rule is
  authorized until measurable transition support exists.

### Codex (GPT-6) — 2026-09-18 — branch-transition feasibility closed

- **Boundary:** one diagnostic-only run used consumed seed 199001, clean unique
  20x and damage unique 10x. Coordinates are 0-based and half-open. Endpoint
  molecules, duplicate placements and molecules supporting competing pairings
  were excluded. Synthetic origins were not used for discovery.
- **Measured opportunity:** clean data had 314 physical branch nodes, 1,661
  possible transitions, 75 with at least one independent molecule and 41 with at
  least two; 41 were uniquely resolved. Damaged data had 139 branch nodes, 514
  transitions, 27 with at least one molecule and 14 with at least two; 14 were
  uniquely resolved.
- **Projection check:** the no-transition projection reproduced the existing N50
  exactly. The optimistic no-node-reuse projection changed clean N50 132 to 152
  (1.15x) and damage N50 145 to 156 (1.08x).
- **Verification:** indexed placement enumeration matched an exhaustive
  orientation/offset oracle on a deterministic fixture; a two-molecule branch
  fixture resolved only its supported transition; instrumentation preserved
  assembly output and diagnostics. Focused tests passed 14/14.
- **Disposition:** closed without production path emission. Despite measurable
  evidence, the clean projection cannot reach the preregistered 1.5x benefit gate.
  Do not tune support, replace seeds or run actual FASTQ for this mechanism.
- **Evidence:** `results/branch-transition-opportunity-dev-v2/` records the frozen
  manifest, source hashes and complete case diagnostics.

### Codex (GPT-6) — 2026-09-18 — existing paired mechanisms rejected

- **Overlapping-mate merge:** a bounded clean seed-199001 20x comparison used 356
  genuine fragment pairs with 75-bp mates and 90–180-bp inserts. The unchanged
  exact string graph assembled unmerged mates to one 2,399-bp contig. The existing
  merger combined 138 pairs, retained 218 as separate mates, and reduced N50 to
  99 across 217 contigs. This fails clean preservation and is closed without a
  damaged-data run or threshold tuning.
- **Existing empirical scaffold evidence:** the preserved TAF016 100k run merged
  25,459 pairs, but only 591 pairs mapped uniquely to progressive contigs. Insert
  calibration was invalid, no pair crossed contigs, and no scaffold link was
  emitted; N50 remained 78. The 500k run likewise recorded zero supported links.
- **Tests:** focused paired merger/scaffolder tests passed 7/7. They verify code
  behavior but do not override the negative opportunity measurements.
- **Disposition:** pre-assembly mate merging and post-contig scaffolding are
  closed for this prototype. Only a bounded graph-level mate-transition audit is
  open; it must use synchronized paired reads as linkage evidence and cannot use
  EMN001 because that input contains already merged single fragments.

### Codex (GPT-6) — 2026-09-18 — paired branch audit ready

- **Implementation:** `paired_branch_opportunity.py` builds the unchanged exact
  damage-string graph from interleaved synchronized mates and audits the reduced
  graph through the existing read-only callback. It does not merge mates, emit a
  candidate assembly or use reference truth.
- **Evidence rules:** a pair must place uniquely with opposite inward orientations
  and jointly span the central node's two 0-based half-open junctions. Endpoint
  molecules, multi-placement pairs and molecules supporting multiple transitions
  abstain. A transition requires two independent molecules; projection never
  reuses a physical node.
- **Resource bound:** the raw index stores eight deterministic 15-mers per read
  orientation and verifies candidates base-for-base. The first empirical audit is
  capped at exactly 25,000 TAF016 pairs/50,000 reads.
- **Verification:** focused paired/branch tests passed 12/12. Full suite passed
  347/348 with one external-aligner skip. No actual paired FASTQ was run by the
  agent.

### User-run result — 2026-09-19 — paired threading rejected

- **Cohort:** synchronized ordered TAF016 prefix, exactly 25,000 pairs/50,000
  reads. Input SHA256 was `828fc4e4dbfb345343b9727fe6940e6c69538271f5df3445e36477ab6c8dbb14`
  for R1 and `a1b03c0891ba6b3d3d1f0b2524ebbe44a7db6e22b2348579258959a62dd1e721`
  for R2. Coordinates were 0-based and half-open.
- **Graph:** 26,417 sequence nodes, 31,625 input physical edges, 19,733 after
  transitive reduction, 2,769 physical branch nodes and 17,499 branch-rejected
  physical edges. Baseline output had 23,465 contigs and N50 76.
- **Opportunity:** 30,114 possible branch transitions; zero with one eligible
  independent mate pair, zero with two pairs and zero uniquely resolved. The
  projection remained N50 76 and its no-transition baseline reproduced N50 76.
- **Zero-result cross-check:** synchronized names/counts passed; the positive
  synthetic fixture detects a supported transition; the empirical run exited
  successfully; and projection equivalence passed. Endpoint molecules, duplicate
  placements, same-orientation mates and cross-transition conflicts remained
  excluded by the frozen rule.
- **Resources:** 44:25.22 wall time, 2,659.27 seconds user CPU and 1,056,460 KiB
  peak RSS. The route fails both opportunity and engineering response.
- **Disposition:** graph-level paired threading is closed without production code,
  larger cohorts, damage-tolerant placement or threshold changes. Evidence is in
  `results/paired-branch-taf016-25000-v1/` and the adjacent run log.

## Implementer findings — Codex (GPT-5) — 2026-09-20 — Phase A

- **Hypothesis:** the clean mate-merging collapse is caused by mixing exact
  containment with extension topology. A containment-free exact suffix-prefix
  graph should recover the attainable clean interval without deleting raw
  evidence. All coordinates below are 0-based and half-open.
- **Implementation:** added `src/anvaya/exact_string_graph.py`, outside the
  production route. It groups identical canonical strings, records every raw
  evidence placement, removes strictly contained strings from topology, discovers
  proper bidirected dovetails with an exact suffix seed, removes shift-equivalent
  transitive paths and emits deterministic maximal non-branching unitigs. Multiple
  maximal placements remain explicit and unassigned. Existing assembly code and
  CLI behavior were not changed.
- **Independent checks:** added `tests/test_exact_string_graph.py`. A test-only
  exhaustive oracle agrees with indexed containment, proper-dovetail discovery
  and transitive reduction across 12 randomized variable-length cases. Tests also
  cover duplicates, reverse orientation, multiple containers, prefix/suffix/
  internal containment, branches, cycles, isolated nodes, input permutation,
  reverse-complement invariance and cumulative 0-based layout shifts.
- **Consumed merge regression:** 138 pairs merged correctly, producing 574 input
  records. The candidate has 144 topology strings, 393 contained groups, 239 raw
  reads with multiple admissible placements, 1,166 candidate/286 reduced directed
  edges and ONE exact 2,399-bp unitig. It equals reference interval `[1, 2400)`;
  reference position 0 has no read coverage. The previous graph produced N50 99
  on the same correctly merged records. This is a regression result, NOT fresh
  validation.
- **Synthetic resources:** 1,000 reads: 0.052 s kernel time, 33,724 KiB peak
  process RSS, three unitigs, longest 11,157 bp. 5,000 reads: 0.374 s, 86,536 KiB,
  five unitigs, longest 24,471 bp. These seeded clean simulations test scaling and
  geometry only; they are not biological accuracy estimates.
- **Verification:** Python 3.13.15. Focused tests RUN 8, PASSED 8, SKIPPED 0.
  Full suite RUN 356, PASSED 355, SKIPPED 1 in 2.957 s. Compile check and
  `git diff --check` passed. No actual FASTQ dataset was run.
- **Source identities:** `exact_string_graph.py`
  `203cd8c4d43177009ad623def2ae8bb309fd9d9fcc4b9879a3a2b9cb8cf325d7`;
  `test_exact_string_graph.py`
  `62fd9b698bca32bc3a8dc68c9fce9bb15f3bb63c2c505b568b51fe24b98a740c`;
  existing `reads.py` `571d21ae8d832f91b144575153eafdecc12089ef4cf88f5ea925c96dc90f29d9`,
  `sequences.py` `b7fcd95a590cc03ee8b462024cbf4a01537f01440abc03a0b38c6aae188bba1d9`,
  and `paired_reads.py` `ab5487bdd282debb1e8c38f771c58358611d2f858427050795a79f1d4e00bd5e`.
- **Literature used:** the Readjoiner paper defines exact suffix-prefix matches
  and irreducible string-graph edges; SGA's official README corroborates overlap
  generation followed by graph assembly. The implementation is independent; no
  source was copied and no external assembler was introduced as a dependency.
- **Disposition:** **PASS Phase A**. The guide permits Phase B: integrate immutable
  raw evidence and existing damage likelihood once, preserving all safety gates.
  No actual-file command is needed yet, and the new kernel must remain outside the
  production route until Phase B's frozen matrix passes.

## Implementer findings — Codex (GPT-6) — 2026-09-22 — Phase B

- **Hypothesis and boundary:** one constrained candidate combined the Phase-A
  exact topology with the existing terminal-damage likelihood and raw consensus.
  It remained outside the CLI/default assembly route. Coordinates are 0-based and
  intervals half-open. No reference truth entered assembly, and no actual FASTQ
  dataset was run.
- **Implementation:** `exact_string_graph.py` now supports correction-bearing
  bidirected edges, correction-aware path spelling/transitive reduction and one
  damage-compatible edge layer using the existing >=3 independent-molecule,
  Q20, supplied-profile gate. `evidence_pool` composes immutable raw placements
  through canonical unitig layouts. A raw read is assigned only if all proposed
  compositions collapse to one global `(unitig, orientation, offset)`; otherwise
  it remains unresolved. Contained reads can therefore contribute without an
  arbitrary parent choice or duplicate molecule vote. The production route is
  unchanged.
- **Invariant tests:** the damage candidate reconnects a terminal C>T/G>A chain,
  produces the same canonical sequence under input permutation and whole-input
  reverse complementation, and projects each molecule once. Separate tests cover
  globally ambiguous branch endpoints, unique contained-read offset composition
  and a transitive path whose corrected spelling differs from the direct edge.
  The historical four-read baseline test was not changed.
- **Runner/gate changes:** `variable_length_validation.py` gained explicit
  `--implementation baseline|candidate`, fixed development seeds 199001/199002,
  default depths 2x/5x/10x/20x and assembly-only timing. Both outputs use the same
  generator/evaluator source and input hashes. `paired_coordinate_gate.py` now
  rejects candidate assembly time above 2x its matched baseline. Scenario names,
  profile settings and evaluator logic were unchanged.
- **Frozen main matrix:** RUN 48 cases/192 correlated mode-stage rows. The benefit
  criteria passed: clean unique 20x N50 was 132->2,398 for seed 199001 and
  135->2,382 for seed 199002 in both modes and stages; each seed also had a
  damage-aware benefit case. There were zero N50 regressions and zero >2x timing
  failures. Safety **FAILED in 140/192 rows**: unresolved burden increased in
  132, mismatch fraction increased in 60, evaluated span fell in 38, the aggregate
  discordant-reference burden increased in 8, mosaic burden increased in 6 and
  incompatible-sequence burden increased in 4. Coordinate transitions included
  9,318 previously correct losses, 558 new discordant events and 26 newly
  evaluated incorrect events; rows/modes are correlated and are not independent
  biological replicates.
- **Consumed failure checks:** in seed 199001, damage-error strains 80:20 at 2x,
  every position in reference-1 interval `[2069,2120)` changed from
  concordant-correct to not uniquely evaluated in exact/damage-aware layout and
  consensus (51/51 positions in each of four correlated rows). In seed 199002,
  profile-mismatch unique 20x, reference-0 position 1876 stayed
  concordant-correct in all four rows; position 2333 stayed correct in exact mode
  but became not uniquely evaluated in damage-aware layout and consensus. Evidence
  preservation therefore failed even though the earlier C/T-discordance symptom
  changed to abstention.
- **Frozen error-only controls:** RUN 24 cases/96 correlated rows; safety **FAILED
  in 80/96 rows**. Unresolved burden increased in 72 rows, mismatch fraction in
  44, evaluated span fell in 28, aggregate discordant-reference burden increased
  in 16 and incompatible-sequence burden increased in 4. There were 4,452
  previously correct losses, 788 new discordant events and 32 newly evaluated
  incorrect events. There were zero N50 and zero >2x timing regressions.
- **Resources:** main baseline 30.51 s wall/1,468,904 KiB peak process RSS;
  candidate 14.70 s/1,248,092 KiB. Error-only baseline 15.31 s/781,508 KiB;
  candidate 8.11 s/697,204 KiB. These process RSS values include evaluator state;
  the matched assembly timing gate is recorded per case and had zero failures.
- **Verification:** focused Phase-B/history tests RUN 38, PASSED 38, SKIPPED 0
  before the final invariant extension. Final unittest suite RUN 361, PASSED 360,
  SKIPPED 1 in 2.649 s. Syntax and whitespace checks passed. An initial unscoped
  `pytest -q` command was invalid because it collected archived test snapshots
  under `results/` and used system Python without `pysam`; it ended with 28
  collection errors and was replaced by the repository-scoped Conda unittest
  command above.
- **Source identities:** `exact_string_graph.py`
  `baa769f420ebec37079fbdb6872e105094943ceef3e3fd6c55196d465381a02e`;
  `test_exact_string_graph.py`
  `308cfbe64dc282c2bfea62043d0e6dd4796b22c57cc177f1fd2d9ae8d9841a84`;
  `variable_length_validation.py`
  `7e70b8fbb2ab4d5cd74b5a2bec03b6b7a6b21b6931af1e0e4f2e3ca201678d97`;
  `paired_coordinate_gate.py`
  `b90aef45aa24dfcee09494649772112f140e12d4a45bd012caa5be4d1cb314b5`.
  Main baseline/candidate/gate summary SHA256 values are respectively
  `bfcbdf8fa98f0f6b1f0ec8a56557fba6a49768d983d55821568fbd436b77686d`,
  `1e10e1feab85c18f6d3c0d9688f2e235bb836c3ac54ea52eabb7a863838832ae`
  and `c8fc8f50b44cf5f5d487ccf05aa1b00b98be97557c870babcaf818894bfa96d9`.
  Error-only baseline/candidate/gate values are
  `f2264b131d0d9ac6e1738cb9988815613e13c0ceb33e29100e24a36cd5c543fa`,
  `1c5981325fb68e1741fca7b30cfa4b52d1acc5d274d6bb076aefddc2cc1348b4` and
  `fd47d8d51923d398e60a63d670b3a477b9bc2044878606284027e6dc1c4ec8f9`.
- **Literature cross-check:** CarpeDeam applies supplied damage likelihood only
  to unextended fragments and reports a sensitivity/misassembly/duplication
  tradeoff even for its safe mode. That supports retaining raw-end coordinates
  and rejecting long but unsafe unitigs; it does not justify treating damage
  compatibility as strain identity. PEAR and AdapterRemoval were consulted for
  mate overlap/orientation context; no code or external assembler was copied.
- **Disposition:** **REJECTED Phase B**. The candidate passed benefit and runtime
  gates but failed mandatory coordinate safety in both main and error-only
  controls. The exact unitig core and candidate remain isolated evidence; the
  safe baseline remains the production route. Under the active guide, Phase C,
  an actual-file command, threshold tuning and a second architecture are not
  permitted. The next action requires a new reviewer decision.
