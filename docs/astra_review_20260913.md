# Astra implementation review — 2026-09-13

## Verdict and scope

The implementation is a useful bounded experimental assembler. It is not yet a
validated damaged-data contiguity improvement. The clean coverage ladder reveals
a real layout problem, but the latest maximum-overlap proposal precedes a causal
separation of graph-input, reduction and path-selection failures.

Reviewed current damage graph, shared projection/bidirected edges, candidate
index, phase adapter, raw consensus, variable-length driver, paired gate, cohort
proxy, relevant tests, guide/findings, saved smoke/proxy/ladder evidence and primary
literature. Older progressive modes were checked through shared interfaces and
regression tests, not exhaustively re-proved. Graphify lacks newer modules; source
is authoritative. The dirty tree cannot establish individual authorship. Recent
findings explicitly include GPT-6 continuations as well as Sol's earlier work.

Evidence: `results/astra-review-20260913/audit.py`, `audit.json`, and pre-edit
guide/findings snapshots. Only synthetic reads were assembled. All coordinates
are 0-based, half-open.

## Prioritized findings

| Priority | Source and finding | Required response |
|---|---|---|
| P1 | `overlap_graph.py::_project_master_edges` checks two-hop witnesses only | Tiny strict-dovetail probe demonstrates a missed redundant edge blocking correct assembly; quantify its contribution before generalizing. |
| P1 | `damage_string_graph.py::assemble` admits boundary containments as extensions | These can have zero progress in one orientation. Separate relation semantics; blanket child deletion is not a safe fix. |
| P1 | `paired_coordinate_gate.py::_case_index` keys only `(seed, scenario)` | The new multi-depth ladder produces duplicate keys. Make comparison depth-aware before testing another candidate. |
| P1 | 50k proxy interpretation extrapolates from 578 eligible reads | The 178-bp component maximum applies to 1.156% of reads, not the dataset. Excluded reads can bridge components. |
| P1 | Next plan prescribes maximum-overlap reciprocal selection | Length alone does not resolve strain/repeat ambiguity. Suspend pending a measured branch classification. |
| P2 | 50k findings copied incorrect reduction counts | Correct totals: 4609−213=4396; 4396−628=3768 reciprocal (81.75%). |
| P2 | Per-node raw/node ID list allocation in `assemble` | Avoidable O(nodes × reads) work; hoist in a separate profiled change. No speedup measured yet. |
| P2 | Raw candidate search runs before knowing whether a mismatch needs evidence | Consider lazy evaluation; require sequence/provenance equivalence and document changed counters. |
| P2 | `phase_graph.py` repeatedly scans `min(unseen)` | Potential quadratic component traversal; profile before changing. |
| P2 | Validation exclusions omit development seeds 199001/199002 | Keep a unified seed-role ledger so inspected seeds cannot be labelled fresh. |
| P2 | Strain fixtures share RNG streams across references | Depth changes alter minor-strain sampling: these are not nested paired depth experiments. Version fixture changes. |
| P2 | Driver retains complete case evidence and only checks `coverage <= 0` | Stream completed evidence if memory blocks iteration; reject nonfinite depths and preflight the assembler read cap. |
| P2 | Cohort filter checks names, not sequence identity, and retains all seen BAM names | ID equality is not content provenance; memory is not exclusively cohort-sized. |
| P2 | Proxy hashes repo paths rather than necessarily imported modules | Record actual import paths/hashes, inputs and environment. Current `PYTHONPATH=src` mitigates the recorded-run risk. |
| P3 | `fastq_names` indexes `split()[0]` before validating an empty identifier | Empty `@` header raises IndexError. Add a malformed-input regression when touching the parser. |

These are source behaviors or labelled performance risks, not a claim that every
one drives real-data fragmentation. Do not delete older experiments or combine
these fixes into an unrelated refactor.

## Reproduced results and interpretation

Four 100-bp reads start at reference positions 0, 20, 40, 60. Strict edges A→B→C→D
and A→D spell the same sequence, but A→D has no two-hop witness. Current reduction
removes zero edges; branch filtering rejects three physical edges and outputs
120/100/100 bp. Removing only the redundant edge yields one 160-bp contig with
the SAME degree-one policy. This proves a reduction limitation, not its share of
the ladder losses or the safety of arbitrary longer-path reduction.

Clean seed 199001 at depths 2/5/10/20 presents 46/324/1346/6072 directed edges to
projection. Of these, 0/4/22/124 violate
`0 < overlap < min(source_length, target_length)`. The audit only counts these
relations; it does not delete them or claim that deleting them improves assembly.
Oriented counts must not be mixed with physical-edge counts.

The saved 50k before/consensus FASTAs share SHA256
`5e83212b909b7d52db9451a1445ed44f3e7552c82b106491a5dd415d6991ae1e`.
The 373 overlap-base corrections occur before that comparison. Identical FASTA
means the additional consensus pass changed nothing, not that the damage-edge
mechanism was inactive. No repeated polishing or identical-sequence QUAST run is
warranted.

The 50k output is 43,478 contigs, N50 65 versus input N50 62, longest 287 bp.
The N50 response gate failed. External time 62.65 seconds is less than application
time 67.30 seconds: inconsistent for the same enclosing wall interval. Both are
under budget, but retain the discrepancy instead of fitting precise scaling.
Missing newer instrumentation and later source modification times alone do not
prove an older run used a stale installed package.

The selected proxy's 65/70 directed overlaps are correlated and strongly biased
toward uniquely mapped reads. They do not measure full-index recall or establish
taxonomic ordering, a global coverage ceiling, or a maximum possible benefit from
candidate changes. High accepted-edge survival similarly does not label rejected
true overlaps. The failed response gate still justifies stopping blind scaling.

The clean ladder's mean N50 decline 237→180→149→133.5 bp is useful evidence.
Two seeds and correlated mode/stage rows are not independent replicates. Means
cannot pass per-case gates. Its 1,467,144-KiB process RSS includes full evaluator
evidence and cannot be assigned solely to the assembler. Uniform lengths 31–200,
Q10-labelled substitution errors and two strains are simplified fixtures.

## Primary literature cross-check

- [Myers (2005)](https://academic.oup.com/bioinformatics/article/21/suppl_2/ii79/227189),
  **10.1093/bioinformatics/bti1114**: formal overlap/string graphs and transitive
  reduction. It does not endorse the current sparse two-hop filter as complete
  or longest-edge selection as strain-safe.
- [Simpson and Durbin (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2881401/):
  efficient exact string-graph construction. Check containment/input assumptions
  before transferring its invariants to retained ambiguous ancient fragments.
- [Kamath et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610600/): contained-read
  removal can affect connectivity. This motivates auditing evidence loss, not
  declaring all containment removal safe or all multiple parents distinct strains.
- [Gargammel](https://pmc.ncbi.nlm.nih.gov/articles/PMC5408798/),
  **10.1093/bioinformatics/btw670**: fragment, damage and sequencing simulation.
  The prior **btx013 DOI is incorrect**. The repo's custom uniform-fragment fixture
  is not a Gargammel implementation or validated empirical model.
- [PyDamage](https://pubmed.ncbi.nlm.nih.gov/34395085/), **10.7717/peerj.11845**:
  damage assessment of assemblies. The prior **10.1093/nargab/lqab096 is incorrect**.
  Authentication evidence does not validate an extension heuristic.
- [metaSPAdes](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411777/),
  **10.1101/gr.213959.116**: coverage variation and strain/repeat ambiguity are
  relevant concepts, not a direct specification for this overlap assembler.

Graph compatibility uses `rate > 0`, so rate scaling with an unchanged support
mask cannot test graph calibration. The half-strength fixture can affect numerical
consensus, but does not validate a quantitative edge-damage model. Keep this claim
boundary without launching another scoring redesign during topology work.

## Verification and disposition

Full unit discovery with isolated pysam 0.24.0 ran **321 tests: 320 passed, one
external-aligner integration test skipped**. WSL Python was 3.12; package metadata
requires >=3.13. Supported-version verification remains open. The new standalone
audit passed its reduction and chain-control assertions. No production code was
changed and no real FASTQ benchmark was run.

Retain provenance, indexed containment, useful diagnostics, and rejected-candidate
evidence. Next: bounded graph-invariant trace, one causally supported topology
repair, per-case safety/benefit gates, fresh validation, then a user-run matched
50k comparison. The active guide specifies branches and stop conditions. Biological
improvement remains an open milestone, not a promised outcome.
