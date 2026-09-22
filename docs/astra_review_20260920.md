# Reviewer audit — Astra, 2026-09-20

## Decision

Anvaya is an executable research prototype with useful provenance and conservative
consensus machinery. It is not yet a demonstrated competitive damage-aware
assembler. The recent experiments produced useful failures, but several guide
interpretations overclaimed what those failures establish. Shared workspace
history does not identify every edit as Sol's; earlier reviewer directions also
contributed to the loop.

The next investment is a conventional, independently checked exact string-graph
unitig core for variable-length reads, followed by one raw-evidence integration.
This is a NEW reviewer-authorized architecture decision, not permission to revive
the rejected containment-suppression patches. See the
[active implementation contract](sol_implementation_guide.md).

## Evidence and scope

Reviewed the active guide/findings, graph construction and projection, raw-read
consensus, paired input/merging, CLI/output contracts, opportunity evaluators,
recent synthetic/real result artifacts and their interpretation. This is not a
proof that every line is bug-free. Graphify's index was stale; current files and
executable probes were used to verify relevant behavior.

New evidence: `experiments/review_20260920.py`, with source hashes and outputs in
`results/reviewer-20260920/probes.json`. These are synthetic reviewer probes, not
independent validation. The merge probe reuses consumed seed 199001 and placement
seed 199002. Coordinates throughout this review are **0-based, half-open**;
single failing positions are 0-based. No actual FASTQ assembly was rerun.

Full unit suite: Python 3.13.15, 348 tests run, 347 passed, one skipped, 4.132 s.
Log: `results/reviewer-20260920/unit-tests.log`. Passing these tests does not
establish assembly benefit; the missing cases below explain that distinction.

## Findings, in priority order

### P1: correct merging was blamed for a graph failure

The consumed clean fixture contains 356 pairs. The merger emits 138 merged
sequences; every one equals its known original fragment. Nevertheless:

| Exact assembly input | N50 | Contigs | Ambiguous containments | Full-containment directed edges after reduction |
|---|---:|---:|---:|---:|
| Unmerged mates | 2,399 | 1 | 0 | 0 |
| Correct merges plus unmerged mates | 99 | 217 | 265 | 64 |

`damage_string_graph.py:232` classifies right extension without left extension
as a dovetail, including some full-read overlaps. At line 255, only children with
exactly one candidate container leave topology. Multiple containers can describe
one locus covered by several reads; their count alone is not repeat evidence.
`overlap_graph.py:318` then projects a mixed topology through degree rules.

The reproduction establishes correct merged sequence and severe downstream
fragmentation. It does NOT establish that deleting those 64 edges or suppressing
all contained nodes is sufficient or safe. Earlier strain and damage failures
explicitly disprove that shortcut. Separate geometric sequence redundancy from
raw evidence ownership, and test the exact graph against exhaustive overlaps.

### P1: opportunity projections are not upper bounds

`branch_transition_opportunity.py:131` implements `_projected_n50` using greedy
node-disjoint path selection. It does not solve an optimal assembly problem.
Selected transitions are inserted in one canonical orientation without symmetric
mirror insertion; traversal and first-successor choices need invariance checks.
Matching the baseline N50 alone does not establish identical baseline paths,
sequences or molecule ownership. Physical graph nodes are not independent raw
molecules. These outputs cannot certify that no alternative graph method can
reach the benefit gate. Preserve them as configuration-specific projections.

### P1: paired conflict detection rejects compatible uses

`paired_branch_opportunity.py:116` rejects any molecule associated with more than
one transition key globally. A molecule can consistently span multiple nearby
branch locations. A synthetic predicate probe gives one resolved transition with
one branch, but rejects both molecules when a second consistent branch location
is included. This is an abstract graph test, not a real-data rescue result.

Crucially, the real TAF016 audit reported ZERO conflicting molecules, so this
bug does **not** explain its measured zero-support result. That result remains a
null under exact full-mate placements inside short three-node paths, stringent
orientation/start rules and endpoint-group exclusion. It is not a biological
connectivity ceiling. Equal-start overlapping pairs are excluded by the current
strict start ordering. The audit lacks a rejection funnel to attribute its null.

An earlier scaffold calibration failure also returned before measuring later
cross-contig support. Default zero counters after that return are not evidence
that no cross-contig placements exist. Do not use failed calibration to justify
an impossibility claim.

### P2: the 44-minute paired audit has avoidable global scans

The real audit examined 30,114 transitions for 25,000 pairs / 50,000 reads.
It loops over all reads to canonicalize them and over all pairs per transition:
approximately 1.506 billion read visits plus 753 million pair visits. The N50
projector also nests edge scans. This is static complexity evidence, not a
profiler attribution. User timing: 44:25.22 elapsed, 99% CPU, 1,056,460 KiB RSS.

If this experiment is ever reopened, build read/molecule and placement indexes
once and query local candidates. Scope conflicts to mutually incompatible
placements, retain mirrors, validate pair-name uniqueness and mate roles, hash
large inputs incrementally, include all imported source hashes, and add bounded
runtime plus measured rejection stages. Currently FASTA acceptance/name handling
also falls short of the stated paired-FASTQ contract. Do not spend the next phase
optimizing a nonessential audit while the assembler core remains broken.

### P2: prototype contract is ahead of algorithmic evidence

`damage-assemble` has useful assembled/unresolved partitioning, but output paths
can still overwrite existing files. Before promotion, validate path collisions
and existing destinations before assembly and stage outputs before publication.
The public command does not yet provide a validated complete paired assembly
workflow. Do not advertise general paired-end support from an experiment alone.

Containment correction dependencies and raw-end coordinates are load-bearing.
A sidecar with missing primary evidence is not preservation. Retain the original
four-read regression and every earlier coordinate safety gate.

## What worked and where we stand

- The raw provenance/consensus work, coordinate regressions, restored failed
  candidates, output partition, source manifests and matched comparator are
  useful progress. They expose failures instead of concealing them.
- Current 50k output partition: 5,125 assembled contigs / 478,252 bp / N50 101;
  38,353 unresolved records. All-record N50 65 mixes two output classes and must
  not be compared with the competitor's assembled-only N50.
- Matched CarpeDeam assembled output: 4,962 contigs / 510,650 bp / N50 112.
  Anvaya aligned length 408,753 versus 455,125; mismatches per 100 kbp 815.65
  versus 1,066.74. Anvaya is more conservative in this comparison but recovers
  less aligned sequence. Both report zero misassemblies on short contigs; this
  is weak safety evidence. Reference mapping is a proxy, not read-origin truth.
- Prior consensus experiments showed promising reference agreement changes, but
  ambiguous mappings limit attribution. They do not prove that the present
  default graph produces a damage-specific contiguity gain.

We are one major graph-correctness milestone plus integrated validation away
from a defensible minimal prototype. Competitive contiguity remains uncertain;
neither a release date nor a guaranteed N50 gain follows from current evidence.

## Literature and repository checks

1. [Readjoiner, original paper](https://link.springer.com/article/10.1186/1471-2105-13-82)
   describes conventional suffix-prefix overlap/string-graph construction and
   contained-read treatment. This motivates a containment-free exact topology,
   not discarding damaged reads' evidence. The
   [GenomeTools repository](https://github.com/genometools/genometools) is the
   official implementation context.
2. [SGA's official assembly README](https://github.com/jts/sga/blob/master/src/README)
   documents overlap generation, irreducible edges and assembly. Use it to check
   graph semantics and optionally as a small offline oracle. No implementation
   body or pinned commit was audited here; Sol must inspect/pin specific code
   before relying on implementation details. SGA is GPLv3; conceptual independent
   implementation does not authorize copying code into Anvaya.
3. [CarpeDeam paper](https://link.springer.com/article/10.1186/s13059-025-03839-5)
   and [official repository](https://github.com/LouisPwr/CarpeDeam) support staged
   extension and attention to terminal damage. Safe mode trades contiguity for
   caution; it is not a universal safety guarantee. Remains a comparator only.
4. [PEAR paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3933873/) and
   [AdapterRemoval v2 documentation](https://adapterremoval.readthedocs.io/en/2.3.x/manpage.html)
   motivate validated overlap/quality handling. Correct mate merging should not
   be rejected because a downstream graph mishandles varying lengths. Raw mate
   ends and one-molecule identity must survive merging. Do not interpret merged
   consensus qualities as two independent molecular observations.

The useful change of perspective is to make geometry trustworthy first and keep
damage evidence as a separate constrained layer. None of these sources proves
that a generic exact graph alone solves repeat/strain ambiguity in short aDNA.

## Changes made in this review

Added the reproducible synthetic audit, replaced the accumulating active guide,
archived the previous guide, and separated reviewer corrections from implementer
findings. No production algorithm was changed, staged or committed. Next: Sol
implements Phase A of the active guide; no additional real FASTQ run yet.
