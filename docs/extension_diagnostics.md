# Extension bottleneck experiment

The current task is to explain fragmentation before relaxing assembly thresholds.
Consensus and guard tuning are paused. Candidate selection and extension behaviour
are unchanged by this instrumentation.

The progressive raw phase emits `overlap_progressive_round_*` counters. Each
attempted round receives one outcome: `no_candidates`, `quality_filtered`,
`no_dovetail`, `ranking_rejected`, `boundary_support_rejected`,
`reciprocal_rejected`, `no_growth`, or `grew`. These count round outcomes, not unique
reads or exclusively terminal stopping events: internal correction may allow a
later round even when a round does not grow. `no_candidates` is after validated
candidate discovery, not proof that no true overlaps exist. Ranking rejection
combines the ranking function's reasons; it is not exclusively ambiguity.

`overlap_progressive_recruit_search_*` exposes the existing offset-level candidate
counters for subsequent recruitment searches in this phase. It excludes initial
cluster discovery and reciprocal reverse searches. One read can generate multiple
offsets and be evaluated in multiple rounds; do not interpret these counts as read
recall. Existing later progressive-iteration rejection counters remain separate.

For truth-based candidate recall, run:

```bash
PYTHONPATH=src python experiments/extension_candidate_recall.py --output NEW_DIR
```

This uses 640 small synthetic cases with known extending read/strand/offset triples,
50/100 bp reads, 30/40 bp overlaps, damage 0/.4 and independent sequencing-error
probability 0/.01, both orientations and 20 seeds. Each target has one true partner
and three random decoys. Positions/offsets are **0-based**, overlaps half-open.
The trace identifies the first failed stage at the true placement: absent anchor
vote, insufficient anchor support, short overlap, DNA identity, RY identity or
molecule selection. Trace collection is optional and used only by the benchmark.

Settings match the current comparison's k=15, eight read anchors, occurrence cap
100, one required anchor match, minimum overlap 30, DNA identity .90 and RY identity
.99. These fixtures do not benchmark repeated genomes, alternative strains,
gapped overlaps, clustering/ownership, ranking or full assembly. Random decoys are
easy negatives; zero false selections is not a precision estimate for metagenomes.
Damage and error simulation here is a small mechanism fixture, not a replacement
for external biological simulation. Configurations reuse seeds, so cases are not
independent replicates.

Verified results are in `results/extension-candidate-recall-v2`: 463/640 true
placements selected, 78 without an anchor vote at the true offset and 99 rejected
by RY identity. No other placements were selected. With no sequencing errors,
297/320 were selected; with 1% sequencing errors, 166/320 were selected. The 99 RY
rejections occurred in the sequencing-error cases. One R/Y-changing mismatch in a
30–40 bp overlap fails .99 identity. This identifies a candidate bottleneck under
the synthetic assumptions, not authorization to relax it on real data.

The v1 report had an incorrect hard-coded case total (1280); its raw stage counts
were 640. The benchmark now counts emitted cases and asserts stage accounting.
Use v2. The full unit suite passed 215 tests including output-neutral candidate
tracing, true-offset RY rejection and a contained-read no-dovetail round.

The subsequent guarded run in `results/extension-diagnostics-100k` retained the
previous FASTA. Recruitment reported 11,945 offsets: 2,672 short, 2,183 DNA
rejections, 2,951 RY rejections and 4,139 accepted alignments. Later iterations
accepted 1,249 candidates but added no bases. These counts motivate testing
quality-aware recruitment, rather than further consensus-guard tuning.

Experimental recruitment allowance (`--progressive-low-quality-ry-rescue`):
when the RY identity gate fails, allow exactly one R/Y-changing mismatch only
when the candidate's original raw base has Q <= 15. Missing qualities, ambiguous
overlap bases and multiple RY mismatches cannot use the allowance. DNA identity,
anchors, minimum overlap, support and ranking still apply. Initial clustering
and reciprocal searches remain strict. The option applies to subsequent raw
cluster rounds and main derived-center raw recruitment, including priority
searches; it does not enable rescue in the optional support-three pass.

This follows CarpeDeam's separation of DNA and RY matching and its use of
multi-read extension support ([methods](https://link.springer.com/article/10.1186/s13059-025-03839-5)),
with base-quality evidence motivated by
[Li's likelihood framework](https://doi.org/10.1093/bioinformatics/btr509).
CarpeDeam uses RY identity .999; Anvaya currently uses .99. The single-error Q15
cutoff is our provisional heuristic, not a threshold established by either
paper. It is not a likelihood ratio and cannot distinguish a low-quality real
strain variant from an error. It cannot recover overlaps without anchor votes.

New `ry_rescued_offsets` and `ry_rescued_alignments` search counters count gate
passes and unique molecule selections, respectively, not successful extensions
or unique reads across searches. Priority-cache counts precede ownership
filtering. Assess actual growth using added bases and aligned assembly metrics.

Run `experiments/raw_consensus_comparison.sh` with mode `ry-rescue` into a new
directory; this also enables the frozen linked guard. Compare against mode
`guarded`, using identical input and profiles. Because layout can change, remap
both assemblies for evaluation; the old fixed-placement substitution audit is
not suitable. Retain the option only if aligned recovery improves without a
material increase in mismatches or misassemblies. N50 alone is insufficient.

Synthetic unit controls include both strands (raw coordinates 0-based), the
Q15/Q16 boundary, multiple errors, missing qualities, ambiguous bases, unchanged
DNA/length/anchor gates, and a correct 20-base extension needing two rescued
reads in addition to three exact reads. They also explicitly demonstrate that a
low-quality biological difference can pass; these are mechanism tests, not an
estimate of metagenomic precision.
