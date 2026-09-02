# Benchmark plan

## Development sequence

1. Toy sequence with an exactly known graph.
2. Freeze the versioned P0 ladder in `experiments/benchmark_ladder.json`.
3. One simulated bacterial genome with clean, damage, error, and combined
   controls.
4. Two related strains at controlled abundances, plus contamination.
5. A small simulated bacterial and archaeal community.
6. Selected empirical ancient metagenomes.

Each simulation will separate sequencing error, postmortem damage, and genuine variation before combining them.

Controlled damage tests report both introduced events and damage-derived k-mer edges that survive graph filtering. Sensitivity is calculated against retained ground-truth edges so events removed by `min_count` are not incorrectly counted as detector failures. Empirical non-UDG/UDG comparisons use equal read counts, identical graph parameters, normalized topology rates, and replicate-level statistics.

## Frozen P0 ladder

The current acceptance baseline is five deterministic seeds over six scenarios:
`clean`, `damage`, `sequencing_error`, `damage_error`, `rare_strain`, and
`contamination`. It uses a 5,000 bp random reference, 100 bp reads, 20× major
coverage, `k=31`, and `min_count=2`. The runner is:

```bash
PYTHONPATH=src:tests python3 experiments/18_validation_ladder.py \
  --output-dir experiments/validation/generated/ladder
```

Every run writes a TSV/JSON summary and a manifest containing the configuration
checksum, Git state, runtime, and output checksums. The baseline is frozen until
an explicit benchmark-plan change; new assembly methods must be compared with
the same seeds and scenarios. `--dry-run` validates the matrix without running
assembly.

Primary acceptance metrics are N50, largest contig, total bases, exact
reference recovery, false contigs, rare-strain recovery, and damage-event
counts. Runtime and peak memory are recorded as operational diagnostics, not as
substitutes for accuracy. The synthetic ladder supplies exact generating truth;
the authentic public library supplies realism and non-regression only.

## Comparisons

- Anvaya with damage decisions disabled;
- Anvaya with internally inferred damage evidence;
- MEGAHIT;
- metaSPAdes, including `--only-assembler` where appropriate;
- CarpeDeam safe and unsafe;
- PenguiN if practical.

The damage-disabled mode is an experimental control, not a separate assembler.

## Primary outcomes

- misassemblies and mismatches;
- genome and reference-path recovery;
- correctly aligned contig length;
- genuine-variant retention;
- damage-derived path removal;
- low-abundance strain preservation;
- gene and protein recovery.

Runtime and memory will be recorded but will not guide the first research decision.

## Go/no-go criterion

An assembly change is accepted only if it improves continuity on the frozen
truth-labelled ladder without introducing false contigs or unacceptable loss of
reference or rare-strain recovery. Public-data N50 changes alone are not a
go/no-go signal. Optimization begins only after the same change passes the
accuracy-recovery trade-off on held-out seeds.

## P1/P2 diagnostic baseline

P1 attributes retained compacted-graph dead ends by rescanning only their
`k-1` contexts in the raw reads. P2 audits low-quality substitutions against
the full solid-k-mer spectrum. Both stages are report-only and therefore have
no expected N50 effect by themselves.

Across the frozen five-seed ladder, uniquely filtered continuations explained
64 of 97 dead ends in damage runs, 101 of 135 in combined damage/error runs,
42 of 48 in sequencing-error runs, 146 of 203 in rare-strain runs, and 303 of
339 in contamination runs. P2 found 1,677 candidates in combined runs and
1,868 in sequencing-error runs, but none in damage-only runs because simulated
damage bases retained high quality. These figures identify candidate pools;
they do not establish that restoring an edge or changing a base is correct.

On the unlabeled `SRR32866683` fixed-`k=31`, `min_count=2` graph, P1 attributed
566,116 dead ends: 85,190 read boundaries, no coverage gaps, 254,761 uniquely
filtered continuations, 226,159 filtered conflicts, and six contexts retained
elsewhere. Only 8,299 dead ends were classified as uniquely filtered with
predominantly low-quality support, while 20,902 were terminal-only. This
supports investigating the uniquely filtered pool but shows that Phred-only
correction addresses a small fraction.

A separate evenly spaced 20,000-read P2 audit against the full 1,409,072-read
spectrum reported 2,926 `would_correct`, 50 protected damage-compatible
reversals, 58 ambiguous candidates, and 10,863 incomplete rescues. The later P1
stage in that combined invocation failed on an ambiguous-base edge case that
was subsequently fixed, so the correction TSV is a valid stage artifact but
the invocation is not treated as a successful end-to-end run.

## Filtered-edge oracle decision

The frozen ladder compares baseline assembly with two truth-only upper bounds:
`oracle_unique` rescues only one-way, truth-valid, globally sub-threshold
continuations at retained dead ends; `oracle_all_truth` rescues every observed
truth-valid k-mer below `min_count`. Both rebuild the graph without removing
existing damage or error edges.

The completed 30-run `oracle_unique` comparison rescued 355 k-mers. Mean N50
changed by +0.2 bp for damage, 0 bp for sequencing error, 0 bp for combined
damage/error, -0.2 bp for rare strains, and +12.2 bp for contamination. No
false-contig count changed, but exact primary-reference recovery declined by
1,230 bp in aggregate, including a 1,103 bp loss after two truth-valid edges
connected into an unsafe sequencing-error path.

The singleton dead-end rescue hypothesis therefore fails the go/no-go gate and
must not become a production assembly policy. Damage-aware read consensus was
then evaluated directly: it produced strong synthetic continuity gains, but
lost 8.98% rare-strain recovery and left public `SRR32866683` N50 unchanged at
38 bp. It remains opt-in research infrastructure rather than a production
policy.

The next continuity backend is iterative whole-fragment cluster-and-extend
assembly with R/Y-aware overlap validation and safe extension consensus. It
must use this same frozen ladder and additionally report overlap-cluster purity,
accepted extension support, chimeric-contig counts, rare-strain retention,
runtime, and peak memory. Direct CarpeDeam, MEGAHIT, and metaSPAdes comparisons
must use identical prepared inputs; empirical N50 without truth remains a
resource and realism observation, not a go decision.

## Direct-overlap P0 checkpoint

The separate overlap backend crossed the DBG continuity baseline on the
CarpeDeam EMN001 100k and 500k subsets. On 500k fragments, direct overlap
produced N50 109 versus 58 for the DBG, completed in 34.82 versus 54.97 seconds,
and used 1.00 versus 2.85 GB peak RSS. Reusing reads in a second layout phase
increased recovered output from 3.52 to 4.32 Mb but reduced N50 to 104. No
contig reached 1,000 bp, so this is proof that the replacement architecture is
viable, not evidence that its layout problem is solved.

Frozen-layout terminal-damage polishing is accepted as an opt-in accuracy
stage, not as a continuity stage. In the exact paired 500k comparison it kept
45,545 contigs, total length, N50, and largest-contig length unchanged; retained
zero MetaQUAST misassemblies; increased aligned length by 6,172 bp and NA50 by
one base; and achieved 99.69% precision among truth-resolved edits (1,591 fixes,
five harmful edits). Aggregate mismatch density increased slightly from 493.45
to 494.22 per 100 kbp while additional sequence aligned, so edit-level causal
auditing is reported alongside assembly-level alignment metrics.

The layout checkpoint produced a safe improvement:

1. Candidate diagnostics showed that identity rejection was not the dominant
   bottleneck; conservative per-base extension support prevented useful layout.
2. Ranked extension substantially improved continuity but introduced one
   misassembly at 500k.
3. Boundary support reduced continuity without removing the error.
4. Reciprocal-best filtering retained N50/NA50 136/134, but the same 240 bp
   translocation survived.
5. Disabling contig merging removed the translocation while retaining N50/NA50
   133/132 and increasing aligned length to 4.90 Mb.

The accepted research configuration is therefore ranked read extension with
extension consensus, reciprocal-best filtering, and
`--max-contig-iterations 0`.

Damage-aware beta-posterior ranking passed as a modest accuracy and recovery
improvement when used with no additional confidence margin. On EMN001 500k it
changed 12,433 length-only winners while preserving N50 133 and zero
misassemblies. Aligned length increased from 4,896,280 to 4,913,114 bp,
mismatches fell from 766.79 to 753.64 per 100 kbp, and indels fell from 1.70 to
1.51; NA50 decreased from 132 to 131. A 0.01 margin was rejected at the 100k
gate because 2,094 abstentions reduced N50 from 128 to 123. The validated
default is therefore zero: reciprocal-best filtering remains the structural
abstention gate.

This result does not materially improve continuity. A report-only recruitment
audit then tested both unassigned reads and reads borrowed from other clusters
on EMN001 100k. Reciprocal validation retained 61 deferred and 73 cross-cluster
assignments, but no frozen contig reached five-read consensus support; projected
extension was zero bases. The output checksum was unchanged, runtime was 12.11
seconds, and peak RSS was 234,532 kB. This closes read recruitment as the
current bottleneck.

The next go/no-go experiment is a read-supported contig-link audit. It must
preserve the assembly checksum while reporting candidate overlaps, reciprocal
links, independently supported linear chains, ambiguous branches, cycles, and
projected merged bases. Only a meaningful population of unambiguous supported
chains justifies an active join stage and MetaQUAST validation.

Further polishing, looser identity thresholds, and indiscriminate read reuse
are not current priorities: the experiments show that they cannot create long
paths and can trade strain fidelity for small local gains.
