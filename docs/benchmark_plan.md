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

The audit completed without crossing that gate. At cluster support five, only
two of 152 reciprocal overlaps had two-read support and none had support three;
at cluster support three, only six of 741 reciprocal overlaps had two-read
support, two had support three, and none had support five. The support-three
projection contained six joins and 1,063 merged bases but did not exceed a
270 bp longest contig. Output checksums matched their respective non-audit
assemblies. Contig linking remains report-only and further threshold tuning is
closed.

The cluster-support sweep exposed a more useful Pareto boundary. Support three
recovered 1,312,891 aligned bases with 0.463% genome fraction and zero detected
misassemblies, approaching CarpeDeam safe at 1,431,468 aligned bases and 0.512%,
but reached only N50/NA50 113/112 versus 134/133. Support five retained the best
Anvaya continuity at N50/NA50 128/127 but only 591,906 aligned bases. The next
experiment is therefore a two-tier redundancy audit: freeze support-five
primary contigs, generate support-three rescue candidates, and project only
rescue representatives that are not contained by or redundant with a longer
primary or rescue contig. It must report primary checksum, candidate and
retained bases, containment and redundancy rejection counts, projected N50,
and projected longest contig before any output-changing mode is considered.

The two-tier audit completed and preserved the primary checksum. From 14,024
support-three candidates totaling 1,461,605 bp, it classified 4,766 as contained
by a primary, 71 as redundant with a longer rescue representative, 259 as
extending a primary, and 8,928 as novel. Adding only novel candidates projected
1,446,636 bp across 13,891 contigs, but N50 fell from 128 to 113 and the longest
contig stayed at 331 bp. The recovery gain therefore does not pass the
continuity gate.

Unique primary replacement also failed the material-improvement gate. Of 259
primary-extending candidates, three matched multiple primaries and abstained;
249 primaries received a longest unique replacement. The projection added
9,100 bp, increased N50 only from 128 to 129, and left the longest contig at
331 bp. Output remained byte-identical, runtime was 28.58 seconds, and peak RSS
was 276,772 kB. Replacement remains report-only.

The extension-round cap was tested independently. Six rounds increased total
length from 611,844 to 612,810 bp and the longest contig from 331 to 353 bp, but
left N50 at 128. Ten rounds produced the same contig count, total length, N50,
and longest-contig length as six rounds. Further round tuning is closed.

The next go/no-go experiment is global iterative reclustering with controlled
read reuse. Each iteration must report reused fragments, multi-contig conflicts,
unique confidence assignments, added bases, redundancy removals, projected N50,
projected longest contig, runtime, and memory. The 100k gate requires a material
continuity improvement over N50 129 and longest contig 353 before an active mode
or 500k/MetaQUAST evaluation is justified.

Further polishing, looser identity thresholds, and indiscriminate read reuse
are not current priorities: the experiments show that they cannot create long
paths and can trade strain fidelity for small local gains.

## Progressive overlap lifecycle checkpoint

The first progressive implementation corrected sequence ownership but extended
each admitted cluster only once, consumed all of its members, and deferred
further growth to a global unused-read pass. On EMN001 100k with one-anchor
candidate discovery, that version emitted 5,388 contigs totaling 673,486 bp,
with N50/NA50 130/129, a 266 bp longest contig, and 657,323 aligned bases. Its
global pass inspected 3,436 candidates but extended zero centers; 2,197 centers
failed the five-read consensus requirement.

The accepted multiround implementation keeps each center alive, reindexes it
after extension, recruits newly reachable unassigned raw reads, shifts retained
member alignments after left extension, and consumes members only after local
convergence. Cluster membership discovered for another center remains
unavailable, preserving deterministic longest-first ownership. On the same
100k input it recruited 3,012 additional reads across 8,123 local extension
rounds and retained the same 5,388 output representatives while increasing:

- total sequence from 673,486 to 717,381 bp;
- aligned sequence from 657,323 to 702,018 bp;
- N50 from 130 to 139 and NA50 from 129 to 137;
- auNA from 129.9 to 141.6;
- the longest contig and alignment from 266 to 325 bp; and
- genome fraction from 0.256% to 0.272%.

MetaQUAST reported zero misassemblies for both versions. Runtime increased from
93.64 to 100.76 seconds and peak RSS from 307,628 to 307,780 kB. The trade-off
was an increase in mismatches from 1,163.51 to 1,272.76 and indels from 2.74 to
3.99 per 100 kbp. The primary FASTA remained byte-identical, confirming that
the projection is isolated. The result passes the 100k continuity gate but must
remain experimental until it scales at 500k and its error trade-off is bounded.

The next go/no-go experiment starts from this multiround projection and audits
raw-supported continuation between exhausted centers. It must preserve exact
or reciprocal-best precedence, abstain at conflicting strain evidence, and
beat N50/NA50 139/137 and the 325 bp longest alignment without introducing a
misassembly or materially increasing the 1,272.76/3.99 mismatch/indel rates.

## Global layout and master-graph checkpoint

Iterative reclustering crossed its 100k gate, but immutable-evidence and
quality-aware mismatch restrictions alone did not produce the best continuity.
The accepted layout therefore preserves immutable raw reads as evidence and
separates sequence generation from graph projection. An exact master overlap
graph first removes transitive edges and spells reciprocal non-branching paths.
A second report-only graph then permits at most one near-exact mismatch when a
Q20 raw-molecule majority has support at least three, the competing allele has
support below two, and the support margin is at least two. Exact edges always
win; independently supported competing alleles abstain as strain conflicts;
terminal C/T and G/A observations inside five bases of a raw-read end do not
count as confirmation.

The 100k raw-confirmed projection improved the exact graph from N50/NA50
133/131 to 135/134, increased the largest alignment from 325 to 412 bp, reduced
duplication from 1.103 to 1.077, and retained zero misassemblies. Mismatch and
indel rates improved from 1,191.23/3.30 to 1,174.71/3.07 per 100 kbp. Genome
fraction decreased from 0.230% to 0.228%.

The result scaled at 500k. Of 6,042 near-exact candidates, 4,211 were accepted;
181 were rejected as strain conflicts and 276 lost to an exact edge. The graph
collapsed 8,239 inputs into 3,807 linear paths and resolved 1,760 overlap bases.
Compared with the exact graph, contigs fell from 49,099 to 44,667, N50/NA50
rose from 135/133 to 141/139, auNA rose from 134.8 to 142.9, and the largest
alignment rose from 474 to 603 bp. Misassemblies remained zero, mismatch and
indel rates fell from 824.01/1.50 to 814.60/1.38, and duplication fell from
1.307 to 1.244. The cost was a small reduction in genome fraction from 1.169%
to 1.157% and aligned sequence from 5,610,113 to 5,305,903 bp.

This passes the safe-continuity gate and is the current balanced projection.
The higher-continuity ranked-consensus result remains an upper bound at
N50/NA50 158/156, but it contains one misassembly and less aligned sequence.
The next go/no-go experiment must address the 767 ambiguous raw-confirmed graph
ends without relaxing mismatch evidence. It should rank competing end edges
with immutable raw molecule support, preserve exact-edge precedence and strain
conflicts, and require zero misassemblies plus an NA50 improvement at 500k.

## Evidence-gated recovery checkpoint

Progressive extension, raw-supported link audits, adaptive and selective rescue,
paired-read merging and scaffolding, and exact or raw-confirmed graph projections
are implemented as isolated projections: the primary overlap FASTA remains
unchanged unless a projected output is explicitly requested. Ambiguous ends,
competing molecule assignments, strain conflicts, insufficient base quality,
and non-reciprocal links abstain rather than forcing a join.

On the EMN001 100k input, the high-confidence support-two rescue is the current
recovery-oriented projection. It produced 16,137 contigs and 1,567,417 bases,
with N50/NA50 106/105, 1,320,332 aligned bases, 0.478% genome fraction, and zero
MetaQUAST misassemblies. The CarpeDeam safe comparator produced 12,491 contigs,
N50/NA50 134/133, 1,431,468 aligned bases, and 0.512% genome fraction. Anvaya
had fewer indels (2.80 versus 3.21 per 100 kbp), but more mismatches (988.31
versus 873.79), 2.57 times as much unaligned sequence, and substantially greater
fragmentation. Exact graph projection of the support-two output improved the
longest contig from 325 to 439 bp but only changed N50 from 106 to 107, so it is
not the selected recovery output.

The 100k support-two rerun reproduced the primary SHA-256 checksum and all
reported progressive and support-two assembly statistics. The full test suite
passed 353 tests after the implementation was assembled on this branch.

TAF016 R1 and R2 were also run independently through the single-input CLI as a
mate-ablation execution check. Both completed, but this is not a genuine
single-end benchmark: the files are mates from a paired-end library, and
discarding the other mate changes the available evidence. Progressive clustering
recruited only 2,916 R1 reads and 661 R2 reads from 100,000 inputs; support-two
N50 remained 76 bp. An anchor-length reduction from 15 to 11 on R1 increased
candidate offsets but reduced validated alignments from 331,088 to 89,786 and
clustered reads from 2,916 to 1,753. The shorter anchor is rejected, and these
mate-ablation measurements must not be presented as biological single-end
validation.

The next validation step is a genuine single-end ancient library or a directly
simulated single-end truth set, followed by the 500k support-two scaling run.
No further overlap-threshold or graph-link tuning is justified until those tests
separate insufficient sampling depth from candidate-discovery failure.
