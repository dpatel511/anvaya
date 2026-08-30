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
must not become a production assembly policy. The result also establishes that
edge truth is insufficient as a contig-safety criterion. Further continuity
work should first measure perfect damage/error normalization and then evaluate
damage-aware consensus projection or multi-k construction.
