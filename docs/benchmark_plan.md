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
