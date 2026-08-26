# Benchmark plan

## Development sequence

1. Toy sequence with an exactly known graph.
2. One simulated bacterial genome.
3. Two related strains at controlled abundances.
4. A small simulated bacterial and archaeal community.
5. Selected empirical ancient metagenomes.

Each simulation will separate sequencing error, postmortem damage, and genuine variation before combining them.

Controlled damage tests report both introduced events and damage-derived k-mer edges that survive graph filtering. Sensitivity is calculated against retained ground-truth edges so events removed by `min_count` are not incorrectly counted as detector failures. Empirical non-UDG/UDG comparisons use equal read counts, identical graph parameters, normalized topology rates, and replicate-level statistics.

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

Optimization begins only if damage-aware decisions reproducibly improve the accuracy-recovery trade-off without unacceptable strain collapse.
