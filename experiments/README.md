# Simulation experiments

Realistic benchmark data will be generated with established tools such as Gargammel or NGSNGS. Generated FASTQ files and result directories must not be committed.

## Development smoke test

`01_graph_disturbance.py` is a deterministic implementation check, not a scientific benchmark. It uses a random 500 bp reference, 60 bp reads placed every 10 bp, `k=21`, and terminal damage rates of 0.50, 0.30, 0.15, 0.08, and 0.04 from each end.

Run it with:

```bash
PYTHONPATH=src:tests python3 experiments/01_graph_disturbance.py
```

Results for seed 42:

| Condition | Nodes | Edges | Branches | Unitigs | Longest unitig |
|---|---:|---:|---:|---:|---:|
| Clean | 481 | 480 | 0 | 1 | 500 bp |
| Terminal damage | 531 | 530 | 24 | 49 | 85 bp |

The clean reads exactly reconstruct the reference. Twenty-seven terminal damage events fragment the graph and reduce the longest unitig from 500 bp to 85 bp.

This single-reference, single-seed result only validates the experimental plumbing. Publication-level experiments will use multiple replicates and established simulators with recorded versions, commands, configurations, seeds, and reference checksums.

## Graph traversal optimization smoke benchmark

The CLI was run with `k=21` on a local gzipped FASTQ containing 19,643 reads. This is a performance smoke test, not an assembly-quality benchmark.

| Measurement | Before | After |
|---|---:|---:|
| Total runtime | 1324.06 s | 8.74 s |
| Unitig extraction | 383.00 s | 0.24 s |
| Peak memory after optimization | — | 132,724 KiB |

The optimized run was approximately 151 times faster. Both runs produced 68,607 nodes, 72,918 edges, 8,082,223 k-mer observations, 5,183 branching nodes, and 9,704 unitigs. The resulting FASTA files were byte-identical.

The improvement came from calculating all node degrees in one graph pass, using constant-time degree lookups during traversal, and avoiding duplicate unitig extraction. The original graph-stage timing included expensive statistics and an implicit unitig extraction, so its stage-level time is not directly comparable with the new separated timings.
