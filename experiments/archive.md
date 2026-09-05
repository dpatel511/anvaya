# Research history

The repository originally explored directed and bidirected de Bruijn graphs,
terminal-evidence diagnostics, tip/bubble classification, graph-event likelihoods,
paired unitig extension and read threading. These implementations and experiments
01 through 18 were retired when the project became overlap-only.

The complete implementation, commands and chronological notebook are preserved
in Git at `d411b31` and earlier. For example:

```bash
git show d411b31:experiments/archive.md
git show d411b31:docs/benchmark_plan.md
```

## Conclusions retained from the retired work

- Coverage and sequence similarity alone could not reliably separate erroneous
  paths from low-abundance biological alternatives.
- Conservative damage-aware cleaning protected ambiguity but produced little
  public-data continuity improvement.
- Paired extension and read threading produced limited continuity gains and
  required explicit repeat protection.
- Singleton filtered-edge rescue did not establish useful safe recovery.
- Unrestricted read correction threatened rare-strain retention. Historical
  exact-recovery numbers should not be reused as overlap acceptance evidence;
  the old evaluator also mixed strand coordinate systems.
- Selected graph-locus damage profiles were not whole-library damage estimates.
  Cross-fitting and quality-aware evidence remain useful design principles,
  not a validated model for the current overlap assembler.

The retirement is a project-scope decision. It does not establish that every
DBG algorithm is unsuitable for ancient DNA.

## Overlap decisions retained

- Reciprocal ranked extension with contig merging disabled became the
  conservative comparison path.
- Contig merging introduced a reported 240 bp translocation for a small
  continuity gain; it remains experimental.
- An added confidence margin of 0.01 caused excessive abstention and was rejected.
- Frozen-layout polishing permits sequence-accuracy evaluation independently
  of layout changes; its audit and restoration scripts remain executable.
- Support-two increased breadth but requires further quality and accuracy
  validation. Exact graph-link tuning made only small continuity gains.
- Anchor length 11 was rejected on the TAF016 R1 diagnostic.

Current quantitative checkpoints and next decisions are maintained in the
[benchmark plan](../docs/benchmark_plan.md), not inferred from superseded runs.
