# Overlap architecture

Anvaya assembles whole fragments using bounded overlap-layout-consensus stages.
The initial scope is ancient bacterial and archaeal metagenomes sequenced from
untreated double-stranded libraries, with merged fragments as the initial input
model. The assembler does not require a reference genome.

## Module responsibilities

| Module | Responsibility |
|---|---|
| `reads.py`, `sequences.py` | FASTA/FASTQ input, quality values, sequence normalization and orientation |
| `overlap_index.py` | Packed canonical anchor sketches and DNA/R/Y identity helpers |
| `overlap_assembly.py` | Initial clustering, overlap ranking, consensus, polishing and compatibility entry points |
| `overlap_progressive.py` | Immutable raw observations, representative lifecycle and raw-read recruitment |
| `overlap_adaptive_rescue.py` | Adaptive, selective and support-two recovery projections |
| `overlap_reclustering.py` | Fixed-pool iterative reclustering audits |
| `overlap_graph.py` | Exact, raw-confirmed and containment projections over assembled contigs |
| `overlap_progressive_links.py` | Raw-molecule-supported links between progressive contigs |
| `paired_reads.py`, `overlap_scaffolding.py` | Experimental paired merging and conservative scaffolding |
| `cli.py`, `output.py` | Workflow composition, diagnostics and deterministic FASTA output |
| `damage_likelihood.py` | Standalone experimental candidate-conditioned damage fitter; not used in assembly |

## Assembly and evidence

Candidate discovery sketches canonical anchors, proposes orientations and
ungapped offsets, and checks full overlaps in DNA and purine/pyrimidine space.
Repeated anchors and ambiguous molecule placements can be rejected. Candidate
discovery and overlap validation are separate from layout decisions.

Ranked extension can discount terminal C/T and G/A differences using a fixed
penalty. The resulting beta mean-minus-standard-deviation score is a ranking
heuristic, not a calibrated correctness probability. Consensus uses counts.
The standalone damage fitter does not provide a sample-specific profile to
these decisions.

Progressive records retain immutable raw reads separately from corrected and
extended representatives. Consumed records remain available as raw evidence.
The initial extension stage can require reciprocal support; later progressive
recruitment uses a different support gate. These policies must not be described
as a uniform reciprocal guarantee across every stage.

Raw reads and qualities are currently loaded in memory. Representative updates
copy the pool, and the index is rebuilt during some stages. Those are known
scaling targets, not changes made by the structural retirement.

## Projection and output contracts

Optional projections write separate FASTAs and must preserve the primary
output. Exact and raw-confirmed contig graphs are overlap graphs, not k-mer
assembly graphs, and remain part of the implementation.

Frozen-layout polishing changes bases without changing contig membership,
orientation, order or length. Its edit report supports one-edit-at-a-time
reference auditing and exact restoration of the pre-polish assembly.

The writer retains legacy `unitig_N` identifiers to preserve existing output
checksums and report associations. These names do not imply a retained DBG
implementation.

Read and contig offsets are 0-based. Internal sequence intervals and PAF
intervals are 0-based, half-open. External tools may report other conventions;
conversions must be explicit in any new evaluator.

## Current limitations

- Quality checks on support-two seed admission do not cover later recruited
  overhangs or every retained seed boundary.
- Flat contributing-molecule sets do not replace per-base raw-read placements.
- A read end or derived contig end is not necessarily a physical molecule end.
- Greedy read ownership and limited anchors may lose valid overlaps; their
  contribution to fragmentation has not been quantified independently.
- There is no complete truth-labelled overlap acceptance ladder yet.
- Paired inputs are associated by record order; callers must synchronize them.

## Structural retirement

The earlier directed/bidirected DBG assemblers, graph diagnostics, graph-event
calibration command and graph-only experiment runners are removed. Their code
and complete original notebook remain in Git history through `d411b31`.
Only the independently tested damage-likelihood fitter and active overlap
helpers were retained from that research infrastructure.

The cleanup changes imports and the available commands, not overlap thresholds,
consensus, selection policy or FASTA serialization. Future algorithm changes
must pass the [benchmark plan](benchmark_plan.md).
