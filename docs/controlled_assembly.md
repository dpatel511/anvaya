# Controlled progressive assembly benchmark

Run `PYTHONPATH=src python experiments/controlled_assembly.py --output NEW_DIR`.
This uses small synthetic fixtures only. Each case stores private origin truth
and injected damage events separately from the assembler. All intervals are
0-based, half-open. Molecules have unique IDs; no duplicate reads are simulated.

The fixed matrix has 36 cases: two seeds, unique/repeated/related-strain panels,
3x/10x/20x coverage per reference, and clean/damaged conditions. References are
1,200 bases, reads 80 bases, strands random. The repeat is 150 bases and strains
differ every 50 bases. Damage follows the supplied five-position C>T/G>A model.
Clean and damaged counterparts share origins. There are no sequencing errors,
indels, coverage bias or PCR effects; Q35 is a nominal consensus input. These
are mechanism controls, not realistic metagenome simulations.

The benchmark measures true candidate-offset recall, truth-pair ownership after
clustering, progressive extension diagnostics and origin-consistent contig
accuracy/recovery before and after one consensus pass. Consensus receives the
generating damage profile, an optimistic assumption. Optional rescues and links
are disabled. Origins never enter assembly/consensus; only evaluation uses them.
Conflicting origins remain unresolved for unique-reference recovery because repeat
placements and strain-shared reads do not establish a unique reference transform.
The scorer also reports per-contig candidate transforms using 0-based, half-open
coordinates; provenance and sequence-compatibility categories; diagnostic allele
positions; and `accuracy_evaluated_bases`, `accuracy_mismatches`, and
`accuracy_unresolved_bases`. Identical candidate sequences and exact matches to one
candidate can be accuracy-evaluable without resolving origin. A sequence that
matches diagnostic alleles from multiple references but no full candidate is
reported as `diagnostic_mosaic` and remains outside the accuracy denominator.
Every sequence category is reported in both contigs and bases. The base totals
matter because a stable count can conceal a longer unsafe contig.

The final repaired evaluator was run on development seeds 95101–95104 with matched
exact and damage-aware layouts in `results/controlled-string-evaluator-v5-dev-control`
and `results/controlled-damage-string-evaluator-v5-dev`. It projects recorded
damage events through raw placements, including reverse strands, and distinguishes
damage-compatible candidates from unexplained substitutions. Candidate mismatch
bounds do not select an origin. A separate 1% substitution-error control confirms
that ordinary errors are not relabeled as injected damage.

Across all 36 damaged paired cases, damage-aware layout matched or improved exact
layout for N50 and evaluated mismatch rate, added no diagnostic-mosaic or
unexplained-incompatible contigs or bases, and preserved all sequence-resolved coverage
outside explicitly identical repeat intervals. The earlier 185 bp seed-95101
repeat contig has one damage-compatible transform with one mismatch; its other
transform has 25 unexplained mismatches. Whole-reference copy coverage remains
reported, but it is not evidence that an unspanned identical repeat was resolved.

The one-time held-out evaluation on seeds 95201–95204 did not pass. Thirty-five
of 36 damaged cases preserved sequence-resolved nonambiguous coverage. In
`95202-strains-3-0.4`, exact layout produced one 135 bp diagnostic mosaic with
two conflicting diagnostic sites; damage-aware layout produced one 254 bp
diagnostic mosaic with five conflicting sites. The mosaic count remained one,
while nonambiguous resolved coverage fell from 2,189 to 2,062 bases. These seeds
are consumed validation evidence and must not be rerun as an independent gate.

Verified v2 results (`results/controlled-assembly-v2`, 36 cases): clean unique
20x controls have single truth-overlap components of 1,198 and 1,200 bases, but
produce 14/15 contigs with N50 92/88 and longest lengths 138/136. Candidate recall
is 7,556/7,684 and 7,985/8,080 directed overlaps. All resolved output bases match
truth in these controls. Clustering assigns 292/300 and 299/300 reads; 1,041 and
1,277 unordered true overlapping pairs cross cluster boundaries. Later raw
recruitment adds zero bases in both cases. v1 lacks ownership/connectivity counts;
use v2 for diagnosis.

This establishes a clean-input layout failure and makes cross-cluster evidence
the next intervention to test. It does not isolate ownership from every other
extension rule. Keep candidate and consensus settings fixed while testing a
layout change against these same controls, then validate on new seeds. Longer
outputs must remain origin-consistent and preserve base accuracy. Do not tune
on the related-strain cases and claim independent validation from them.

## Assigned-read reuse experiment

`--reuse-raw-evidence` enables an API-only experimental option in later
progressive extension. Candidate search sees all immutable raw observations;
support, identity, ranking and iteration limits stay fixed. Reads are not claimed
across centers and no other center is retired. Candidate search still selects
at most one placement per molecule for a center. Existing and new raw placements
are retained for consensus. There is no assembly CLI option or real-data wrapper
for this experiment.

Verified `results/controlled-reuse-v1` covers the same 36 cases. At clean unique
20x coverage, N50 increases from 92/88 to 273/273, longest contigs from 138/136 to
347/312, and resolved unique bases from 1155/1126 to 1181/1180. Base mismatches
remain zero. The clean related-strain controls, however, increase conflicting-
origin contigs from 11/12 to 28/26. Conflicting origins are not automatically
false sequence joins, but their increase prevents claiming improved accuracy.
Outputs also retain overlapping redundant contigs: summed added bases are not
novel recovery. Three extension rounds are a fixed experimental bound, not a
claim of convergence.

Conclusion: ownership limits growth in the unique controls, but unrestricted
reuse is not ready for promotion. The next layout intervention needs compatible
placement/allele constraints and a way to consolidate redundant paths. Do not
enable this mode on the real dataset based solely on the N50 increase.

## Exact reuse and compatible containment

`--constrained` is another synthetic-only comparison. It requires DNA identity
1.0 during later evidence reuse, then consolidates a contained contig only when
its exact sequence has one placement (considering both orientations), at least
two distinct shared molecules support the transform, and shared read placements
agree. Ambiguous placements or conflicting offsets block retirement. Original
raw intervals remain unchanged; reverse placement offsets are transformed in
0-based coordinates. The helper is quadratic and rejects pools over 1,000 active
contigs; it is not a production-scale algorithm. This is an exact-match control,
not a damage likelihood or full allele-phasing model.

`results/controlled-constrained-v1` ran all 36 cases. Clean unique 20x N50 remains
273/273, while clean strain conflicting-origin counts remain 25/24 (baseline
11/12). Damaged unique N50 is 173/202 versus unrestricted reuse's 233/214, and
post-consensus mismatches are 2/1 versus 0/0. Containment retires only 0/1 contigs
in the clean unique cases. This variant does not meet the accuracy/contiguity
acceptance criteria and remains experimental. The v1 hash manifest predates
adding the containment helper to the manifest; subsequent runs include it.

The result argues against continued threshold tuning or post-hoc containment
as the main solution. A next intervention must represent competing overlap
paths and preserve placement compatibility during growth; independent cluster
extension cannot supply that merely by demanding exact local matches.

## Damage-aware string-layout comparison

Run the fixed development matrix with:

```bash
PYTHONPATH=src python experiments/controlled_assembly.py \
  --damage-string-layout \
  --output results/controlled-damage-string-phase2-v1
```

Use `--seeds` for a separately declared seed set. The held-out Phase 2 run used
`--seeds 94011 94012`. `--damage-string-layout` cannot be combined with the
progressive reuse modes or `--string-layout`. It remains a bounded synthetic
experiment. Its source implementation now supports variable lengths, while this
fixed matrix deliberately retains 80-base reads as the Phase 1/2 oracle.
