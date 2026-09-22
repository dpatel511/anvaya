# Raw-fragment consensus experiment

This opt-in projection recalls bases on the **fixed progressive assembly layout**.
It does not change candidate discovery, overlap ranking, read ownership, extension,
contig membership or length. The primary and rescue projections continue to use
their existing consensus. This separation lets us test base calling before
changing assembly decisions. Trimming is not part of this experiment.

## Literature basis and adaptation

CarpeDeam uses a supplied damage profile for raw-fragment correction and extension,
and separates fragment-end information from contig merging. Its model and the
distinction between damage and sequencing error motivate this implementation;
Anvaya is not a port of CarpeDeam. [CarpeDeam Methods](https://doi.org/10.1186/s13059-025-03839-5).

mapDamage models positional damage and supports quality rescaling. Here we compose
a supplied damage transition with the **original** FASTQ error probability; we do
not first rescale quality and then apply damage again.
[mapDamage2.0](https://doi.org/10.1093/bioinformatics/btt193).

For each proposed true base t, observed raw base o, original position p and Phred Q:

`P(o | t,p,Q) = sum_b D(t,b,p) E(b,o,Q)`

`E` has diagonal `1 - 10^(-Q/10)` and off-diagonal `10^(-Q/10)/3`.
For the current untreated double-stranded model, `D` permits C→T at the original
5' end and G→A at the original 3' end. Other bases remain unchanged by damage.
Reverse placements complement the proposed contig base into the original read
orientation before evaluating the observation. Damage precedes sequencing error.
This is a conditional observation model, not a complete strain mixture model.

The model has a uniform prior over A/C/G/T and multiplies evidence across molecule
IDs, using log likelihoods for stability. It excludes missing qualities, Q<3 and
ambiguous bases. A bounded cache shares likelihood vectors across observations;
columns are accumulated and discarded one contig at a time.

## Conservative calling rules

These are provisional Anvaya rules, not thresholds claimed from the papers:

- Require at least three usable molecule IDs at a site.
- A change requires two Q20 observations of the proposed base and model posterior
  at least 0.99. Do not infer a base that lacks this observed support.
- If two alleles each have two Q20 observations for which every alternative true
  base gives likelihood below 0.01, flag the site and retain the existing base.
  This preserves ambiguity in the report; it does not reconstruct haplotypes or
  undo an earlier collapsed allele.
- Count a molecule once per column. Agreeing observations use the highest quality,
  then least damage-ambiguous observation. Discordant observations exclude that
  molecule at the site. A read placed at multiple offsets/orientations excludes
  its molecule from this contig's inference.

The nominal posterior is **not calibrated correctness**: the profile is treated
as known, molecules are treated as independent and the overlap selection process
is conditioned on. These assumptions require held-out evaluation. The single-input
loader currently assigns a distinct molecule ID per record; it does not detect
PCR duplicates or infer UMI families. Identical biological duplicates can therefore
still inflate support. Reference genomes are not used in consensus.

## Input profile and physical-end scope

`--raw-consensus-profile-prefix PREFIX` reads `PREFIX5p.prof` and `PREFIX3p.prof`
using the existing 12-column CarpeDeam profile format. The first row is distance
zero from the respective original read end. Values must be finite in [0,1].
This first model accepts only C>T in the 5p file and G>A in the 3p file; other
nonzero transitions are rejected rather than silently discarded. Outside supplied
rows, damage probability is zero. No damage parameters are fitted automatically.
Raw mismatch-frequency tables are not automatically equivalent to damage-only
transition probabilities; validate their provenance before conversion.

The CLI currently accepts this experiment only with one input file, intended to
contain merged untreated double-stranded fragments whose original molecule ends
have been retained. It cannot infer library preparation or whether upstream
quality trimming already removed those ends. Unmerged reads, single-stranded
libraries, UDG-treated libraries and already-rescaled qualities need explicit
model/preprocessing treatment; they are not validated by this experiment.

The existing EMN001 CarpeDeam benchmark log used its bundled `example/dhigh`
profile. It can be reused for a controlled comparison, but it is an **example
profile**, not an EMN001-specific estimate. The comparison script copies the
profile and records profile/source SHA256 hashes plus the Git revision.

## Placements and outputs

Tracking is enabled only when the consensus projection is requested. Each derived
record retains ungapped placements referring to immutable raw reads. Placements
include read index, orientation, original `[read_start, read_stop)` and the offset
of the full oriented read. All intervals are **0-based, half-open**:

- Forward: `contig_position = offset + original_position`.
- Reverse: `contig_position = offset + original_read_length - 1 - original_position`.

Left extension shifts offsets; trimming clips retained intervals without changing
the original read length used for damage distances. Later progressive extension
retains existing placements and adds the raw candidates used by its internal
consensus. Gapped placements are not supported yet. Untracked transformations clear
provenance, and this consensus refuses derived records without tracked placements.

`--raw-consensus-projection` writes the damage-model FASTA.
`--raw-consensus-quality-control` optionally writes the same caller with zero damage,
so we can separate quality weighting/guards from the damage profile's effect.
Both match the progressive projection's membership, ordering and lengths.

`--raw-consensus-report` writes changed, conflicting or rejected candidate decisions
with model posterior and 0-based position. `--raw-consensus-placements` writes the
source layout, with eligibility status for ambiguous mappings. `unitig_N` identifiers
match the length-filtered progressive FASTA. Placement eligibility does not imply
every base contributed: quality, ambiguity and duplicate rules still apply.

## Evaluation

Run `experiments/raw_consensus_comparison.sh` with merged FASTQ, profile prefix and
a new output directory. The user runs actual datasets; unit tests use small
synthetic data. Compare `before-contigs.fasta`, `quality-only.fasta` and
`damage-consensus.fasta` against the same reference/settings. Their N50 and lengths
must match by construction. Assess substitutions, unique coverage, per-contig
alignment changes and strain alleles; do not treat a changed alignment denominator
as proof of improved base accuracy. The scope is progressive contigs, not the full
support-two output evaluated in the trimming experiments.

### Fixed-alignment correction audit

`experiments/raw_consensus_audit.sh CONSENSUS_DIR REFERENCE_FASTA NEW_OUTPUT_DIR`
maps only the baseline contigs and runs `raw_consensus_audit.py` against those
fixed coordinates. It compares before→quality, before→damage and quality→damage
separately. No FASTQ is read and no assembly is rerun. Dependencies are minimap2,
samtools and optional evaluation-only `pysam==0.24.0`; the assembler does not
require pysam. Override binaries through `MINIMAP2`, `SAMTOOLS` and
`ANVAYA_AUDIT_PYTHON`. Run the wrapper from the repository root.

This follows QUAST's separation of reference-based accuracy and contiguity, with
an additional paired-site comparison specific to this experiment. Changes in
alignment can change the denominator of aggregate mismatch rates; holding the
baseline alignment fixed removes that particular confounder.
[QUAST](https://doi.org/10.1093/bioinformatics/btt086).

The wrapper uses minimap2's short-read preset with smaller seeds and lower chaining
and alignment score thresholds for these short contigs. These are exploratory
settings, not validated exhaustive mapping or a reproduction of MetaQUAST. It
retains secondary alignments and rejects any contig with multiple reported mapped
records, supplementary/secondary-only records, SA/XA alternatives, MAPQ below 30
or unavailable MAPQ (255). The MAPQ threshold is provisional; absence of a reported
alternative does not prove uniqueness. The 1 Gb index batch accommodates the
current 475 Mb reference; larger references need a single-index configuration
before interpreting MAPQ. [minimap2 options](https://lh3.github.io/minimap2/minimap2.html),
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

`audit/sites.tsv` contains one row per changed site per comparison, including
reference coordinates and reference base oriented to the original contig.
Positions are **0-based**, intervals **half-open**. Reverse alignments complement
the reference base and reverse the query position. Insertions and soft-clipped
positions have no reference base and remain unresolved. Hard-clipped alignments
are excluded. Ambiguous bases, missing alignments and ambiguous mappings remain
explicitly unresolved. `audit/summary.json` reports gains/losses of reference
agreement, substitutions where neither base matches, and unresolved counts.
The resolved-change fraction includes the neither-matches category; it is not
calibrated correction precision or a statistical significance test.

The audit checks FASTA ID/length joins, duplicate IDs, reference header lengths and
baseline alignment sequence identity. It streams alignment groups and fetches
FASTA sequences one contig at a time through pysam; reference fetches cover aligned
blocks, not potentially large skipped spans. Memory scales with contig/record IDs
plus one contig and its reported alignments. A new output directory is required.
Input/source hashes, dependency versions and mapping logs accompany the results.
[pysam API](https://pysam.readthedocs.io/en/latest/api.html).

The first EMN001 progressive comparison reported 2,290 damage-model substitutions
versus 2,263 quality-only substitutions, with 576 versus 853 conflict flags.
MetaQUAST mismatch rates were 1,272.76 before, 952.29 quality-only and 937.39 with
damage per 100 kb. These aggregate improvements motivate the audit; they do not
establish the correctness of individual substitutions.

Synthetic tests exercise mixed internal alleles across support ratios and profile
strengths, and explicitly document a counterexample: three terminal T observations
plus two internal C observations can represent either damaged C molecules or a
true strain mixture. A supplied 0.4 C→T profile changes the call to C and removes
the conflict flag even when the mixture interpretation is true. No single-column
caller can distinguish those identical observations without additional evidence.
These are deterministic mechanism tests, not held-out biological validation.
Reference agreement also cannot establish minor-strain preservation. Evaluate
held-out simulated haplotypes and profile misspecification before promoting this
caller to rescue or extension decisions.

### Alternative-mapping site audit

The strict `comparisons`, `baseline_alignment_status` and `sites.tsv` remain
available. `alternative_comparisons` and `alternative-sites.tsv` add two separate
conditional reference-agreement views:

- `all_reported`: every reported primary/secondary alignment must cover the changed
  position with the same unambiguous, contig-oriented reference base.
- `near_best_95pct`: apply the same rule to alignments with positive-best-score
  `AS >= 0.95 * best_AS`. This is an explicitly provisional sensitivity analysis,
  not a probability threshold or a validated definition of a credible mapping.
  Comparing it with `all_reported` exposes results dependent on excluding weaker
  mappings. Scores are from baseline alignments and can favor the baseline allele.

No majority vote is taken across reference genomes. Different loci can agree on a
base without identifying a unique source strain. A selected alignment that does
not cover the position leaves that site unresolved. Reverse-strand reference bases
are complemented; positions remain **0-based**, intervals **half-open**. Secondary
records may omit SEQ; their full-query CIGAR length is checked against the baseline
FASTA and their coordinates use the recorded strand. Present sequences must match.

Missing scores, split/supplementary alignments, unexpanded SA/XA alternatives and
missing/multiple primaries remain unresolved. A lone low/unknown-MAPQ alignment is
not rescued. Multiple low-MAPQ mappings can support only the stated conditional
agreement across reported candidates. Reaching the mapper's secondary-output cap
also leaves these views unresolved; `--secondary-limit` must equal the mapping
run's `-N` (20 for the existing wrapper). Minimap2 also filters candidate chains
and secondary score ratios before reporting, so even below that cap these sets
are not exhaustive. [minimap2 mapping options](https://lh3.github.io/minimap2/minimap2.html).

`alignments.tsv` exposes every record's AS, NM, MAPQ, flags, reference locus, strand
and CIGAR for inspection. Reference projection is shared across the two alternative
policies and computed only for contigs with substitutions. Existing name-grouped
BAMs can be passed directly to the Python audit, without remapping or reassembly;
use a new output directory. The summary records the current audit source hash and
pysam version, including when the BAM comes from an older run.

The local EMN001 download contains renamed reads, reference sequences, comparator
assemblies and subsets. No separate origin-coordinate or haplotype metadata file
was found in that directory inventory. Renamed IDs have not been interpreted as
truth coordinates. These reference-agreement reports therefore remain distinct
from simulator-origin validation.

### Competing-allele and neighbouring-evidence diagnostics

`--raw-consensus-mixture-report mixture.tsv` and
`--raw-consensus-linkage-report linkage.tsv` inspect the same deduplicated molecule
columns used by consensus. They do not change bases, thresholds, placements,
membership or extension. The comparison shell script now requests both reports;
use a new run directory. Existing input/profile/source hashes still apply.

The motivation is to compare competing sequence-source explanations, as in
[schmutzi](https://doi.org/10.1186/s13059-015-0776-0), and investigate overlap-linked
variation, as in [StrainXpress](https://doi.org/10.1093/nar/gkac543). This diagnostic
is an Anvaya adaptation, not either method's algorithm: it neither estimates human
contamination nor reconstructs strain haplotypes.

For raw observation likelihoods L_i(a), the single-base model evidence is
`Z1 = mean_a product_i L_i(a)`. The two-base evidence is
`Z2 = mean_(a<b) mean_f product_i [f L_i(a) + (1-f) L_i(b)]`.
The within-model priors are uniform over four single bases, or six unordered pairs
and 19 discrete fractions f = .05, .10, ..., .95. This finite prior is explicit:
it is not a continuous integral or maximum-likelihood estimate, and does not model
fractions below .05 or above .95. Averaging accounts for the additional parameter
choices; using only the best mixture fit would unfairly favour the mixture model.

`log_bf_mixture_vs_single` is natural-log(Z2/Z1). Positive values favour the mixture
under these assumptions; negative values favour a single allele. No between-model
prior, posterior probability, significance threshold or multiple-testing procedure
is applied. The reports are selected for columns with at least two observed alleles
among usable molecules, so their score distribution is also affected by selection.
All molecules share the supplied profile and are conditionally independent, as in
the caller. Misrecruitment, systematic error, correlated molecules and profile
misspecification remain alternative explanations; a positive score does not prove
two strains. Missing/Q<3 bases and duplicate/conflicting molecules retain the
caller's exclusions. Sites with fewer than three usable molecules are marked
`insufficient_support`; larger counts alone do not establish reliability.

The report includes baseline and projected bases, molecule counts, observed and
robust counts in A,C,G,T order, and the highest-evidence single base and pair.
Those best labels do not assert uniqueness; ties use deterministic base order.
An uncertain retained base is therefore distinguishable from an introduced change.
Pure observed columns are omitted even if a model might assign nonzero mixture
evidence. Tests compare the evidence against independent direct enumeration and
exercise terminal damage, internal mixtures, low support and high depth. These
are mechanism checks, not biological calibration or a held-out assembly benchmark.

Neighbouring sites are selected independently of the focal base: at least two
robust Q20 molecules must support each of at least two alleles. Robust here retains
the caller's rule that every alternative base likelihood is below .01. For each
focal site, `linkage.tsv` reports counts of observed Q20 focal alleles paired with
robust neighbouring alleles on the same molecule. The focal position itself is
excluded. Both positions must vary among the shared molecules. Counts, including
discordant combinations, are exported without a phasing decision or statistical
association claim. A nonzero `linked_neighbour_sites` count means eligible
co-occurrence evidence exists, not that linkage is established. Zero may simply
reflect short reads or inadequate support. All positions are **0-based**.

Likelihoods are evaluated in log space with identical vectors grouped by count.
For U distinct vectors at a conflicting column, the finite model costs 114*U
mixture evaluations, plus four single-base evaluations per vector. Neighbour
checks operate within one contig and only on independently selected variable
neighbours; their cost grows with conflicting sites, eligible neighbours and
shared molecule depth. Neither diagnostic runs unless a report is requested.
No runtime improvement on actual data is claimed.

### Experimental linked-allele guard

`--raw-consensus-linked-allele-guard` enables a conservative guard on the damage
projection only. It defaults off. The comparison script accepts an optional fourth
argument, `guarded`, to enable it; `baseline` preserves the previous behaviour.
Quality-only output is deliberately unguarded so it remains a stable control.

Only a proposed substitution that passes all existing consensus checks is eligible
for blocking. A witness must be a different position with independently selected
robust alleles (at least two robust molecules for each of at least two alleles).
All Q20 observations at the focal site must support exactly the baseline and
proposed alleles, each with at least two molecules and at least one robust
observation. Every one of these Q20 focal molecules must have a robust observation
at the witness. Their joint allele table must contain exactly two combinations,
different at both positions, each supported by at least two molecules. Missing
neighbour coverage, additional combinations, singleton groups and exclusively
damage-ambiguous focal support cannot trigger the guard.

The first qualifying witness in contig order blocks the change. Decisions record
`linked_allele_conflict` and `guard_neighbour_0based`; the summary counter is
`raw_consensus_linked_allele_rejections`. The witness is **0-based**, including a
valid witness at zero. Neighbours come from immutable raw columns before any
consensus edits, so calling order cannot manufacture linkage. Molecule counting
and exclusions are identical to consensus; this still does not identify PCR
duplicates. Only guarded runs build the guard's neighbour index, once per contig;
searching witnesses is limited to substitutions that would otherwise be accepted.

This is a provisional deterministic guard, not a phasing test or an association
p-value. It can retain errors in the baseline or correlated sequencing errors.
It does not use the mixture Bayes factor, change the posterior threshold, split
contigs or rescue blocked beneficial corrections. Tests reproduce a synthetic
two-group pattern and controls with discordant observations, absent overlap,
singletons, excluded molecules and genuine terminal damage. Actual evaluation
must count both prevented losses and blocked gains across every guarded change,
as well as the previously identified 117 damage-specific sites.

### Frozen synthetic validation v1

`PYTHONPATH=src python experiments/linked_guard_validation.py --output NEW_DIR`
runs generated fixtures only. It records source hashes and refuses an existing
output directory. The guard was not changed or tuned during this experiment.
The grid uses eight new seeds (73001–73008), 6/12/24 fragments, competing-haplotype
probabilities 0/.2/.5/.8, terminal damage 0/.2/.4, matched or overestimated profiles,
variant separations 6/35 bp, and Q20/Q30: 2,304 cases. Reused seeds and related
parameter settings are correlated; these are not 2,304 independent replicates.

This is a small, purpose-built fixed-layout mechanism benchmark, not a realistic
replacement for an external simulator or end-to-end assembly evaluation. It uses
known placements, random 31–65 bp fragments in both orientations, independent
damage followed by sequencing errors, and a 120 bp designated source haplotype.
A competing haplotype differs at two sites. The baseline contains independent 2%
substitutions so indiscriminately preserving it is not optimal. The target source
is designated explicitly, not inferred from a seed read. No indels, repeats,
placement errors, correlated artefacts or PCR families are simulated. Source-allele
retention does not measure recovery of both strains. Counts include only positions
covered by at least one generated fragment. All coordinates are 0-based/half-open.

Results in `results/linked-guard-synthetic-v1/summary.json` and `cases.tsv`:

- 178,428 covered positions across cases; baseline errors 3,456.
- Unguarded errors 1,644; guarded errors 1,634.
- Ten blocked changes prevented source-relative errors; zero blocked beneficial
  corrections and zero changes where both alternatives were wrong.
- Source alleles retained at covered variant sites: 2,856 → 2,866 of 3,456.

All ten blocks occur for seed 73002, Q20, competitor probability .8, nearby variants
and depth 6/12, across related damage/profile settings. None occur for the other
seven seeds. The apparently clean zero-loss result therefore provides only narrow
evidence, not a calibrated error rate or a reason to enable the guard by default.
Independent tests check forward/reverse truth coordinates, profile-only changes,
reproducibility and the equality between error reduction and prevented errors minus
blocked beneficial corrections. The full suite passed 212 tests. Further validation
needs independent templates and realistic recruitment/artefact errors; the rule
remains frozen and opt-in.
