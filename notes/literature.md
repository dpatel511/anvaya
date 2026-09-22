# Literature and method boundaries

These references motivate the overlap research direction. They do not establish
equivalence between Anvaya's implementation and any published method.

See [the 2026-09-06 survey and prototype roadmap](../docs/prototype_roadmap.md)
for primary-source checks, pinned CarpeDeam/PenguiN source inspection, reuse
decisions, benchmark limitations and the next implementation milestones.

## Assembly comparators

- Kraft et al. (2025), *CarpeDeam: a de novo metagenome assembler for heavily
  damaged ancient datasets*. The principal damage-aware comparator uses
  iterative fragment clustering and sample-specific damage likelihoods.
  Anvaya currently uses fixed-penalty overlap ranking and count-based consensus.
  DOI: <https://doi.org/10.1186/s13059-025-03839-5>.
- PenguiN motivates joint consideration of overlap length and mismatch rate:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11443906/>.
- StrainXpress motivates conservative resolution of competing overlap paths:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC9508831/>.
- MEGAHIT and metaSPAdes remain possible external metagenomic comparators;
  their inclusion does not require retaining an internal DBG assembler.

## Damage inference

- Borry et al. (2021), *PyDamage*, separates damage testing from reliability
  calibration on mapped contigs. DOI: <https://doi.org/10.7717/peerj.11845>.
- Jonsson et al. (2013), *mapDamage2.0*, models position-dependent molecular
  damage separately from sequencing qualities.
  DOI: <https://doi.org/10.1093/bioinformatics/btt193>.
- Skoglund et al. (2014), *PMDtools*, motivates comparison of damage and
  alternative explanations while accounting for quality and polymorphism.
  DOI: <https://doi.org/10.1073/pnas.1318934111>.
- Zhao et al. (2025), *ngsBriggs*, motivates uncertainty-aware damage estimation
  and separation of fitting from downstream classification.
  DOI: <https://doi.org/10.1111/1755-0998.14029>.

The retained `damage_likelihood.py` fits candidate-conditioned count profiles.
It is experimental standalone infrastructure, not a fitted profile for overlap
assembly. Importing reference-alignment or graph-selected models into a de novo
overlap workflow requires independent ascertainment and calibration validation.

## Working research gap

The active question is whether auditable raw-fragment overlap assembly can
improve recovery while protecting low-abundance strain alleles. Candidate
recruitment, per-base evidence and damage-aware consensus must be evaluated
separately. The targeted survey does not establish novelty; comparative performance
still requires a held-out benchmark assessment before publication.

## String-graph ambiguity boundary

- Myers (2005) defines transitive removal through equivalent spelled strings;
  Anvaya therefore compares complete direct and alternate path spellings before
  removing a damage-aware edge. DOI:
  <https://doi.org/10.1093/bioinformatics/bti1114>.
- Simpson and Durbin (2010) treat contained reads as redundant in an exact string
  graph. For damaged metagenomic reads, Anvaya applies the narrower rule that a
  contained sequence must have one exact parent/orientation/offset placement.
  DOI: <https://doi.org/10.1093/bioinformatics/btq217>.
- Baaijens et al. (2017), SAVAGE, resolves strain paths using strongly supported
  overlap groups and co-occurring mutations, with paired reads as additional
  linkage. This supports abstention when short merged fragments do not span a
  repeat or link diagnostic alleles. DOI: <https://doi.org/10.1101/gr.215038.116>.

An identical unspanned repeat copy cannot be assigned from sequence alone. The
controlled benchmark therefore reports whole-reference repeat coverage, but its
recovery gate uses sequence-resolved coverage outside explicitly recorded repeat
intervals. This prevents fragmented singleton reads from being scored as repeat
resolution while still detecting losses in inferable sequence.

## Phase before strain-aware path spelling

SAVAGE constructs overlap edges from statistically compatible haplotypic
sequence and uses co-occurring mutations to phase strains. StrainXpress separates
read clustering and local strain-aware assembly from its later master overlap
graph. Strainy calls informative variants, builds a read connection graph,
locally reassembles strain groups, and represents low-heterozygosity sequence as
unphased. These methods support a phase-first boundary for Anvaya: a conserved
overlap without a molecule spanning its flanking variants cannot select a strain
path.

Damage likelihood must precede variant selection. Terminal C→T/G→A observations
expected under post-mortem deamination must not become phasing markers merely
because both observed bases have support. PyDamage authenticates assembled
contigs; it does not provide missing molecule linkage for assembly.

- Baaijens et al. (2017), SAVAGE: <https://doi.org/10.1101/gr.215038.116>.
- Vicedomini et al. (2022), StrainXpress: <https://doi.org/10.1093/nar/gkac631>.
- Kazantseva et al. (2024), Strainy:
  <https://doi.org/10.1038/s41592-024-02424-1>.
- Borry et al. (2021), PyDamage: <https://doi.org/10.7717/peerj.11845>.
