# Literature and method boundaries

These references motivate the overlap research direction. They do not establish
equivalence between Anvaya's implementation and any published method.

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
separately. Novelty and comparative performance require a fresh literature and
held-out benchmark assessment before publication.
