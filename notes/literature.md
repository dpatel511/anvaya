# Essential literature notes

## CarpeDeam

Kraft et al. (2025), *CarpeDeam: a de novo metagenome assembler for heavily damaged ancient datasets*. The closest assembler comparator. It uses sample-specific damage information in a greedy, PenguiN-derived overlap workflow rather than a de Bruijn graph. Its nucleotide likelihood marginalizes over the latent post-damage base and then applies a symmetric sequencing-error emission with fixed error rate 0.001. Anvaya adopts that latent-base structure but substitutes each observation's Phred probability and scores held-out graph loci. CarpeDeam can improve recovery but has an accuracy-sensitivity and misassembly trade-off.

DOI: <https://doi.org/10.1186/s13059-025-03839-5>

## PyDamage

Borry et al. (2021), *PyDamage: automated ancient damage identification and estimation for contigs in ancient DNA de novo assembly*. PyDamage estimates damage after reads have been aligned to assembled contigs. It combines a damage-versus-null likelihood-ratio test, multiple-testing correction, and a simulated-data accuracy model driven by damage, coverage, and contig length. Anvaya borrows the separation between evidence testing and reliability calibration, but calibrates local graph events and remains non-destructive.

DOI: <https://doi.org/10.7717/peerj.11845>

## mapDamage 2.0

Jónsson et al. (2013), *mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters*. The model estimates sample-level, position-dependent damage parameters and can downscale base qualities at likely damaged sites. Anvaya borrows the separation between the original sequencing quality and the damage process, but cannot directly reuse its reference-alignment likelihood in a de novo graph.

DOI: <https://doi.org/10.1093/bioinformatics/btt193>

## PMDtools

Skoglund et al. (2014), *Separating endogenous ancient DNA from modern day contamination in a Siberian Neandertal*. PMDtools compares postmortem-damage and null explanations while accounting for base quality and biological polymorphism. This motivates Anvaya's three-way damage/error/variation likelihood. Unlike PMDtools' reference-aligned per-read score, Anvaya uses graph-linked molecule observations and class-conditional conformal calibration.

DOI: <https://doi.org/10.1073/pnas.1318934111>

## ngsBriggs

Zhao et al. (2025), *Revisiting the Briggs Ancient DNA Damage Model: A Fast Maximum Likelihood Method to Estimate Post-Mortem Damage*. ngsBriggs combines sample-level damage estimation with posterior ancient-versus-modern read classification and fragment-length evidence. It reinforces the need to separate profile fitting from downstream scoring and to expose uncertainty. Anvaya now propagates cross-fit curve instability and calibrates decisions empirically, but still does not claim whole-library parameters or posterior probabilities.

DOI: <https://doi.org/10.1111/1755-0998.14029>

## Current aDNA assembly workflow

Fellows Yates et al. (2026), *De novo assembly and authentication of ancient DNA metagenomes with nf-core/mag*. Describes current assembly, authentication, and post-assembly correction. It does not introduce damage-aware DBG decisions.

DOI: <https://doi.org/10.1371/journal.pcbi.1014591>

## General DBG comparators

- MEGAHIT: succinct, iterative de Bruijn graph metagenomic assembler.
- metaSPAdes: metagenomic multi-k assembly and graph simplification.

## Overlap-ranking comparators

- PenguiN compares overlap error-rate distributions so overlap length
  contributes statistical evidence rather than winning independently of the
  mismatch rate. This motivates Anvaya's beta-posterior lower-confidence score:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11443906/>.
- StrainXpress balances overlap identity and length and joins only unambiguous
  paths after branch removal. This supports retaining Anvaya's reciprocal
  structural gate separately from its candidate-ranking score:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC9508831/>.
- Anvaya's implementation is a clean-room approximation rather than a copy of
  CarpeDeam or PenguiN. It discounts terminal C/T and G/A differences, applies
  full penalties elsewhere, and ranks the resulting beta-posterior lower bound.
  EMN001 testing showed that an additional 0.01 confidence margin was overly
  conservative; zero margin improved alignment accuracy while reciprocal-best
  filtering continued to enforce structural abstention.

## Damage-aware consensus implementation choices

- CarpeDeam clusters overlapping fragments using shared k-mers, filters them at
  90% sequence identity and 99.9% R/Y identity, and uses a sample-specific
  damage likelihood for correction and extension. Its safe consensus mode
  reduces but does not eliminate misassemblies. Anvaya's next backend will test
  the published iterative cluster, correction, and safe-consensus extension
  architecture through a clean-room implementation rather than copying the
  GPL-3 source. Unsafe extension is not an acceptance target:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12557918/
- BayesHammer shows why global abundance thresholds fail under nonuniform
  coverage and instead separates local Hamming-neighborhood clusters using base
  qualities and Bayesian penalties. Anvaya consequently requires local overlap
  consensus and quality-complete supporting observations rather than treating
  low-count k-mers globally:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3549815/
- Ancient-DNA work predating specialized assemblers notes that randomly
  distributed deamination can usually be excluded when at least three
  overlapping molecules support a consensus. Anvaya uses three independent
  internal molecules as its initial minimum and adds a 4:1 dominance gate:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2847228/
- Multi-assembly ancient-DNA workflows map reads back to candidate contigs and
  remove unsupported contigs before overlap-based merging. This supports the
  decision not to restore the stashed zero-read-support trusted high-k edges:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5384568/

## Working gap

The DBG evidence program remains a useful diagnostic contribution, but its
cleaning and correction paths have not produced a major public-data continuity
gain. The active gap is a scalable, auditable whole-fragment overlap backend
that preserves low-abundance strain paths while applying sample-specific damage
evidence during correction and safe extension.

This is a working novelty assessment and must be rechecked before publication.
