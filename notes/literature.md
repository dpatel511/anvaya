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

## Working gap

Existing tools either assess damage after assembly or use it in overlap-based assembly. The targeted search has not identified an ancient-metagenomic DBG assembler that retains original-molecule positional evidence and uses it directly for graph construction, simplification, and path selection.

This is a working novelty assessment and must be rechecked before publication.
