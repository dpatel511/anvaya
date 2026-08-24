# Essential literature notes

## CarpeDeam

Kraft et al. (2025), *CarpeDeam: a de novo metagenome assembler for heavily damaged ancient datasets*. The closest assembler comparator. It uses sample-specific damage information in a greedy, PenguiN-derived overlap workflow rather than a de Bruijn graph. It can improve recovery but has an accuracy-sensitivity and misassembly trade-off.

DOI: <https://doi.org/10.1186/s13059-025-03839-5>

## PyDamage

Borry et al. (2021), *PyDamage: automated ancient damage identification and estimation for contigs in ancient DNA de novo assembly*. PyDamage estimates damage after reads have been aligned to assembled contigs. It authenticates contigs but does not use damage evidence during graph construction or simplification.

DOI: <https://doi.org/10.7717/peerj.11845>

## Current aDNA assembly workflow

Fellows Yates et al. (2026), *De novo assembly and authentication of ancient DNA metagenomes with nf-core/mag*. Describes current assembly, authentication, and post-assembly correction. It does not introduce damage-aware DBG decisions.

DOI: <https://doi.org/10.1371/journal.pcbi.1014591>

## General DBG comparators

- MEGAHIT: succinct, iterative de Bruijn graph metagenomic assembler.
- metaSPAdes: metagenomic multi-k assembly and graph simplification.

## Working gap

Existing tools either assess damage after assembly or use it in overlap-based assembly. The targeted search has not identified an ancient-metagenomic DBG assembler that retains original-molecule positional evidence and uses it directly for graph construction, simplification, and path selection.

This is a working novelty assessment and must be rechecked before publication.
