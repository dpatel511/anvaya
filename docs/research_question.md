# Research question

Can a reference-free, evidence-annotated assembler distinguish
postmortem-damage-induced differences and repeat-induced overlaps from genuine
biological variation during de novo assembly of ancient bacterial and archaeal
metagenomes, improving continuity and genome recovery without introducing
misassemblies or collapsing low-abundance strains?

## Hypothesis

Damage-derived paths should be supported mainly by damage-compatible substitutions near molecule ends. Genuine variation should receive support across molecule positions, orientations, and independent fragments.

Anvaya retains this evidence in both its established de Bruijn graph diagnostics
and its active whole-fragment overlap backend. Competing paths are treated as
damage-derived, likely genuine, repeat-ambiguous, or unresolved; ambiguous paths
are preserved rather than joined or corrected aggressively.
