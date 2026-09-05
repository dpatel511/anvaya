# Research question

Can a reference-free overlap assembler distinguish postmortem damage and
repeat-induced overlaps from genuine biological variation, improving ancient
bacterial and archaeal metagenome recovery without false joins or loss of
low-abundance strains?

Damage-compatible substitutions near physical molecule ends should receive
different support from ordinary errors and alleles observed across independent
fragments. Original sequencing qualities and molecular positions should remain
available after layout and correction.

The current implementation explores this with bounded overlap discovery,
fixed-penalty damage-aware ranking, count-based consensus and experimental
recovery projections. Calibrated damage inference integrated into overlap
decisions and general biological validation remain research goals.
