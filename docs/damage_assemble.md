# Damage-aware string graph command

`anvaya damage-assemble` is an experimental Phase 3 interface for merged ancient
DNA FASTQ fragments. It uses one damage-aware string-graph layout and one final
raw-evidence consensus pass. Coordinates in the placement and decision reports
are 0-based; raw intervals are half-open.

Required inputs are a FASTQ file and a CarpeDeam-format profile prefix resolving
to `<prefix>5p.prof` and `<prefix>3p.prof`. Outputs must be distinct paths:

```bash
anvaya damage-assemble \
  --input small-test.fq.gz \
  --profile-prefix profile \
  --output-before before.fasta \
  --output-consensus consensus.fasta \
  --output-assembled assembled.fasta \
  --output-unresolved unresolved.fasta \
  --diagnostics diagnostics.json \
  --placements-report placements.tsv \
  --consensus-report decisions.tsv
```

The assembled and unresolved outputs are optional but must be requested together.
They partition the consensus FASTA without changing sequences: assembled records
have at least two distinct, uniquely placed raw molecules; all others are written
to the unresolved sidecar. Stable `unitig_N` identifiers link both subsets to the
all-record outputs and placement report. Diagnostics record both partitions.

The command accepts variable 31–200 bp fragments and preserves unique exact
containments. Its default maximum is 1,000 reads. This cap is intentional: the
localized containment candidate failed independent strain-safety validation and
was reverted. The reader stops at the first record beyond the cap instead of
retaining the entire input. Bounded synthetic profiling passed at 10,000 reads after
exact-containment discovery was indexed and verified against the exhaustive
implementation. An exactly 10,000-record EMN001 smoke then passed at 5.57 seconds
and 83,924 KiB, but remained mostly read-through: 9,343 contigs, N50 61 bp versus
60 bp for the input reads, and zero consensus changes. The next allowed operational
gate was an exactly 25,000-record subset with `--max-reads 25000`; it passed at
19.26 seconds and 195,328 KiB but remained 90.81% read-through with only a 3.3%
N50 gain. A current-source reference-proxy audit selected 565 of 610 expected
directed dovetails (92.62%) on its small high-confidence subset. An unchanged,
exactly 50,000-record run is now allowed with `--max-reads 50000`, stopping above
120 seconds or 1,048,576 KiB RSS. It must not yet be used for the full 100k EMN001
comparison, and a smoke result is not evidence of improved assembly accuracy or
contiguity.

The 50k run passed resources and output/read count but failed the N50-response
gate (62 to 65 bp). Its exact-cohort mapping proxy found 535 overlap components
among 578 eligible reads, with a longest optimistic component of 178 bp and
65/70 candidate recall. Do not scale this ordered prefix further. Current
contiguity development uses the truth-known coverage ladder documented in the
implementation guide.

The Phase B unresolved-containment candidate passed development but failed its
preregistered independent synthetic validation. It is preserved as evidence under
`results/phase-c-unresolved-containment-validation-v1/`; production behavior was
restored to the Phase A baseline. No validated contiguity improvement is claimed.
