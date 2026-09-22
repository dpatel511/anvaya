# Audited CarpeDeam comparator

`anvaya carpedeam-assemble` runs an external CarpeDeam executable as a pinned
benchmark comparator. Anvaya validates inputs and the damage profile,
forces CarpeDeam safe mode, records the exact binary and inputs, and summarizes
the emitted FASTA. It does not copy, modify or vendor CarpeDeam source.

This command measures a strong external assembler after internal containment and
post-layout linking mechanisms failed their frozen gates. It is not Anvaya's
assembly engine and must not be presented as an Anvaya assembly result. The
internal `damage-assemble` command remains the prototype under development.

## Requirements

- merged FASTA or FASTQ fragments, plain or gzip compressed;
- every fragment at least 20 bp, checked by a streaming scan;
- a CarpeDeam profile prefix resolving to `<prefix>5p.prof` and
  `<prefix>3p.prof`;
- an installed CarpeDeam executable exposing the documented
  `ancient_assemble` safe-mode options;
- new output, temporary-directory and diagnostics paths.

The executable is not a Python dependency. CarpeDeam is distributed separately
under GPL-3.0. Pin the external release and retain the binary SHA-256 written by
Anvaya. The repository's historical comparator used CarpeDeam v1.0.1 source at
commit `171e5c7a2968023361b75af09353ec2487c8772c`; a new binary must be identified
by its own recorded hash rather than assumed to match that checkout.

## Command

First verify the installed executable yourself:

```bash
carpedeam ancient_assemble -h
```

Then run through Anvaya:

```bash
anvaya carpedeam-assemble \
  --input merged.fastq.gz \
  --profile-prefix /path/to/sample_damage_ \
  --output results/carpedeam-safe/contigs.fasta \
  --temporary-directory results/carpedeam-safe/tmp \
  --diagnostics results/carpedeam-safe/diagnostics.json \
  --executable /absolute/path/to/carpedeam \
  --threads 2 \
  --min-contig-length 31
```

The command always supplies `--unsafe 0`. It does not expose an unsafe-mode
switch. It explicitly pins five raw-read iterations, ten total iterations,
merge identity 0.99 and safe coverage five instead of inheriting executable
defaults. Output and temporary paths must not exist, preventing accidental resume
or overwrite of an earlier run.

The diagnostics JSON records:

- the exact argument vector and elapsed time;
- executable path, executable SHA-256 and help-text SHA-256;
- compressed input-file SHA-256 plus streamed record/length statistics;
- both profile-file SHA-256 values;
- safe-mode settings, stdout/stderr log paths and return status;
- output FASTA SHA-256, contigs, bases, N50 and longest contig after success.

Coordinates do not appear in this interface. CarpeDeam's FASTA output does not
expose raw-read placements, so Anvaya must not claim molecule provenance,
damage-aware repolishing, strain retention or coordinate safety from this run.
Those require a separately validated mapping/evaluation stage. Reference-backed
evaluation remains external to assembly.

## Acceptance checkpoint

On the same frozen input cohort, compare this safe-backend output with the
unchanged internal baseline using identical minimum-contig and evaluator settings.
Report N50, NA50, auN/auNA, longest alignment, aligned length, genome fraction,
mismatches, indels, unaligned length and misassemblies. Hash equality of inputs,
profiles and references is mandatory. A longer N50 alone is not acceptance.

Actual FASTQ execution is user-run under the repository workflow. Unit tests use
a mocked executable and small synthetic FASTQ/profile fixtures.

For the frozen 50k checkpoint, run from WSL Bash exactly once:

```bash
set -euo pipefail
cd /mnt/e/Projects/anvaya
RUN=results/carpedeam-safe-emn001-50000-v1
mkdir "$RUN"
/usr/bin/time -v -o "$RUN/time.txt" \
  env PYTHONPATH=src /root/miniforge3/envs/anvaya/bin/python -m anvaya \
  carpedeam-assemble \
  --input results/damage-string-graph-emn001-smoke-50000-v1/reads-50000.fq.gz \
  --profile-prefix results/raw-consensus-100k/profile \
  --output "$RUN/contigs.fasta" \
  --temporary-directory "$RUN/tmp" \
  --diagnostics "$RUN/diagnostics.json" \
  --executable /root/miniforge3/envs/carpedeam/bin/carpedeam \
  --threads 2 \
  --min-contig-length 31 \
  > "$RUN/driver-stdout.txt" \
  2> "$RUN/driver-stderr.txt"
```

`mkdir` intentionally fails if this run directory already exists. Do not delete
or reuse an earlier run to make the command pass. After this completes, inspect
the saved diagnostics and timing before running the matched reference evaluation.
