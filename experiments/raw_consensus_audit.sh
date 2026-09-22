#!/usr/bin/env bash
# User-run evaluation; never reads original FASTQ or reruns assembly.
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo "Usage: bash experiments/raw_consensus_audit.sh CONSENSUS_DIR REFERENCE_FASTA NEW_OUTPUT_DIR" >&2
    exit 2
fi
input=$1
reference=$2
output=$3
python_bin=${ANVAYA_AUDIT_PYTHON:-python3}
minimap2_bin=${MINIMAP2:-minimap2}
samtools_bin=${SAMTOOLS:-samtools}
"$python_bin" -c 'import pysam'
"$minimap2_bin" --version
"$samtools_bin" --version
for path in "$input/before-contigs.fasta" "$input/quality-only.fasta" "$input/damage-consensus.fasta" "$reference"; do
    [[ -f "$path" ]] || { echo "Missing input: $path" >&2; exit 2; }
done
mkdir -- "$output"
for name in before-contigs quality-only damage-consensus; do
    cp -- "$input/$name.fasta" "$output/$name.fasta"
    "$samtools_bin" faidx "$output/$name.fasta"
done
[[ -f "${reference}.fai" ]] || "$samtools_bin" faidx "$reference"
sha256sum "$output/"*.fasta "$reference" > "$output/input.sha256"
sha256sum experiments/raw_consensus_audit.py experiments/raw_consensus_audit.sh > "$output/audit-source.sha256"
{
    "$minimap2_bin" --version
    "$samtools_bin" --version
    "$python_bin" -c 'import sys, pysam; print(sys.version); print("pysam", pysam.__version__)'
} > "$output/versions.txt"
# Sensitive exploratory settings for short contigs, not a calibrated mapping oracle.
# Keep secondary records; -I 1G avoids a split index for the current 475 Mb reference.
"$minimap2_bin" -a -x sr -k 9 -w 5 -m 20 -s 20 -n 2 -N 20 --secondary=yes -I 1G -t 2 \
    "$reference" "$output/before-contigs.fasta" 2> "$output/mapping.log" |
    "$samtools_bin" sort -n -@ 2 -o "$output/before.name.bam" -
"$python_bin" experiments/raw_consensus_audit.py \
    --before "$output/before-contigs.fasta" \
    --quality "$output/quality-only.fasta" \
    --damage "$output/damage-consensus.fasta" \
    --reference "$reference" --alignments "$output/before.name.bam" \
    --output-dir "$output/audit"
