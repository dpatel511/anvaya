#!/usr/bin/env bash
# Run from the repository root. The user runs actual FASTQ experiments.
set -euo pipefail
if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: bash experiments/raw_consensus_comparison.sh MERGED_FASTQ PROFILE_PREFIX NEW_OUTPUT_DIR [baseline|guarded|ry-rescue|ry-rescue-audit|ry-rescue-links]" >&2
    exit 2
fi
input=$1
profile=$2
output=$3
python_bin=${ANVAYA_PYTHON:-python3}
guard_args=()
case ${4:-baseline} in
    baseline) ;;
    guarded) guard_args+=(--raw-consensus-linked-allele-guard) ;;
    ry-rescue) guard_args+=(--raw-consensus-linked-allele-guard --progressive-low-quality-ry-rescue) ;;
    ry-rescue-audit) guard_args+=(--raw-consensus-linked-allele-guard --progressive-low-quality-ry-rescue --progressive-extension-failure-audit) ;;
    ry-rescue-links) guard_args+=(--raw-consensus-linked-allele-guard --progressive-low-quality-ry-rescue --raw-supported-progressive-link-audit --progressive-link-min-read-support 2 --raw-supported-progressive-link-projection "$output/linked-contigs.fasta") ;;
    *) echo "Expected baseline, guarded, ry-rescue, ry-rescue-audit, or ry-rescue-links mode" >&2; exit 2 ;;
esac
for path in "$input" "${profile}5p.prof" "${profile}3p.prof"; do
    [[ -f "$path" ]] || { echo "Missing input: $path" >&2; exit 2; }
done
mkdir -- "$output"
cp -- "${profile}5p.prof" "$output/profile5p.prof"
cp -- "${profile}3p.prof" "$output/profile3p.prof"
sha256sum "$output/profile5p.prof" "$output/profile3p.prof" > "$output/profile.sha256"
git rev-parse HEAD > "$output/revision.txt"
sha256sum src/anvaya/*.py > "$output/source.sha256"
PYTHONPATH=src /usr/bin/time -v -o "$output/timing.txt" \
    "$python_bin" -m anvaya overlap-assemble -i "$input" \
    --min-cluster-size 5 --max-rounds 3 --min-anchor-matches 1 \
    --min-output-length 31 --max-contig-iterations 0 \
    --ranked-extension --damage-aware-ranking \
    --min-overlap-confidence-margin 0 --reciprocal-best-extension \
    --progressive-raw-phase-audit --max-progressive-raw-iterations 3 \
    --progressive-raw-phase-projection "$output/before-contigs.fasta" \
    --raw-consensus-profile-prefix "$output/profile" \
    --raw-consensus-projection "$output/damage-consensus.fasta" \
    --raw-consensus-quality-control "$output/quality-only.fasta" \
    --raw-consensus-report "$output/decisions.tsv" \
    --raw-consensus-placements "$output/placements.tsv" \
    --raw-consensus-mixture-report "$output/mixture.tsv" \
    --raw-consensus-linkage-report "$output/linkage.tsv" \
    "${guard_args[@]}" \
    -o "$output/primary-contigs.fasta" > "$output/run.log" 2> "$output/progress.log"
grep -E '^raw_(consensus|quality_control)_' "$output/run.log"
cat "$output/timing.txt"
