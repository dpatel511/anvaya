#!/usr/bin/env bash
# Run from the repository root. Dataset experiments are run by the user.
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: bash experiments/support_two_quality_ablation.sh INPUT_FASTQ NEW_OUTPUT_DIR [quality|trimming|fixed]" >&2
    exit 2
fi
input=$1
output=$2
python_bin=${ANVAYA_PYTHON:-python3}
case "${3:-quality}" in
    quality) modes=(neither seed_only recruit_only both) ;;
    trimming) modes=(recruit_only trimmed) ;;
    fixed) modes=(fixed) ;;
    *) echo "Comparison must be quality, trimming or fixed" >&2; exit 2 ;;
esac
[[ -f "$input" ]] || { echo "Input does not exist: $input" >&2; exit 2; }
mkdir -- "$output"

for mode in "${modes[@]}"; do
    case "$mode" in
        neither) flags=(--no-support-two-seed-quality-filter --no-support-two-recruit-quality-filter) ;;
        seed_only) flags=(--support-two-seed-quality-filter --no-support-two-recruit-quality-filter) ;;
        recruit_only) flags=(--no-support-two-seed-quality-filter --support-two-recruit-quality-filter) ;;
        both) flags=(--support-two-seed-quality-filter --support-two-recruit-quality-filter) ;;
        trimmed) flags=(--no-support-two-seed-quality-filter --support-two-recruit-quality-filter --support-two-trim-seed-tails) ;;
        fixed) flags=(--no-support-two-seed-quality-filter --support-two-recruit-quality-filter --support-two-trim-seed-tails --support-two-fixed-membership-before "$output/fixed/before-contigs.fasta") ;;
    esac
    run_dir="$output/$mode"
    mkdir -- "$run_dir"
    echo "Running $mode"
    PYTHONPATH=src /usr/bin/time -v -o "$run_dir/timing.txt" \
        "$python_bin" -m anvaya overlap-assemble -i "$input" \
        --min-cluster-size 5 --max-rounds 3 --min-anchor-matches 1 \
        --min-output-length 31 --max-contig-iterations 0 \
        --ranked-extension --damage-aware-ranking \
        --min-overlap-confidence-margin 0 --reciprocal-best-extension \
        --progressive-raw-phase-audit --max-progressive-raw-iterations 3 \
        --progressive-raw-phase-projection "$run_dir/progressive-contigs.fasta" \
        --selective-support-rescue-audit --adaptive-rescue-min-support 3 \
        --selective-support-rescue-projection "$run_dir/selective-rescue-contigs.fasta" \
        --high-confidence-support-two-rescue-audit \
        --high-confidence-support-two-rescue-projection "$run_dir/support-two-contigs.fasta" \
        --raw-confirmed-master-min-base-quality 20 \
        --raw-confirmed-master-damage-end-window 5 \
        "${flags[@]}" -o "$run_dir/contigs.fasta" \
        > "$run_dir/run.log" 2> "$run_dir/progress.log"
    grep '^overlap_support_two_rescue_' "$run_dir/run.log"
done

if [[ "${3:-quality}" == trimming ]]; then
    for name in contigs.fasta progressive-contigs.fasta selective-rescue-contigs.fasta; do
        cmp "$output/recruit_only/$name" "$output/trimmed/$name"
    done
    echo "Primary, progressive and selective FASTAs match between configurations."
fi
