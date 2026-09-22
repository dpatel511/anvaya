"""Audit substitutions against a fixed baseline alignment; never infer strain truth.

Requires pysam, indexed FASTAs and query-name-grouped SAM/BAM. Coordinates are
0-based; reference intervals are half-open. Run separately from the assembler.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter
from contextlib import ExitStack
from itertools import groupby
from pathlib import Path

import pysam

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
COMPARISONS = ((0, 1, "before_to_quality"), (0, 2, "before_to_damage"),
               (1, 2, "quality_to_damage"))


def reference_columns(records, sequence, reference):
    """Return original-query-position -> (reference position, oriented base)."""
    mapped = [record for record in records if not record.is_unmapped]
    if not mapped:
        return "unmapped", {}, "", ""
    if (len(mapped) != 1 or mapped[0].is_secondary or mapped[0].is_supplementary
            or mapped[0].has_tag("SA") or mapped[0].has_tag("XA")):
        return "ambiguous_alignment", {}, "", ""
    record = mapped[0]
    if record.mapping_quality == 255 or record.mapping_quality < 30:
        return "low_or_unknown_mapq", {}, "", ""
    return alignment_columns(record, sequence, reference)


def alignment_columns(record, sequence, reference):
    """Project one full-query CIGAR; secondary SAM records may omit SEQ."""
    if record.is_paired or not record.cigartuples or any(op == 5 for op, _ in record.cigartuples):
        return "unsupported_alignment", {}, "", ""
    expected = sequence.translate(COMPLEMENT)[::-1] if record.is_reverse else sequence
    if record.query_sequence is None and not record.is_secondary:
        raise ValueError(f"{record.query_name}: baseline primary sequence is missing")
    if record.query_sequence is not None and record.query_sequence.upper() != expected:
        raise ValueError(f"{record.query_name}: alignment sequence does not match baseline FASTA")
    if record.infer_query_length() != len(sequence):
        raise ValueError(f"{record.query_name}: alignment query length mismatch")
    name = record.reference_name
    start, stop = record.reference_start, record.reference_end
    if not 0 <= start < stop <= reference.get_reference_length(name):
        raise ValueError(f"{record.query_name}: alignment outside reference")
    # Fetch aligned blocks only: large reference skips must not allocate their span.
    bases = {}
    for block_start, block_stop in record.get_blocks():
        bases.update(enumerate(reference.fetch(name, block_start, block_stop).upper(), block_start))
    columns = {}
    for query_position, reference_position in record.get_aligned_pairs(matches_only=True):
        base = bases[reference_position]
        position = query_position
        if record.is_reverse:
            position = len(sequence) - 1 - position
            base = base.translate(COMPLEMENT)
        columns[position] = (reference_position, base)
    return "aligned", columns, name, "-" if record.is_reverse else "+"


def alternative_sites(records, sequence, reference, positions, secondary_limit):
    """Conditional agreement among reported mappings, not unique mapping or truth.

    All candidates must cover the site. Near-best AS >=95% of best is an
    exploratory sensitivity analysis alongside agreement across all records.
    """
    mapped = [r for r in records if not r.is_unmapped]
    reason = None
    if not mapped:
        reason = "unmapped"
    elif any(r.is_supplementary or r.has_tag("SA") or r.has_tag("XA") for r in mapped):
        reason = "split_or_unexpanded_alternatives"
    elif sum(not r.is_secondary for r in mapped) != 1:
        reason = "missing_or_multiple_primary"
    elif sum(r.is_secondary for r in mapped) >= secondary_limit:
        reason = "secondary_limit_reached"
    elif any(not r.has_tag("AS") for r in mapped):
        reason = "missing_alignment_score"
    elif max(r.get_tag("AS") for r in mapped) <= 0:
        reason = "nonpositive_best_score"
    elif len(mapped) == 1 and (mapped[0].mapping_quality < 30 or mapped[0].mapping_quality == 255):
        reason = "low_or_unknown_mapq"
    policies = {"all_reported": {}, "near_best_95pct": {}}
    if reason:
        return {name: {p: (reason, "", 0) for p in positions} for name in policies}
    projected = [alignment_columns(r, sequence, reference) for r in mapped]
    best = max(r.get_tag("AS") for r in mapped)
    for policy in policies:
        selected = [projection for r, projection in zip(mapped, projected)
                    if policy == "all_reported" or r.get_tag("AS") >= 0.95 * best]
        for position in positions:
            bases = {columns[position][1] for status, columns, _, _ in selected
                     if status == "aligned" and position in columns}
            if any(status != "aligned" for status, _, _, _ in selected):
                status = "unsupported_alignment"
            elif any(position not in columns for _, columns, _, _ in selected):
                status = "alternative_does_not_cover_site"
            elif any(base not in "ACGT" for base in bases):
                status = "ambiguous_reference_base"
            elif len(bases) != 1:
                status = "reference_base_disagreement"
            else:
                status = "reported_reference_agreement"
            policies[policy][position] = (status, next(iter(bases)) if len(bases) == 1 else "", len(selected))
    return policies


def classify(before, after, reference_base):
    if any(base not in "ACGT" for base in (before, after, reference_base)):
        return "ambiguous_base"
    if after == reference_base:
        return "gains_reference_agreement"
    if before == reference_base:
        return "loses_reference_agreement"
    return "neither_matches_reference"


def audit(before, quality, damage, reference, alignments, output_dir, secondary_limit=20):
    if secondary_limit < 1:
        raise ValueError("secondary limit must match a positive mapper -N setting")
    paths = [Path(path) for path in (before, quality, damage, reference)]
    for path in paths:
        if not Path(f"{path}.fai").is_file():
            raise ValueError(f"Missing FASTA index: {path}.fai; run samtools faidx")
    with ExitStack() as stack:
        fastas = [stack.enter_context(pysam.FastaFile(str(path))) for path in paths]
        baseline, _, _, ref = fastas
        lengths = dict(zip(baseline.references, baseline.lengths))
        for path, fasta in zip(paths[:3], fastas[:3]):
            # Validate actual IDs/lengths too: an index alone can hide duplicate IDs.
            seen = set()
            with pysam.FastxFile(str(path)) as records:
                for record in records:
                    if record.name in seen or lengths.get(record.name) != len(record.sequence):
                        raise ValueError(f"{path}: duplicate ID or membership/length mismatch: {record.name}")
                    seen.add(record.name)
            if seen != set(lengths) or dict(zip(fasta.references, fasta.lengths)) != lengths:
                raise ValueError(f"{path}: FASTA membership/index mismatch")
        sam = stack.enter_context(pysam.AlignmentFile(str(alignments)))
        reference_names = set(ref.references)
        for name, length in zip(sam.references, sam.lengths):
            if name not in reference_names or ref.get_reference_length(name) != length:
                raise ValueError(f"Alignment reference header mismatch: {name}")
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=False)
        handle = stack.enter_context((output_dir / "sites.tsv").open("w", newline="", encoding="utf-8"))
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig_id", "position_0based", "comparison", "before", "after",
                         "reference_id", "reference_position_0based", "strand",
                         "reference_base_contig_oriented", "classification"])
        counts = {label: Counter() for _, _, label in COMPARISONS}
        policies = ("all_reported", "near_best_95pct")
        alternative_counts = {policy: {label: Counter() for _, _, label in COMPARISONS} for policy in policies}
        site_handle = stack.enter_context((output_dir / "alternative-sites.tsv").open("w", newline="", encoding="utf-8"))
        site_writer = csv.writer(site_handle, delimiter="\t", lineterminator="\n")
        site_writer.writerow(["contig_id", "position_0based", "comparison", "before", "after",
                              "policy", "selected_alignments", "reference_base_contig_oriented", "classification"])
        alignment_handle = stack.enter_context((output_dir / "alignments.tsv").open("w", newline="", encoding="utf-8"))
        alignment_writer = csv.writer(alignment_handle, delimiter="\t", lineterminator="\n")
        alignment_writer.writerow(["contig_id", "flag", "reference_id", "reference_start_0based",
                                   "reference_stop_exclusive", "strand", "mapq", "AS", "NM", "cigar"])
        statuses = Counter()
        seen = set()

        def evaluate(name, records):
            sequences = [fasta.fetch(name).upper() for fasta in fastas[:3]]
            status, columns, ref_name, strand = reference_columns(records, sequences[0], ref)
            statuses[status] += 1
            positions = {p for p, bases in enumerate(zip(*sequences)) if len(set(bases)) > 1}
            alternatives = alternative_sites(records, sequences[0], ref, positions, secondary_limit) if positions else {}
            for r in records:
                alignment_writer.writerow([name, r.flag, r.reference_name or "", r.reference_start,
                    r.reference_end, "-" if r.is_reverse else "+", r.mapping_quality,
                    r.get_tag("AS") if r.has_tag("AS") else "",
                    r.get_tag("NM") if r.has_tag("NM") else "", r.cigarstring or ""])
            for left, right, label in COMPARISONS:
                for position, (old, new) in enumerate(zip(sequences[left], sequences[right])):
                    if old == new:
                        continue
                    ref_position, base = columns.get(position, ("", ""))
                    result = status if status != "aligned" else (
                        classify(old, new, base) if base else "unaligned_position")
                    counts[label][result] += 1
                    writer.writerow([name, position, label, old, new, ref_name,
                                     ref_position, strand, base, result])
                    for policy in policies:
                        alternative_status, alternative_base, selected = alternatives[policy][position]
                        classification = (classify(old, new, alternative_base)
                            if alternative_status == "reported_reference_agreement" else alternative_status)
                        alternative_counts[policy][label][classification] += 1
                        site_writer.writerow([name, position, label, old, new, policy, selected,
                                              alternative_base, classification])

        for name, group in groupby(sam, key=lambda record: record.query_name):
            if name not in lengths or name in seen:
                raise ValueError(f"Unknown or non-grouped alignment query: {name}")
            seen.add(name)
            evaluate(name, list(group))
        missing = len(lengths) - len(seen)
        for name in lengths:
            if name not in seen:
                evaluate(name, [])
        summary = {"contigs": len(lengths), "alignment_query_groups": len(seen),
                   "contigs_absent_from_alignment": missing,
                   "baseline_alignment_status": dict(statuses), "comparisons": {}}
        def summarize(values):
            gains = values["gains_reference_agreement"]
            losses = values["loses_reference_agreement"]
            resolved = gains + losses + values["neither_matches_reference"]
            return {
                "changed_sites": sum(values.values()), "resolved_sites": resolved,
                "unresolved_sites": sum(values.values()) - resolved,
                "net_reference_agreement_gain": gains - losses,
                "gains_fraction_of_resolved_changes": gains / resolved if resolved else None,
                "classifications": dict(values),
            }
        summary["comparisons"] = {label: summarize(values) for label, values in counts.items()}
        summary["secondary_limit"] = secondary_limit
        summary["audit_source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        summary["pysam_version"] = pysam.__version__
        summary["alternative_comparisons"] = {
            policy: {label: summarize(values) for label, values in comparisons.items()}
            for policy, comparisons in alternative_counts.items()}
        (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("before", "quality", "damage", "reference", "alignments", "output-dir"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--secondary-limit", type=int, default=20,
                        help="Mapper -N limit; reaching it leaves alternative-site evidence unresolved (default 20)")
    args = parser.parse_args()
    print(json.dumps(audit(**vars(args)), indent=2))


if __name__ == "__main__":
    main()
