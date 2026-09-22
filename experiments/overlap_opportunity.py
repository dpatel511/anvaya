"""Reference-based overlap opportunity diagnostic, not simulation ground truth.

Input is query-name grouped SAM/BAM. Internal coordinates are 0-based half-open.
Only single, MAPQ >=30, full-length ungapped mappings enter the proxy benchmark.
"""

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from anvaya.overlap_assembly import _candidate_alignments
from anvaya.overlap_index import _anchor_index
from anvaya.reads import Read


def fastq_names(path):
    """Return unique FASTQ query names without retaining sequences."""
    opener = gzip.open if path.suffix == ".gz" else open
    names = set()
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if not sequence or not separator or not quality:
                raise ValueError("truncated FASTQ record")
            if not header.startswith("@") or not separator.startswith("+"):
                raise ValueError("invalid FASTQ record")
            if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
                raise ValueError("FASTQ sequence and quality lengths differ")
            name = header[1:].split()[0]
            if not name or name in names:
                raise ValueError("FASTQ query names must be nonempty and unique")
            names.add(name)
    return names


def interval_metrics(intervals, minimum_overlap=30):
    """Coverage union and components connected by overlaps of at least the cutoff."""
    covered = 0
    union_end = -1
    components = []
    component_start = component_end = None
    for start, end in sorted(intervals):
        covered += max(0, end - max(start, union_end))
        union_end = max(union_end, end)
        if component_end is None or min(component_end, end) - start < minimum_overlap:
            if component_end is not None:
                components.append(component_end - component_start)
            component_start, component_end = start, end
        else:
            component_end = max(component_end, end)
    if component_end is not None:
        components.append(component_end - component_start)
    return covered, components


def eligible(group):
    primary = [r for r in group if not r.is_secondary and not r.is_supplementary]
    if len(primary) != 1:
        return "non_single_primary", None
    record = primary[0]
    if record.is_unmapped:
        return "unmapped", None
    if len(group) != 1 or record.has_tag("SA") or record.has_tag("XA"):
        return "alternative_or_split", None
    if record.mapping_quality == 255 or record.mapping_quality < 30:
        return "low_or_unknown_mapq", None
    if not record.query_sequence or not record.cigartuples or any(
        op not in (0, 7, 8) for op, length in record.cigartuples
    ):
        return "clipped_or_gapped", None
    if record.reference_length != len(record.query_sequence):
        return "length_mismatch", None
    return "eligible", record


def run(args):
    if args.max_reads < 1 or args.targets < 1:
        raise ValueError("read and target bounds must be positive")
    requested_names = fastq_names(args.read_names_fastq) if getattr(
        args, "read_names_fastq", None
    ) else None
    statuses = Counter()
    loci = defaultdict(list)
    reads, placements = [], []
    seen = set()
    selected_seen = set()
    with pysam.AlignmentFile(str(args.alignments)) as handle:
        lengths = dict(zip(handle.references, handle.lengths))
        for name, records in itertools.groupby(handle, key=lambda r: r.query_name):
            if name in seen:
                raise ValueError("input must be query-name grouped with unique read names")
            seen.add(name)
            records = list(records)
            if requested_names is not None and name not in requested_names:
                continue
            selected_seen.add(name)
            if len(selected_seen) > args.max_reads:
                raise ValueError("input exceeds --max-reads; use a bounded dataset")
            status, record = eligible(records)
            statuses[status] += 1
            if record is None:
                continue
            index = len(reads)
            # SAM SEQ and QUAL already follow the reference orientation on reverse hits.
            reads.append(Read(name, record.query_sequence,
                              tuple(record.query_qualities) if record.query_qualities is not None else None))
            placement = (record.reference_name, record.reference_start, record.reference_end)
            placements.append(placement)
            loci[placement[0]].append((placement[1], placement[2], index))
    if requested_names is not None:
        missing = requested_names - selected_seen
        if missing:
            raise ValueError(f"{len(missing)} FASTQ query names are absent from alignments")
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "coverage.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["reference", "reference_length", "eligible_reads", "aligned_bases",
                         "covered_bases", "mean_depth", "overlap_components", "longest_component"])
        for reference, length in lengths.items():
            intervals = [(s, e) for s, e, _ in loci[reference]]
            covered, components = interval_metrics(intervals)
            bases = sum(e - s for s, e in intervals)
            writer.writerow([reference, length, len(intervals), bases, covered,
                             bases / length if length else 0, len(components), max(components, default=0)])
    stage_counts = Counter()
    selected_outside_proxy = sampled_targets = 0
    if reads:
        anchors = _anchor_index(reads, 15, 0, 8, 100)
        molecule_ids = list(range(len(reads)))
        position_bits = max(len(r.sequence) for r in reads).bit_length()
        window = max(len(r.sequence) for r in reads)
        indices = random.Random(92007).sample(range(len(reads)), min(args.targets, len(reads)))
        with (args.output / "overlaps.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["target", "partner", "reference", "offset_0based", "stage"])
            for index in indices:
                reference, start, end = placements[index]
                expected = {}
                for other_start, other_end, other in loci[reference]:
                    if other == index:
                        continue
                    overlap = min(end, other_end) - max(start, other_start)
                    # Dovetails only: exclude containment and identical starts/ends.
                    if overlap >= 30 and ((start < other_start and end < other_end)
                                          or (other_start < start and other_end < end)):
                        expected[(other, False, other_start - start)] = other
                trace = {}
                candidates = _candidate_alignments(reads[index].sequence, reads, molecule_ids,
                    anchors, {index}, anchor_k=15, anchors_per_read=8, maximum_anchor_occurrences=100,
                    minimum_anchor_matches=1, minimum_overlap=30, minimum_identity=.9,
                    minimum_ry_identity=.99, position_bits=position_bits, target_window=window,
                    stage_trace=trace)
                sampled_targets += 1
                for key, other in expected.items():
                    stage = trace.get(key, "no_anchor_vote")
                    stage_counts[stage] += 1
                    writer.writerow([reads[index].name, reads[other].name, reference, key[2], stage])
                for candidate in candidates:
                    # Only count dovetails; outside-proxy placements are not proven false.
                    if (candidate.offset < 0 and candidate.offset + len(candidate.sequence) < len(reads[index].sequence)) or (
                        candidate.offset > 0 and candidate.offset + len(candidate.sequence) > len(reads[index].sequence)):
                        selected_outside_proxy += (candidate.read_index, candidate.reverse, candidate.offset) not in expected
    summary = dict(evidence="reference_mapping_proxy_not_truth",
                   source_query_groups=len(seen), query_groups=len(selected_seen),
                   requested_query_groups=(len(requested_names)
                                           if requested_names is not None else None),
                   statuses=dict(statuses), sampled_targets=sampled_targets,
                   directed_expected_overlap_stages=dict(stage_counts),
                   selected_dovetails_outside_proxy=selected_outside_proxy,
                   settings=dict(anchor_k=15, anchors_per_read=8, occurrence_cap=100,
                                 minimum_overlap=30, minimum_identity=.9, minimum_ry_identity=.99,
                                 ry_rescue=False, seed=92007),
                   limitations=["Only single MAPQ>=30 full ungapped alignments; selection bias toward easy reads",
                                "Reference mapping is not read-origin truth; unresolved reads excluded",
                                "Candidate index uses eligible reads only; occurrence caps differ from full assembly",
                                "Components assume correct mappings and ignore repeats, strain ambiguity and support gates",
                                "Directed pairs are correlated; counts are not independent statistical replicates",
                                "No clustering, ownership, extension or consensus acceptance measured"])
    root = Path(__file__).resolve().parents[1]
    summary["source_sha256"] = {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (Path(__file__), root / "src/anvaya/overlap_assembly.py", root / "src/anvaya/overlap_index.py")}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignments", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-reads", type=int, default=100000)
    parser.add_argument("--targets", type=int, default=2000)
    parser.add_argument("--read-names-fastq", type=Path,
                        help="restrict the audit to query names in this FASTQ[.gz]")
    run(parser.parse_args())
