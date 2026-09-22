import argparse
import gzip
import importlib.util
import json
import random
import tempfile
import unittest
from pathlib import Path

import pysam

_path = Path(__file__).resolve().parents[1] / "experiments/overlap_opportunity.py"
_spec = importlib.util.spec_from_file_location("overlap_opportunity", _path)
audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(audit)


class OpportunityTests(unittest.TestCase):
    def test_fastq_name_filter_reads_gzip_and_rejects_duplicates(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fq.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("@read1 description\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n")
            self.assertEqual(audit.fastq_names(path), {"read1", "read2"})
            with gzip.open(path, "at") as handle:
                handle.write("@read1\nACGT\n+\nIIII\n")
            with self.assertRaisesRegex(ValueError, "unique"):
                audit.fastq_names(path)

    def test_union_and_overlap_threshold(self):
        self.assertEqual(audit.interval_metrics([(0, 60), (30, 90), (61, 110)]), (110, [90, 49]))
        self.assertEqual(audit.interval_metrics([(0, 60), (30, 90), (60, 110)]), (110, [110]))
        self.assertEqual(audit.interval_metrics([]), (0, []))

    def test_mapping_filters(self):
        record = pysam.AlignedSegment()
        record.query_name = "r"
        record.query_sequence = "A" * 50
        record.reference_id = 0
        record.reference_start = 0
        record.cigarstring = "50M"
        record.mapping_quality = 30
        self.assertEqual(audit.eligible([record])[0], "eligible")
        self.assertEqual(audit.eligible([record, record])[0], "non_single_primary")
        record.mapping_quality = 255
        self.assertEqual(audit.eligible([record])[0], "low_or_unknown_mapq")
        record.mapping_quality = 30
        record.cigarstring = "1S49M"
        self.assertEqual(audit.eligible([record])[0], "clipped_or_gapped")

    def test_end_to_end_known_overlap_both_strands(self):
        rng = random.Random(92008)
        truth = "".join(rng.choice("ACGT") for _ in range(120))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "reads.bam"
            with pysam.AlignmentFile(str(bam), "wb", header={"HD": {"SO": "queryname"},
                                      "SQ": [{"SN": "ref", "LN": 120}]}) as handle:
                for index, start in enumerate((0, 30, 40)):
                    record = pysam.AlignedSegment()
                    record.query_name = f"read{index}"
                    record.query_sequence = truth[start:start + 80]
                    record.flag = 16 if index else 0
                    record.reference_id = 0
                    record.reference_start = start
                    record.mapping_quality = 60
                    record.cigarstring = "80M"
                    record.query_qualities = [35] * 80
                    handle.write(record)
            cohort = root / "cohort.fq.gz"
            with gzip.open(cohort, "wt") as handle:
                handle.write(f"@read0\n{truth[:80]}\n+\n{'I' * 80}\n")
                handle.write(f"@read1\n{truth[30:110]}\n+\n{'I' * 80}\n")
            audit.run(argparse.Namespace(alignments=bam, output=root / "out",
                      max_reads=10, targets=10, read_names_fastq=cohort))
            summary = json.loads((root / "out/summary.json").read_text())
            self.assertEqual(summary["query_groups"], 2)
            self.assertEqual(summary["requested_query_groups"], 2)
            self.assertEqual(summary["source_query_groups"], 3)
            self.assertEqual(summary["directed_expected_overlap_stages"], {"selected": 2})
            self.assertEqual(summary["selected_dovetails_outside_proxy"], 0)
