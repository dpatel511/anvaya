import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


def _load_module():
    path = (
        Path(__file__).resolve().parents[1]
        / "experiments"
        / "19_overlap_correction_truth_audit.py"
    )
    spec = importlib.util.spec_from_file_location("overlap_correction_truth_audit", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class OverlapCorrectionTruthAuditTests(unittest.TestCase):
    def test_reconstructs_one_edit_per_event(self) -> None:
        module = _load_module()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "events.tsv"
            path.write_text(
                "position\toriginal_base\tcorrected_base\tsequence_before\tsequence_after\n"
                "0\tT\tC\tTTA\tCCA\n"
                "1\tT\tC\tTTA\tCCA\n",
                encoding="utf-8",
            )

            events = module._read_events(path)

        self.assertEqual(
            [event["sequence_after"] for event in events],
            ["CTA", "TCA"],
        )

    def test_classifies_edit_distance_improvement_at_stable_locus(self) -> None:
        module = _load_module()
        before = [{
            "target": "ref", "strand": "+", "target_start": 10,
            "target_end": 45, "matches": 34, "block_length": 35,
            "edit_distance": 1,
        }]
        after = [{**before[0], "matches": 35, "edit_distance": 0}]

        self.assertEqual(
            module._classify(before, after),
            ("fixes_mismatch", 1, 0),
        )

    def test_rejects_tied_reference_placements(self) -> None:
        module = _load_module()
        alignment = {
            "strand": "+", "target_start": 10, "target_end": 45,
            "matches": 35, "block_length": 35, "edit_distance": 0,
        }

        self.assertIsNone(module._unique_best([
            {**alignment, "target": "first"},
            {**alignment, "target": "second"},
        ]))

    def test_parses_only_full_length_paf_alignments(self) -> None:
        module = _load_module()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "alignments.paf"
            path.write_text(
                "full\t35\t0\t35\t+\tref\t100\t10\t45\t34\t35\t60\tNM:i:1\n"
                "partial\t35\t2\t35\t+\tref\t100\t10\t43\t33\t33\t60\tNM:i:0\n",
                encoding="utf-8",
            )

            parsed = module._parse_paf(path)

        self.assertEqual(list(parsed), ["full"])
        self.assertEqual(parsed["full"][0]["edit_distance"], 1)


if __name__ == "__main__":
    unittest.main()
