import csv
import tempfile
import unittest
from pathlib import Path

from anvaya.damage_consensus import apply_damage_aware_consensus
from anvaya.reads import Read


class DamageConsensusTests(unittest.TestCase):
    def test_corrects_terminal_damage_from_nonterminal_molecules(self) -> None:
        truth = "CCGATCGTACGTTAGCTACGATCGTACGTA"
        damaged = "T" + truth[1:]
        reads = [Read("damaged", damaged, (35,) * len(damaged))]
        for index in range(3):
            reads.append(Read(f"truth-{index}", "GGGGG" + truth, (35,) * (len(truth) + 5)))

        corrected, summary = apply_damage_aware_consensus(
            reads,
            end_window=5,
            anchor_k=7,
            minimum_overlap=20,
        )

        self.assertEqual(corrected[0].sequence, truth)
        self.assertEqual(summary.corrected_bases, 1)
        self.assertEqual(summary.corrected_reads, 1)

    def test_abstains_without_independent_support(self) -> None:
        sequence = "TCGATCGTACGTTAGCTACGATCGTACGTA"

        corrected, summary = apply_damage_aware_consensus(
            [Read("single", sequence)],
            end_window=5,
            anchor_k=7,
            minimum_overlap=20,
        )

        self.assertEqual(corrected[0].sequence, sequence)
        self.assertEqual(summary.corrected_bases, 0)
        self.assertGreater(summary.insufficient_support, 0)

    def test_report_records_consensus_evidence(self) -> None:
        truth = "CCGATCGTACGTTAGCTACGATCGTACGTA"
        reads = [Read("damaged", "T" + truth[1:], (35,) * len(truth))]
        reads.extend(
            Read(f"truth-{index}", "GGGGG" + truth, (35,) * (len(truth) + 5))
            for index in range(3)
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "consensus.tsv"

            apply_damage_aware_consensus(
                reads,
                report_path=path,
                end_window=5,
                anchor_k=7,
                minimum_overlap=20,
            )
            with path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        decision = next(row for row in rows if row["read_name"] == "damaged")
        self.assertEqual(decision["decision"], "correct")
        self.assertEqual(decision["proposed_base"], "C")


if __name__ == "__main__":
    unittest.main()
