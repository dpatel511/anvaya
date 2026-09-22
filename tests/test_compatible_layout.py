import importlib.util
import sys
import unittest
from dataclasses import replace
from pathlib import Path

from anvaya.overlap_progressive import ProgressiveSequencePool
from anvaya.raw_consensus import RawPlacement
from anvaya.reads import Read
from anvaya.sequences import reverse_complement

path = Path(__file__).resolve().parents[1] / "experiments/compatible_layout.py"
spec = importlib.util.spec_from_file_location("compatible_layout", path)
layout = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = layout
spec.loader.exec_module(layout)


class CompatibleLayoutTests(unittest.TestCase):
    def pool(self, reverse=False, conflict=False, duplicate=False):
        sequence = "ACGTTGCATGACCTGATCGTACCGTAGCTAGGCTACGTAGTCGATCGATGCTAGCTACGATCGGTACCTGATCGTAGCTAGCA"
        child = sequence[10:70]
        if reverse:
            child = reverse_complement(child)
        reads = [Read("parent", sequence), Read("child", child), Read("support", sequence[20:50])]
        pool = ProgressiveSequencePool.from_reads(reads, [0, 1, 0 if duplicate else 2])
        parent = replace(pool.records[0].corrected(reads[0]), raw_placements=(
            RawPlacement(0, 0, False, 0, len(sequence)), RawPlacement(2, 20, False, 0, 30)))
        placements = [RawPlacement(0, -10, False, 10, 70), RawPlacement(2, 10, False, 0, 30)]
        if reverse:
            placements = [replace(p, offset=60 - p.offset - len(reads[p.read_index].sequence), reverse=True) for p in placements]
        if conflict:
            placements[1] = replace(placements[1], offset=placements[1].offset + 1)
        child_record = replace(pool.records[1].corrected(reads[1]), raw_placements=tuple(placements))
        return pool.replace_record(parent).replace_record(child_record)

    def test_containment_both_strands_preserves_raw(self):
        for reverse in (False, True):
            pool = self.pool(reverse=reverse)
            result, retired = layout.consolidate(pool)
            self.assertEqual(retired, 1)
            self.assertEqual(len(result.active_derived), 1)
            self.assertEqual(result.raw_evidence, pool.raw_evidence)
            self.assertEqual({(p.read_index, p.offset, p.reverse) for p in result.records[0].raw_placements},
                             {(0, 0, False), (2, 20, False)})

    def test_conflict_and_duplicate_molecules_block_retirement(self):
        for pool in (self.pool(conflict=True), self.pool(duplicate=True)):
            result, retired = layout.consolidate(pool)
            self.assertEqual(retired, 0)
            self.assertEqual(result, pool)
