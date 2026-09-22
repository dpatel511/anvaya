"""Bounded synthetic-only exact containment consolidation; 0-based coordinates."""
from dataclasses import replace

from anvaya.sequences import reverse_complement


def consolidate(pool):
    """Retire exact contained centers only with two shared, consistent molecules.

    This quadratic prototype is intentionally limited to small synthetic pools.
    It never joins non-contained sequences or uses origin truth.
    """
    if len(pool.active_derived) > 1000:
        raise ValueError("synthetic consolidation is limited to 1000 contigs")
    order = sorted(pool.active_derived, key=lambda r: (-len(r.current.sequence), r.index))
    retired = 0
    for number, initial in enumerate(order):
        keeper = pool.records[initial.index]
        if keeper not in pool.active_derived or not keeper.raw_placements:
            continue
        for initial_child in order[number + 1:]:
            child = pool.records[initial_child.index]
            if child not in pool.active_derived or not child.raw_placements:
                continue
            matches = []
            for reverse in (False, True):
                sequence = reverse_complement(child.current.sequence) if reverse else child.current.sequence
                offset = keeper.current.sequence.find(sequence)
                if offset >= 0:
                    matches.append((offset, reverse))
                    if keeper.current.sequence.find(sequence, offset + 1) >= 0:
                        matches.append((offset + 1, reverse))  # Mark repeated containment ambiguous.
            if len(matches) != 1:
                continue
            offset, reverse = matches[0]
            translated = []
            for p in child.raw_placements:
                raw_length = len(pool.records[p.read_index].raw.sequence)
                translated.append(replace(p,
                    offset=offset + (len(child.current.sequence) - p.offset - raw_length if reverse else p.offset),
                    reverse=p.reverse != reverse))
            def locations(placements):
                result = {}
                for p in placements:
                    result.setdefault(p.read_index, set()).add((p.offset, p.reverse))
                return result
            first, second = locations(keeper.raw_placements), locations(translated)
            shared = first.keys() & second.keys()
            molecules = {pool.records[index].molecule_id for index in shared}
            if len(molecules) < 2 or any(len(values) != 1 for values in [*first.values(), *second.values()]):
                continue
            if any(first[index] != second[index] for index in shared):
                continue
            keeper = replace(keeper,
                raw_placements=tuple(dict.fromkeys((*keeper.raw_placements, *translated))),
                contributing_molecules=keeper.contributing_molecules | child.contributing_molecules)
            pool = pool.replace_record(keeper)
            pool = pool.replace_record(child.consumed())
            retired += 1
    return pool, retired
