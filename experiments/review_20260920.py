"""Bounded reviewer probes; synthetic inputs only, no assembly changes.

Coordinates are 0-based, half-open. These reproduce limitations, not validation.
"""
import hashlib
import json
import random
from pathlib import Path

from anvaya.damage_string_graph import assemble
from anvaya.overlap_assembly import _MasterOverlapEdge, _n50
from anvaya.paired_reads import merge_overlapping_pairs
from anvaya.reads import Read
from anvaya.sequences import reverse_complement
from paired_branch_opportunity import audit_paired_transitions


def merge_probe():
    rng = random.Random(199001)
    reference = ''.join(rng.choice('ACGT') for _ in range(2400))
    placement = random.Random(199002)
    left, right, fragments = [], [], []
    for i in range(round(20 * 2400 / 135)):
        length = placement.randint(90, 180)
        start = placement.randrange(2400 - length + 1)
        fragment = reference[start:start + length]
        fragments.append(fragment)
        left.append(Read(f'pair_{i}/1', fragment[:75], (35,) * 75))
        right.append(Read(f'pair_{i}/2', reverse_complement(fragment[-75:]), (35,) * 75))
    merged, molecules, diagnostics = merge_overlapping_pairs(left, right)
    incorrect = [read.name for read in merged if '|merged' in read.name
                 and read.sequence != fragments[int(read.name.split('_')[1].split('/')[0])]]
    result = {'merged_pairs': diagnostics.merged_pairs, 'incorrect_merged_sequences': incorrect}
    raw = [read for pair in zip(left, right) for read in pair]
    for label, reads, ids in [('unmerged', raw, [i for i in range(len(left)) for _ in range(2)]),
                              ('merged', merged, molecules)]:
        captured = {}
        def capture(nodes, edges, excluded):
            captured['reduced_directed_edges'] = len(edges)
            captured['full_containment_directed_edges'] = sum(
                edge.overlap >= min(len(nodes[edge.source[0]].sequence), len(nodes[edge.target[0]].sequence))
                for edge in edges.values())
        pool, graph = assemble(reads, molecule_ids=ids, maximum_reads=len(reads), reduced_graph_audit=capture)
        result[label] = dict(captured, n50=_n50([len(r.current.sequence) for r in pool.active_derived]),
                             contigs=len(pool.active_derived), ambiguous_containments=graph['ambiguous_containments'])
    return result


def paired_conflict_probe():
    rng = random.Random(91)
    reference = ''.join(rng.choice('ACGT') for _ in range(140))
    nodes, edges = [], {}
    def branch(shift):
        base = len(nodes)
        source, centre, target = reference[:55 + shift], reference[35 + shift:80 + shift], reference[60 + shift:115 + shift]
        nodes.extend([Read('s', source), Read('c', centre), Read('t', target),
                      Read('alt', reference[60 + shift:80 + shift] + 'A' * 35)])
        for edge in [_MasterOverlapEdge((base, False), (base+1, False), 35+shift, 20, ()),
                     _MasterOverlapEdge((base+1, False), (base+2, False), 25, 20, ()),
                     _MasterOverlapEdge((base+1, False), (base+3, False), 25, 20, ())]:
            edges[(edge.source, edge.target)] = edge
    left = [Read(f'p{i}/1', reference[10:40]) for i in range(2)]
    right = [Read(f'p{i}/2', reverse_complement(reference[80:110])) for i in range(2)]
    branch(0)
    single = audit_paired_transitions(nodes, edges, set(), left, right)
    branch(5)
    both = audit_paired_transitions(nodes, edges, set(), left, right)
    return {'one_branch': single, 'two_consistent_branch_locations': both}


if __name__ == '__main__':
    result = {'merge': merge_probe(), 'global_pair_conflict': paired_conflict_probe()}
    assert result['merge']['merged_pairs'] == 138
    assert not result['merge']['incorrect_merged_sequences']
    assert result['merge']['unmerged']['n50'] == 2399
    assert result['merge']['merged']['n50'] == 99
    assert result['global_pair_conflict']['one_branch']['uniquely_resolved_transitions'] == 1
    assert result['global_pair_conflict']['two_consistent_branch_locations']['conflicting_molecules'] == 2
    result['source_sha256'] = {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in [Path(__file__), *Path('src/anvaya').glob('*.py'),
                     Path('experiments/paired_branch_opportunity.py'), Path('experiments/branch_transition_opportunity.py')]
    }
    destination = Path('results/reviewer-20260920')
    destination.mkdir(exist_ok=True, parents=True)
    (destination / 'probes.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'source_sha256'}, indent=2))
