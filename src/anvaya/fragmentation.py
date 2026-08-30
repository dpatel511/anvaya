"""Report-only attribution of compacted-unitig boundaries."""

import csv
from dataclasses import dataclass
from pathlib import Path

from anvaya.unitig_graph import CompactedUnitigGraph


@dataclass(slots=True, frozen=True)
class FragmentationSummary:
    """Counts of physical unitig ends by compacted-graph topology."""

    ends: int = 0
    dead_ends: int = 0
    unique_continuations: int = 0
    ambiguous_branches: int = 0
    isolated_unitigs: int = 0
    one_sided_unitigs: int = 0
    connected_unitigs: int = 0


def _category(out_degree: int) -> str:
    if out_degree == 0:
        return "dead_end"
    if out_degree == 1:
        return "unique_continuation"
    return "ambiguous_branch"


def _unitig_topology(left_degree: int, right_degree: int) -> str:
    if left_degree == 0 and right_degree == 0:
        return "isolated"
    if left_degree == 0 or right_degree == 0:
        return "one_sided"
    return "connected"


def write_fragmentation_report(
    graph: CompactedUnitigGraph, path: Path
) -> FragmentationSummary:
    """Write one deterministic TSV row per physical unitig end.

    The report attributes boundaries using retained graph topology only. It does
    not infer whether absent edges were caused by filtering, damage, or coverage.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = {"dead_end": 0, "unique_continuation": 0, "ambiguous_branch": 0}
    topology_counts = {"isolated": 0, "one_sided": 0, "connected": 0}
    fieldnames = [
        "unitig_id",
        "side",
        "oriented_handle",
        "category",
        "unitig_topology",
        "out_degree",
        "candidate_handles",
        "length",
        "edge_count",
        "minimum_edge_support",
        "mean_edge_support",
        "terminal_observations",
        "internal_observations",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for unitig_id, sequence in enumerate(graph.sequences):
            ends = (
                ("left", (unitig_id << 1) | 1),
                ("right", unitig_id << 1),
            )
            candidates_by_end = tuple(
                graph.outgoing(oriented_handle) for _, oriented_handle in ends
            )
            topology = _unitig_topology(
                len(candidates_by_end[0]), len(candidates_by_end[1])
            )
            topology_counts[topology] += 1
            for (side, oriented_handle), candidates in zip(
                ends, candidates_by_end, strict=True
            ):
                category = _category(len(candidates))
                counts[category] += 1
                support = graph.oriented_support(oriented_handle)
                writer.writerow(
                    {
                        "unitig_id": unitig_id + 1,
                        "side": side,
                        "oriented_handle": oriented_handle,
                        "category": category,
                        "unitig_topology": topology,
                        "out_degree": len(candidates),
                        "candidate_handles": ",".join(map(str, candidates)),
                        "length": len(sequence),
                        "edge_count": support.edge_count,
                        "minimum_edge_support": support.minimum_edge_support,
                        "mean_edge_support": f"{support.mean_edge_support:.6f}",
                        "terminal_observations": support.terminal_observations,
                        "internal_observations": support.internal_observations,
                    }
                )
    return FragmentationSummary(
        ends=graph.unitig_count * 2,
        dead_ends=counts["dead_end"],
        unique_continuations=counts["unique_continuation"],
        ambiguous_branches=counts["ambiguous_branch"],
        isolated_unitigs=topology_counts["isolated"],
        one_sided_unitigs=topology_counts["one_sided"],
        connected_unitigs=topology_counts["connected"],
    )
