"""Unitig extraction from de Bruijn graphs."""

from anvaya.graph import DeBruijnGraph, node_degrees


def _spell_path(path: list[str]) -> str:
    return path[0] + "".join(node[-1] for node in path[1:])


def extract_unitigs(graph: DeBruijnGraph) -> list[str]:
    """Return maximal non-branching sequences from *graph*.

    Each distinct edge is used once. Edge multiplicity does not create
    duplicate unitigs.
    """
    visited: set[tuple[str, str]] = set()
    unitigs: list[str] = []
    in_degrees, out_degrees = node_degrees(graph)

    for start in sorted(graph):
        if in_degrees[start] == 1 and out_degrees[start] == 1:
            continue

        for successor in sorted(graph[start]):
            edge = (start, successor)
            if edge in visited:
                continue

            path = [start, successor]
            visited.add(edge)
            current = successor

            while in_degrees[current] == 1 and out_degrees[current] == 1:
                following = next(iter(graph[current]))
                edge = (current, following)
                if edge in visited:
                    break
                path.append(following)
                visited.add(edge)
                current = following

            unitigs.append(_spell_path(path))

    for start in sorted(graph):
        for successor in sorted(graph[start]):
            edge = (start, successor)
            if edge in visited:
                continue

            path = [start, successor]
            visited.add(edge)
            current = successor

            while current != start:
                following = next(iter(graph[current]))
                edge = (current, following)
                if edge in visited:
                    break
                path.append(following)
                visited.add(edge)
                current = following

            unitigs.append(_spell_path(path))

    return unitigs
