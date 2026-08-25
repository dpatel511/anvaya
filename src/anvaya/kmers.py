"""K-mer utilities."""

from anvaya.sequences import normalize_dna


def kmers(sequence: str, k: int) -> list[str]:
    """Return all overlapping k-mers in *sequence*.

    The sequence is normalized to uppercase. Standard DNA bases and ``N`` are
    accepted; other symbols are rejected.
    """
    if not isinstance(k, int) or isinstance(k, bool):
        raise TypeError("k must be an integer")
    if k < 1:
        raise ValueError("k must be at least 1")

    normalized = normalize_dna(sequence)
    if k > len(normalized):
        raise ValueError("k must not exceed the sequence length")

    return [normalized[start : start + k] for start in range(len(normalized) - k + 1)]
