"""K-mer utilities."""

DNA_ALPHABET = frozenset("ACGTN")


def kmers(sequence: str, k: int) -> list[str]:
    """Return all overlapping k-mers in *sequence*.

    The sequence is normalized to uppercase. Standard DNA bases and ``N`` are
    accepted; other symbols are rejected.
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    if not isinstance(k, int) or isinstance(k, bool):
        raise TypeError("k must be an integer")
    if not sequence:
        raise ValueError("sequence must not be empty")
    if k < 1:
        raise ValueError("k must be at least 1")

    normalized = sequence.upper()
    invalid_bases = set(normalized) - DNA_ALPHABET
    if invalid_bases:
        invalid = ", ".join(sorted(invalid_bases))
        raise ValueError(f"sequence contains unsupported symbols: {invalid}")
    if k > len(normalized):
        raise ValueError("k must not exceed the sequence length")

    return [normalized[start : start + k] for start in range(len(normalized) - k + 1)]

