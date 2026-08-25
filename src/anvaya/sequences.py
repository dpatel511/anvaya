"""DNA sequence utilities."""

DNA_ALPHABET = frozenset("ACGTN")
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def normalize_dna(sequence: str) -> str:
    """Return an uppercase DNA sequence after validating its symbols."""
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    if not sequence:
        raise ValueError("sequence must not be empty")

    normalized = sequence.upper()
    invalid_bases = set(normalized) - DNA_ALPHABET
    if invalid_bases:
        invalid = ", ".join(sorted(invalid_bases))
        raise ValueError(f"sequence contains unsupported symbols: {invalid}")
    return normalized


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return normalize_dna(sequence).translate(_COMPLEMENT)[::-1]


def canonical_sequence(sequence: str) -> str:
    """Return the lexicographically smaller strand representation."""
    normalized = normalize_dna(sequence)
    return min(normalized, reverse_complement(normalized))
