"""Consensus on fixed ungapped raw-read placements, separate from assembly scoring."""

import csv
from dataclasses import dataclass
from functools import lru_cache
from math import exp, isfinite, log
from pathlib import Path
from typing import TYPE_CHECKING, TextIO

from anvaya.reads import Read
from anvaya.mixture_diagnostics import diagnostic_writers, write_diagnostics
from anvaya.linked_allele_guard import blocking_neighbour, variable_neighbours

if TYPE_CHECKING:
    from anvaya.overlap_progressive import ProgressiveSequencePool

BASES = "ACGT"


@dataclass(frozen=True, slots=True)
class RawPlacement:
    """Original read interval [read_start, read_stop), in 0-based coordinates.

    offset places the full oriented read's first base on the contig. Reverse
    placements map original position p to offset + read_length - 1 - p.
    """

    read_index: int
    offset: int
    reverse: bool
    read_start: int
    read_stop: int


@dataclass(frozen=True, slots=True)
class DamageProfile:
    """Untreated double-stranded fragment rates, distance zero at each raw end.

    Only 5' C>T and 3' G>A are modeled; positions beyond each supplied array
    have zero damage probability. Profiles are supplied, not fitted here.
    """

    five_prime_ct: tuple[float, ...]
    three_prime_ga: tuple[float, ...]

    def __post_init__(self):
        for name in ("five_prime_ct", "three_prime_ga"):
            values = tuple(getattr(self, name))
            if any(not isfinite(value) or not 0 <= value <= 1 for value in values):
                raise ValueError("damage probabilities must be finite and in [0, 1]")
            object.__setattr__(self, name, values)

    @classmethod
    def from_prefix(cls, prefix: str | Path) -> "DamageProfile":
        """Read the 12-column CarpeDeam .prof format without silently dropping rates."""
        expected = [f"{a}>{b}" for a in BASES for b in BASES if a != b]
        profiles = []
        for end, transition in (("5p", "C>T"), ("3p", "G>A")):
            path = Path(f"{prefix}{end}.prof")
            values = []
            with path.open(encoding="utf-8") as handle:
                rows = csv.reader(handle, delimiter="\t")
                header = next(rows, [])
                if header != expected:
                    raise ValueError(f"{path}: expected the 12 substitution columns")
                for number, row in enumerate(rows, start=2):
                    if len(row) != 12:
                        raise ValueError(f"{path}:{number}: expected 12 probabilities")
                    rates = [float(value) for value in row]
                    if any(not isfinite(value) or not 0 <= value <= 1 for value in rates):
                        raise ValueError(f"{path}:{number}: invalid probability")
                    if any(value != 0 for label, value in zip(header, rates) if label != transition):
                        raise ValueError(f"{path}:{number}: this model supports only {transition} at this end")
                    values.append(rates[header.index(transition)])
            if not values:
                raise ValueError(f"{path}: profile has no data rows")
            profiles.append(tuple(values))
        return cls(*profiles)


@lru_cache(maxsize=32768)
def _likelihoods(observed: str, quality: int, ct: float, ga: float, reverse: bool) -> tuple[float, ...]:
    error = 10 ** (-quality / 10)
    values = []
    for base in BASES:
        original = BASES[3 - BASES.index(base)] if reverse else base
        damaged, probability = ("T", ct) if original == "C" else (("A", ga) if original == "G" else (original, 0))
        unchanged = 1 - error if original == observed else error / 3
        changed = 1 - error if damaged == observed else error / 3
        values.append((1 - probability) * unchanged + probability * changed)
    return tuple(values)


def observation_likelihoods(observed: str, quality: int, position: int, read_length: int,
                            reverse: bool, profile: DamageProfile) -> tuple[float, ...]:
    """P(raw observed base | each contig-oriented true base), composing damage then error."""
    if observed not in BASES or not 0 <= quality <= 93 or not 0 <= position < read_length:
        raise ValueError("invalid raw base, quality or original read position")
    ct = profile.five_prime_ct[position] if position < len(profile.five_prime_ct) else 0
    distance = read_length - position - 1
    ga = profile.three_prime_ga[distance] if distance < len(profile.three_prime_ga) else 0
    return _likelihoods(observed, quality, ct, ga, reverse)


@dataclass(slots=True)
class RawConsensusDiagnostics:
    contigs: int = 0
    covered_sites: int = 0
    changed_bases: int = 0
    allele_conflict_sites: int = 0
    insufficient_support_sites: int = 0
    posterior_rejections: int = 0
    skipped_quality_observations: int = 0
    ambiguous_base_observations: int = 0
    duplicate_molecule_observations: int = 0
    conflicting_molecule_observations: int = 0
    ambiguous_placement_molecules: int = 0
    linked_allele_rejections: int = 0


def project_raw_consensus(pool: "ProgressiveSequencePool", profile: DamageProfile, *,
                          report: TextIO | None = None,
                          placements_report: TextIO | None = None,
                          mixture_report: TextIO | None = None,
                          linkage_report: TextIO | None = None,
                          linked_allele_guard: bool = False) -> tuple[list[Read], RawConsensusDiagnostics]:
    """Project substitutions only; never alter layout, membership, raw data or length.

    Conditional single-allele model: uniform base prior, independent molecules.
    Require >=3 molecules, >=2 Q20 observations of the proposed base and posterior
    >=0.99. Two alleles each supported by >=2 Q20 observations that are unlikely
    (<0.01) under every other base block correction. This guard is not phasing.
    Missing qualities and Q<3 are uninformative and excluded. Duplicate molecule
    observations agree or the entire molecule is excluded at that position.
    """
    diagnostics = RawConsensusDiagnostics()
    mixture_writer, linkage_writer = diagnostic_writers(mixture_report, linkage_report)
    writer = csv.writer(report, delimiter="\t", lineterminator="\n") if report is not None else None
    placement_writer = csv.writer(placements_report, delimiter="\t", lineterminator="\n") if placements_report is not None else None
    if writer:
        writer.writerow(["contig_id", "position_0based", "before", "after", "model_best", "model_posterior", "molecules", "reason", "guard_neighbour_0based"])
    if placement_writer:
        placement_writer.writerow(["contig_id", "read_index", "read_name", "molecule_id", "strand", "read_start_0based", "read_stop_exclusive", "full_read_offset", "status"])
    output = []
    for index, record in enumerate(pool.active_derived, start=1):
        if record.raw_placements is None:
            raise ValueError("raw placements must be tracked before consensus")
        diagnostics.contigs += 1
        target = record.current.sequence
        columns: list[dict] = [{} for _ in target]
        mappings: dict[int, set[tuple[int, bool]]] = {}
        for placement in record.raw_placements:
            mappings.setdefault(placement.read_index, set()).add((placement.offset, placement.reverse))
        ambiguous_molecules = {pool.records[read_index].molecule_id for read_index, positions in mappings.items()
                               if len(positions) > 1}
        diagnostics.ambiguous_placement_molecules += len(ambiguous_molecules)
        for placement in record.raw_placements:
            source = pool.records[placement.read_index]
            read = source.raw
            length = len(read.sequence)
            if not 0 <= placement.read_start < placement.read_stop <= length:
                raise ValueError("raw placement interval is outside the original read")
            if placement_writer:
                placement_writer.writerow([f"unitig_{index}", source.index, read.name, source.molecule_id,
                    "-" if placement.reverse else "+", placement.read_start, placement.read_stop, placement.offset,
                    "ambiguous_mapping" if source.molecule_id in ambiguous_molecules else "eligible"])
            if source.molecule_id in ambiguous_molecules:
                continue
            for raw_position in range(placement.read_start, placement.read_stop):
                position = placement.offset + (length - 1 - raw_position if placement.reverse else raw_position)
                if not 0 <= position < len(target):
                    continue
                if read.qualities is None or read.qualities[raw_position] < 3:
                    diagnostics.skipped_quality_observations += 1
                    continue
                base = read.sequence[raw_position]
                if base not in BASES:
                    diagnostics.ambiguous_base_observations += 1
                    continue
                quality = read.qualities[raw_position]
                likelihoods = observation_likelihoods(base, quality, raw_position, length, placement.reverse, profile)
                allele = 3 - BASES.index(base) if placement.reverse else BASES.index(base)
                observation = (allele, quality, likelihoods)
                column = columns[position]
                if source.molecule_id in column:
                    diagnostics.duplicate_molecule_observations += 1
                    previous = column[source.molecule_id]
                    if previous is None:
                        continue
                    if previous[0] != allele:
                        column[source.molecule_id] = None
                        diagnostics.conflicting_molecule_observations += 1
                        continue
                    # One observation only: best quality, then least damage ambiguity.
                    rank = lambda item: (item[1], -max(p for b, p in enumerate(item[2]) if b != item[0]))
                    if rank(previous) >= rank(observation):
                        continue
                column[source.molecule_id] = observation
        neighbours = variable_neighbours(columns) if linked_allele_guard else {}
        sequence = list(target)
        for position, column in enumerate(columns):
            observations = [value for value in column.values() if value is not None]
            if not observations:
                continue
            diagnostics.covered_sites += 1
            if len(observations) < 3:
                diagnostics.insufficient_support_sites += 1
                continue
            logs = [0.0] * 4
            counts, robust = [0] * 4, [0] * 4
            for allele, quality, likelihoods in observations:
                for base in range(4):
                    logs[base] += log(max(likelihoods[base], 1e-300))
                if quality >= 20:
                    counts[allele] += 1
                    if max(p for b, p in enumerate(likelihoods) if b != allele) < 0.01:
                        robust[allele] += 1
            best = max(range(4), key=logs.__getitem__)
            posterior = 1 / sum(exp(value - logs[best]) for value in logs)
            reason = "unchanged"
            witness = None
            if sum(count >= 2 for count in robust) >= 2:
                diagnostics.allele_conflict_sites += 1
                reason = "allele_conflict"
            elif BASES[best] != target[position]:
                if counts[best] < 2:
                    diagnostics.insufficient_support_sites += 1
                    reason = "insufficient_observed_support"
                elif posterior < 0.99:
                    diagnostics.posterior_rejections += 1
                    reason = "uncertain"
                else:
                    if linked_allele_guard and target[position] in BASES:
                        witness = blocking_neighbour(position, column, BASES.index(target[position]), best, neighbours)
                    if witness is not None:
                        diagnostics.linked_allele_rejections += 1
                        reason = "linked_allele_conflict"
                    else:
                        sequence[position] = BASES[best]
                        diagnostics.changed_bases += 1
                        reason = "changed"
            if writer and reason != "unchanged":
                writer.writerow([f"unitig_{index}", position, target[position], sequence[position],
                    BASES[best], f"{posterior:.8g}", len(observations), reason,
                    "" if witness is None else witness])
        if mixture_writer is not None or linkage_writer is not None:
            write_diagnostics(columns, f"unitig_{index}", target, sequence, mixture_writer, linkage_writer)
        output.append(Read(record.current.name, "".join(sequence)))
    return output, diagnostics
