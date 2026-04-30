"""PWM and degenerate operator scanning."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import polars as pl

BASES = ("A", "C", "G", "T")
IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}
COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


@dataclass(frozen=True)
class PWM:
    """Position weight matrix with log-odds scoring."""

    motif_id: str
    probabilities: list[dict[str, float]]

    @property
    def width(self) -> int:
        return len(self.probabilities)


def reverse_complement(sequence: str) -> str:
    """Return reverse complement for DNA sequence text."""

    return sequence.translate(COMPLEMENT)[::-1].upper()


def load_pwm_csv(path: Path, *, motif_id: str | None = None) -> PWM:
    """Load a PWM CSV with columns position,A,C,G,T."""

    frame = pl.read_csv(path)
    missing = set(BASES) - set(frame.columns)
    if missing:
        raise ValueError(f"PWM missing base column(s): {', '.join(sorted(missing))}")
    probabilities: list[dict[str, float]] = []
    for row in frame.sort("position").iter_rows(named=True):
        values = {base: float(row[base]) for base in BASES}
        total = sum(values.values())
        if total <= 0.0:
            raise ValueError("PWM rows must have positive mass")
        probabilities.append({base: value / total for base, value in values.items()})
    return PWM(motif_id=motif_id or path.stem, probabilities=probabilities)


def score_pwm_window(window: str, pwm: PWM, *, background: float = 0.25) -> float:
    """Return log2 odds score for one DNA window."""

    if len(window) != pwm.width:
        raise ValueError("window length must match PWM width")
    score = 0.0
    for base, probabilities in zip(window.upper(), pwm.probabilities, strict=True):
        score += math.log2(max(probabilities.get(base, 0.0), 1e-6) / background)
    return score


def scan_pwm(
    sequence: str,
    pwm: PWM,
    *,
    min_score: float,
    scan_reverse: bool = True,
) -> pl.DataFrame:
    """Scan a DNA sequence for PWM hits on one or both strands."""

    rows: list[dict[str, object]] = []
    sequence = sequence.upper()
    strands = [("+", sequence)]
    if scan_reverse:
        strands.append(("-", reverse_complement(sequence)))
    for strand, strand_sequence in strands:
        for start in range(0, len(strand_sequence) - pwm.width + 1):
            window = strand_sequence[start : start + pwm.width]
            if re.search(r"[^ACGT]", window):
                continue
            score = score_pwm_window(window, pwm)
            if score >= min_score:
                rows.append(
                    {
                        "motif_id": pwm.motif_id,
                        "start": start,
                        "stop": start + pwm.width,
                        "strand": strand,
                        "sequence": window,
                        "score": score,
                    },
                )
    return pl.DataFrame(
        rows,
        schema={
            "motif_id": pl.Utf8,
            "start": pl.Int64,
            "stop": pl.Int64,
            "strand": pl.Utf8,
            "sequence": pl.Utf8,
            "score": pl.Float64,
        },
    )


def degenerate_to_regex(pattern: str) -> re.Pattern[str]:
    """Compile an IUPAC degenerate operator pattern into a regex."""

    pieces = []
    for character in pattern.upper().replace("-", ""):
        alphabet = IUPAC.get(character)
        if alphabet is None:
            raise ValueError(f"unsupported IUPAC base: {character}")
        pieces.append(f"[{alphabet}]")
    return re.compile("".join(pieces))


def scan_degenerate_operator(
    sequence: str,
    pattern: str,
    *,
    motif_id: str,
    scan_reverse: bool = True,
) -> pl.DataFrame:
    """Scan for an exact IUPAC-degenerate operator such as CooA dyads."""

    rows: list[dict[str, object]] = []
    regex = degenerate_to_regex(pattern)
    for strand, strand_sequence in [("+", sequence.upper())] + (
        [("-", reverse_complement(sequence))] if scan_reverse else []
    ):
        for match in regex.finditer(strand_sequence):
            rows.append(
                {
                    "motif_id": motif_id,
                    "start": match.start(),
                    "stop": match.end(),
                    "strand": strand,
                    "sequence": match.group(0),
                    "score": float(match.end() - match.start()),
                },
            )
    return pl.DataFrame(
        rows,
        schema={
            "motif_id": pl.Utf8,
            "start": pl.Int64,
            "stop": pl.Int64,
            "strand": pl.Utf8,
            "sequence": pl.Utf8,
            "score": pl.Float64,
        },
    )
