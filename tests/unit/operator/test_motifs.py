from __future__ import annotations

from pathlib import Path

import pytest

from gasregnet.operator.motifs import (
    PWM,
    load_pwm_csv,
    scan_degenerate_operator,
    scan_pwm,
    score_pwm_window,
)


def test_scan_degenerate_operator_finds_cooa_dyad() -> None:
    hits = scan_degenerate_operator(
        "AAATTGTCAGGGGGGTGACAATTT",
        "TTGTCANNNNNNTGACAA",
        motif_id="CooA_operator",
        scan_reverse=False,
    )

    assert hits.height == 1
    assert hits["motif_id"].item() == "CooA_operator"


def test_scan_pwm_scores_forward_hits(tmp_path: Path) -> None:
    pwm_path = tmp_path / "motif.csv"
    pwm_path.write_text(
        "position,A,C,G,T\n"
        "1,0.97,0.01,0.01,0.01\n"
        "2,0.01,0.97,0.01,0.01\n",
        encoding="utf-8",
    )
    pwm = load_pwm_csv(pwm_path, motif_id="test_pwm")

    hits = scan_pwm("TTACGG", pwm, min_score=3.0, scan_reverse=False)

    assert hits["sequence"].to_list() == ["AC"]


def test_score_pwm_rejects_wrong_width() -> None:
    pwm = PWM(motif_id="x", probabilities=[{"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0}])

    with pytest.raises(ValueError):
        score_pwm_window("", pwm)
