"""Deterministic, uncalibrated score bands for regulation hypotheses."""

from __future__ import annotations

import math
from typing import Any, cast

import polars as pl
from scipy.stats import beta  # type: ignore[import-untyped]

from gasregnet.schemas import RegulatorCandidatesSchema, validate
from gasregnet.scoring.candidates import CANDIDATE_SCHEMA

SCORE_BAND_MODEL_NAME = "uncalibrated_sigmoid_beta_score_band_94"


def _sigmoid(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-value))


def _score_to_logit_probability(score: float, midpoint: float, scale: float) -> float:
    return max(0.001, min(0.999, _sigmoid((score - midpoint) / scale)))


def assign_operon_regulation_score_bands(
    candidates: pl.DataFrame,
    *,
    band_mass: float = 0.94,
    concentration: float = 24.0,
    midpoint: float = 6.0,
    scale: float = 2.0,
) -> pl.DataFrame:
    """Convert decomposable evidence scores into deterministic, uncalibrated bands."""

    candidates = validate(candidates, RegulatorCandidatesSchema)
    if candidates.is_empty():
        return candidates

    if scale <= 0.0:
        raise ValueError("scale must be positive")

    tail = (1.0 - band_mass) / 2.0
    rows: list[dict[str, Any]] = []
    for row in candidates.iter_rows(named=True):
        updated = dict(row)
        logit_score = _score_to_logit_probability(
            float(cast(float, row["candidate_score"])),
            midpoint,
            scale,
        )
        alpha = logit_score * concentration
        beta_param = (1.0 - logit_score) * concentration
        updated["regulation_logit_score"] = logit_score
        updated["score_band_low"] = float(
            beta.ppf(tail, alpha, beta_param),
        )
        updated["score_band_high"] = float(
            beta.ppf(1.0 - tail, alpha, beta_param),
        )
        updated["score_band_model"] = SCORE_BAND_MODEL_NAME
        rows.append(updated)

    return validate(
        pl.DataFrame(rows, schema_overrides=CANDIDATE_SCHEMA),
        RegulatorCandidatesSchema,
    )


def assign_operon_regulation_posteriors(
    candidates: pl.DataFrame,
    **kwargs: Any,
) -> pl.DataFrame:
    """Backward-compatible wrapper for the renamed score-band layer."""

    return assign_operon_regulation_score_bands(candidates, **kwargs)
