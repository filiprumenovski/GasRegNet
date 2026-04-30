"""Operon-level posterior probabilities for regulation hypotheses."""

from __future__ import annotations

import math
from typing import Any, cast

import polars as pl
from scipy.stats import beta  # type: ignore[import-untyped]

from gasregnet.schemas import RegulatorCandidatesSchema, validate
from gasregnet.scoring.candidates import CANDIDATE_SCHEMA

POSTERIOR_MODEL_NAME = "baseline_logit_beta_hdi_94"


def _sigmoid(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-value))


def _score_to_posterior(score: float, midpoint: float, scale: float) -> float:
    return max(0.001, min(0.999, _sigmoid((score - midpoint) / scale)))


def assign_operon_regulation_posteriors(
    candidates: pl.DataFrame,
    *,
    hdi_mass: float = 0.94,
    concentration: float = 24.0,
    midpoint: float = 6.0,
    scale: float = 2.0,
) -> pl.DataFrame:
    """Convert decomposable evidence scores into posterior probabilities.

    This is a deterministic baseline posterior layer: the raw candidate score remains
    available, but report-facing outputs can use ``P(regulation | evidence)`` plus a
    beta-approximate HDI. Synthetic-truth calibration can later replace the fixed
    score-to-logit mapping with a fitted calibration curve without changing schemas.
    """

    candidates = validate(candidates, RegulatorCandidatesSchema)
    if candidates.is_empty():
        return candidates

    if scale <= 0.0:
        raise ValueError("scale must be positive")

    tail = (1.0 - hdi_mass) / 2.0
    rows: list[dict[str, Any]] = []
    for row in candidates.iter_rows(named=True):
        updated = dict(row)
        posterior = _score_to_posterior(
            float(cast(float, row["candidate_score"])),
            midpoint,
            scale,
        )
        alpha = posterior * concentration
        beta_param = (1.0 - posterior) * concentration
        updated["regulation_posterior"] = posterior
        updated["regulation_posterior_hdi_low"] = float(
            beta.ppf(tail, alpha, beta_param),
        )
        updated["regulation_posterior_hdi_high"] = float(
            beta.ppf(1.0 - tail, alpha, beta_param),
        )
        updated["posterior_evidence_model"] = POSTERIOR_MODEL_NAME
        rows.append(updated)

    return validate(
        pl.DataFrame(rows, schema_overrides=CANDIDATE_SCHEMA),
        RegulatorCandidatesSchema,
    )
