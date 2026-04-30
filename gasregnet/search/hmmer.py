"""HMMER search helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import polars as pl
from pyhmmer import easel, hmmer, plan7

from gasregnet.errors import MissingInputError


def _decode_name(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _query_name(top_hits: Any) -> str:
    query = getattr(top_hits, "query", None)
    if query is not None and hasattr(query, "name"):
        return _decode_name(query.name)
    return _decode_name(getattr(top_hits, "query_name", ""))


def hmmsearch(
    profile_hmm: Path,
    sequences_faa: Path,
    *,
    e_value: float = 1e-5,
    cpus: int = 0,
) -> pl.DataFrame:
    """Run pyhmmer hmmsearch and return a hit summary table."""

    if not profile_hmm.exists():
        raise MissingInputError(f"HMM profile does not exist: {profile_hmm}")
    if not sequences_faa.exists():
        raise MissingInputError(f"sequence FASTA does not exist: {sequences_faa}")

    with plan7.HMMFile(str(profile_hmm)) as hmm_file:
        queries = list(hmm_file)
    with easel.SequenceFile(str(sequences_faa), digital=True) as sequence_file:
        sequences = list(sequence_file)

    rows: list[dict[str, object]] = []
    search_results = hmmer.hmmsearch(
        queries,
        cast(Any, sequences),
        E=e_value,
        cpus=cpus,
    )
    for top_hits in search_results:
        query_name = _query_name(top_hits)
        for hit in top_hits:
            rows.append(
                {
                    "query_id": query_name,
                    "target_id": _decode_name(hit.name),
                    "evalue": float(hit.evalue),
                    "score": float(hit.score),
                    "bias": float(getattr(hit, "bias", 0.0)),
                    "included": bool(getattr(hit, "included", False)),
                },
            )

    return pl.DataFrame(
        rows,
        schema_overrides={
            "query_id": pl.Utf8,
            "target_id": pl.Utf8,
            "evalue": pl.Float64,
            "score": pl.Float64,
            "bias": pl.Float64,
            "included": pl.Boolean,
        },
    )
