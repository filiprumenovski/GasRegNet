"""HMMER search helpers."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from itertools import islice
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


def _batched(iterable: Iterable[Any], batch_size: int) -> Iterator[list[Any]]:
    iterator = iter(iterable)
    while batch := list(islice(iterator, batch_size)):
        yield batch


def _iter_target_blocks(targets: Any, block_residues: int) -> Iterator[Any]:
    if not hasattr(targets, "read_block"):
        yield targets
        return

    while True:
        block = targets.read_block(residues=block_residues)
        if len(block) == 0:
            break
        yield block


def hmmsearch(
    profile_hmm: Path,
    sequences_faa: Path,
    *,
    e_value: float = 1e-5,
    cpus: int = 0,
    block_residues: int = 100_000,
    query_batch_size: int = 32,
) -> pl.DataFrame:
    """Run pyhmmer hmmsearch and return a hit summary table."""

    if not profile_hmm.exists():
        raise MissingInputError(f"HMM profile does not exist: {profile_hmm}")
    if not sequences_faa.exists():
        raise MissingInputError(f"sequence FASTA does not exist: {sequences_faa}")
    if query_batch_size < 1:
        msg = "query_batch_size must be at least 1"
        raise ValueError(msg)

    rows: list[dict[str, object]] = []
    with plan7.HMMFile(str(profile_hmm)) as hmm_file:
        for queries in _batched(cast(Iterable[Any], hmm_file), query_batch_size):
            with easel.SequenceFile(str(sequences_faa), digital=True) as sequence_file:
                targets = cast(Any, sequence_file)
                for block in _iter_target_blocks(targets, block_residues):
                    search_results = hmmer.hmmsearch(
                        cast(Any, queries),
                        cast(Any, block),
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
