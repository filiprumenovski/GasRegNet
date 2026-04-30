"""Benchmark recovery metrics for corpus discovery outputs."""

from __future__ import annotations

from pathlib import Path

import polars as pl

BENCHMARK_RECOVERY_SCHEMA = {
    "benchmark_id": pl.Utf8,
    "analyte": pl.Utf8,
    "protein_name": pl.Utf8,
    "organism": pl.Utf8,
    "hit": pl.Boolean,
    "rank": pl.Int64,
    "candidate_score": pl.Float64,
}


def _contains(haystack: str, needle: str) -> bool:
    return bool(needle) and needle.lower() in haystack.lower()


def _candidate_ranks(candidates: pl.DataFrame | None) -> dict[str, tuple[int, float]]:
    if candidates is None or candidates.is_empty():
        return {}
    ranked = candidates.sort("candidate_score", descending=True).with_row_index(
        "rank",
        offset=1,
    )
    ranks: dict[str, tuple[int, float]] = {}
    for row in ranked.iter_rows(named=True):
        ranks[str(row["gene_accession"])] = (
            int(row["rank"]),
            float(row["candidate_score"]),
        )
    return ranks


def _anchor_hit(row: dict[str, object], anchor_hits: pl.DataFrame) -> bool:
    analyte = str(row["analyte"])
    protein_name = str(row["protein_name"])
    accession = str(row["uniprot_accession"])
    organism = str(row["organism"])
    subset = anchor_hits
    if analyte != "negative_control":
        subset = subset.filter(pl.col("analyte") == analyte)
    for hit in subset.iter_rows(named=True):
        searchable = " ".join(
            [
                str(hit["protein_accession"]),
                str(hit["locus_tag"]),
                str(hit["gene"]),
                str(hit["product"]),
                str(hit["anchor_family"]),
                str(hit["dataset_name"]),
            ],
        )
        if accession and _contains(searchable, accession):
            return True
        if _contains(searchable, protein_name):
            return True
        if organism and _contains(str(hit["dataset_name"]), organism):
            if _contains(searchable, protein_name):
                return True
    return False


def evaluate_benchmark(
    benchmark_csv: Path,
    anchor_hits: pl.DataFrame,
    candidates: pl.DataFrame | None = None,
) -> pl.DataFrame:
    """Evaluate benchmark recovery against anchor hits and ranked candidates."""

    benchmark = pl.read_csv(benchmark_csv)
    candidate_ranks = _candidate_ranks(candidates)
    rows: list[dict[str, object]] = []
    for benchmark_row in benchmark.iter_rows(named=True):
        is_negative = benchmark_row["analyte"] == "negative_control"
        recovered = _anchor_hit(benchmark_row, anchor_hits)
        hit = not recovered if is_negative else recovered
        rank: int | None = None
        score: float | None = None
        protein_name = str(benchmark_row["protein_name"])
        for gene_accession, rank_score in candidate_ranks.items():
            candidate_rank, candidate_score = rank_score
            if _contains(gene_accession, protein_name):
                rank = candidate_rank
                score = candidate_score
                break
        rows.append(
            {
                "benchmark_id": str(benchmark_row["benchmark_id"]),
                "analyte": str(benchmark_row["analyte"]),
                "protein_name": protein_name,
                "organism": str(benchmark_row["organism"]),
                "hit": hit,
                "rank": rank,
                "candidate_score": score,
            },
        )
    if not rows:
        return pl.DataFrame(schema=BENCHMARK_RECOVERY_SCHEMA)
    return pl.DataFrame(rows, schema_overrides=BENCHMARK_RECOVERY_SCHEMA)
