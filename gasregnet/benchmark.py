"""Benchmark recovery metrics for corpus discovery outputs."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.schemas import BenchmarkSchema, validate

BENCHMARK_RECOVERY_SCHEMA = {
    "benchmark_id": pl.Utf8,
    "analyte": pl.Utf8,
    "protein_name": pl.Utf8,
    "organism": pl.Utf8,
    "hit": pl.Boolean,
    "rank": pl.Int64,
    "candidate_score": pl.Float64,
    "regulation_posterior": pl.Float64,
}
BENCHMARK_LIST_COLUMNS = ("expected_sensory_domains", "pmid")


def _split_list_cell(value: object) -> list[str]:
    text = "" if value is None else str(value).strip()
    if text in {"", "[]"}:
        return []
    if text.startswith("[") and text.endswith("]"):
        text = text[1:-1].replace('"', "").replace("'", "")
    return [item.strip() for item in text.split("|") if item.strip()]


def load_benchmark_csv(path: Path) -> pl.DataFrame:
    """Load a benchmark CSV and coerce list-like columns for schema validation."""

    frame = pl.read_csv(path)
    if frame.is_empty():
        return validate(frame, BenchmarkSchema)
    string_columns = [
        "benchmark_id",
        "analyte",
        "protein_name",
        "uniprot_accession",
        "organism",
        "anchor_family",
        "expected_regulator_class",
        "sensing_evidence_class",
        "notes",
        "first_publication",
    ]
    frame = frame.with_columns(
        pl.col(column).fill_null("").cast(pl.Utf8)
        for column in string_columns
        if column in frame.columns
    )
    for column in BENCHMARK_LIST_COLUMNS:
        if column in frame.columns:
            values = [_split_list_cell(value) for value in frame[column].to_list()]
            frame = frame.with_columns(
                pl.Series(column, values, dtype=pl.List(pl.Utf8)),
            )
    if "verify_pmid" in frame.columns and frame["verify_pmid"].dtype != pl.Boolean:
        frame = frame.with_columns(
            pl.col("verify_pmid")
            .cast(pl.Utf8)
            .str.to_lowercase()
            .is_in(["true", "1", "yes"])
            .alias("verify_pmid"),
        )
    return validate(frame, BenchmarkSchema)


def write_benchmark_csv(frame: pl.DataFrame, path: Path) -> Path:
    """Validate and write a benchmark CSV with list columns serialized."""

    validated = validate(frame, BenchmarkSchema)
    output = validated.with_columns(
        pl.col(column).list.join("|").alias(column)
        for column in BENCHMARK_LIST_COLUMNS
        if column in validated.columns
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    output.write_csv(path)
    return path


def _contains(haystack: str, needle: str) -> bool:
    return bool(needle) and needle.lower() in haystack.lower()


def _candidate_ranks(
    candidates: pl.DataFrame | None,
) -> dict[str, tuple[int, float, float | None]]:
    if candidates is None or candidates.is_empty():
        return {}
    sort_column = (
        "regulation_posterior"
        if "regulation_posterior" in candidates.columns
        else "candidate_score"
    )
    ranked = candidates.sort(sort_column, descending=True).with_row_index(
        "rank",
        offset=1,
    )
    ranks: dict[str, tuple[int, float, float | None]] = {}
    for row in ranked.iter_rows(named=True):
        posterior = row.get("regulation_posterior")
        ranks[str(row["gene_accession"])] = (
            int(row["rank"]),
            float(row["candidate_score"]),
            None if posterior is None else float(posterior),
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

    benchmark = load_benchmark_csv(benchmark_csv)
    candidate_ranks = _candidate_ranks(candidates)
    rows: list[dict[str, object]] = []
    for benchmark_row in benchmark.iter_rows(named=True):
        is_negative = benchmark_row["analyte"] == "negative_control"
        recovered = _anchor_hit(benchmark_row, anchor_hits)
        hit = not recovered if is_negative else recovered
        rank: int | None = None
        score: float | None = None
        posterior: float | None = None
        protein_name = str(benchmark_row["protein_name"])
        for gene_accession, rank_score in candidate_ranks.items():
            candidate_rank, candidate_score, candidate_posterior = rank_score
            if _contains(gene_accession, protein_name):
                rank = candidate_rank
                score = candidate_score
                posterior = candidate_posterior
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
                "regulation_posterior": posterior,
            },
        )
    if not rows:
        return pl.DataFrame(schema=BENCHMARK_RECOVERY_SCHEMA)
    return pl.DataFrame(rows, schema_overrides=BENCHMARK_RECOVERY_SCHEMA)
