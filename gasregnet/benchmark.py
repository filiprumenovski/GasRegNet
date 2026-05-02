"""Benchmark recovery metrics for corpus discovery outputs."""

from __future__ import annotations

from pathlib import Path
from typing import cast

import polars as pl

from gasregnet.schemas import BenchmarkSchema, validate

BENCHMARK_RECOVERY_SCHEMA = {
    "benchmark_id": pl.Utf8,
    "analyte": pl.Utf8,
    "sensing_evidence_class": pl.Utf8,
    "is_negative_control": pl.Boolean,
    "verified_pmid": pl.Boolean,
    "protein_name": pl.Utf8,
    "organism": pl.Utf8,
    "hit": pl.Boolean,
    "rank": pl.Int64,
    "candidate_score": pl.Float64,
    "regulation_logit_score": pl.Float64,
}
BENCHMARK_SUMMARY_SCHEMA = {
    "metric": pl.Utf8,
    "value": pl.Float64,
    "n": pl.Int64,
    "notes": pl.Utf8,
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
        "regulation_logit_score"
        if "regulation_logit_score" in candidates.columns
        else "candidate_score"
    )
    ranked = candidates.sort(sort_column, descending=True).with_row_index(
        "rank",
        offset=1,
    )
    ranks: dict[str, tuple[int, float, float | None]] = {}
    for row in ranked.iter_rows(named=True):
        logit_score = row.get("regulation_logit_score")
        ranks[str(row["gene_accession"])] = (
            int(row["rank"]),
            float(row["candidate_score"]),
            None if logit_score is None else float(logit_score),
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
        logit_score: float | None = None
        protein_name = str(benchmark_row["protein_name"])
        for gene_accession, rank_score in candidate_ranks.items():
            candidate_rank, candidate_score, candidate_logit_score = rank_score
            if _contains(gene_accession, protein_name):
                rank = candidate_rank
                score = candidate_score
                logit_score = candidate_logit_score
                break
        rows.append(
            {
                "benchmark_id": str(benchmark_row["benchmark_id"]),
                "analyte": str(benchmark_row["analyte"]),
                "sensing_evidence_class": str(benchmark_row["sensing_evidence_class"]),
                "is_negative_control": is_negative,
                "verified_pmid": bool(benchmark_row.get("verify_pmid", False)),
                "protein_name": protein_name,
                "organism": str(benchmark_row["organism"]),
                "hit": hit,
                "rank": rank,
                "candidate_score": score,
                "regulation_logit_score": logit_score,
            },
        )
    if not rows:
        return pl.DataFrame(schema=BENCHMARK_RECOVERY_SCHEMA)
    return pl.DataFrame(rows, schema_overrides=BENCHMARK_RECOVERY_SCHEMA)


def _binary_auc(labels: list[int], scores: list[float]) -> float | None:
    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return None
    ordered = sorted(zip(scores, labels, strict=True), key=lambda item: item[0])
    rank_sum = 0.0
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][0] == ordered[index][0]:
            end += 1
        average_rank = (index + 1 + end) / 2.0
        rank_sum += average_rank * sum(label for _, label in ordered[index:end])
        index = end
    return (rank_sum - positives * (positives + 1) / 2.0) / (positives * negatives)


def _average_precision(labels: list[int], scores: list[float]) -> float | None:
    positives = sum(labels)
    if positives == 0:
        return None
    ordered = sorted(
        zip(scores, labels, strict=True),
        key=lambda item: item[0],
        reverse=True,
    )
    recovered = 0
    precision_sum = 0.0
    for rank, (_, label) in enumerate(ordered, start=1):
        if label == 1:
            recovered += 1
            precision_sum += recovered / rank
    return precision_sum / positives


def _metric_row(
    metric: str,
    value: float | None,
    n: int,
    notes: str,
) -> dict[str, object]:
    return {
        "metric": metric,
        "value": None if value is None else float(value),
        "n": n,
        "notes": notes,
    }


def _mean_bool(series: pl.Series) -> float | None:
    value = series.mean()
    return None if value is None else float(cast(float, value))


def summarize_benchmark_recovery(recovery: pl.DataFrame) -> pl.DataFrame:
    """Summarize recovery with recall, FPR, top-k, and ranking metrics."""

    if recovery.is_empty():
        return pl.DataFrame(schema=BENCHMARK_SUMMARY_SCHEMA)
    if "is_negative_control" not in recovery.columns:
        recovery = recovery.with_columns(
            (pl.col("analyte") == "negative_control").alias("is_negative_control"),
        )
    if "sensing_evidence_class" not in recovery.columns:
        recovery = recovery.with_columns(
            pl.lit("unspecified").alias("sensing_evidence_class"),
        )
    if "verified_pmid" not in recovery.columns:
        recovery = recovery.with_columns(pl.lit(False).alias("verified_pmid"))

    positives = recovery.filter(~pl.col("is_negative_control"))
    negatives = recovery.filter(pl.col("is_negative_control"))
    direct = positives.filter(pl.col("sensing_evidence_class") == "direct")
    verified = positives.filter(pl.col("verified_pmid"))

    rows: list[dict[str, object]] = [
        _metric_row(
            "positive_recall",
            _mean_bool(positives["hit"]) if positives.height else None,
            positives.height,
            "All non-negative benchmark rows.",
        ),
        _metric_row(
            "direct_positive_recall",
            _mean_bool(direct["hit"]) if direct.height else None,
            direct.height,
            "Direct evidence positives only.",
        ),
        _metric_row(
            "verified_positive_recall",
            _mean_bool(verified["hit"]) if verified.height else None,
            verified.height,
            "Positive rows with verify_pmid=true.",
        ),
        _metric_row(
            "negative_false_positive_rate",
            _mean_bool(~negatives["hit"]) if negatives.height else None,
            negatives.height,
            "Negative-control rows where a recovery was observed.",
        ),
    ]
    for k in (1, 5, 10):
        ranked = positives.filter(pl.col("rank").is_not_null())
        rows.append(
            _metric_row(
                f"top_{k}_recall",
                (
                    _mean_bool(ranked["rank"] <= k)
                    if ranked.height and positives.height
                    else None
                ),
                ranked.height,
                f"Positive benchmark rows with candidate rank <= {k}.",
            ),
        )

    scored = recovery.filter(pl.col("candidate_score").is_not_null())
    if scored.height:
        labels = [
            0 if bool(row["is_negative_control"]) else 1
            for row in scored.iter_rows(named=True)
        ]
        scores = [float(row["candidate_score"]) for row in scored.iter_rows(named=True)]
        rows.extend(
            [
                _metric_row(
                    "candidate_score_auroc",
                    _binary_auc(labels, scores),
                    scored.height,
                    "AUROC over benchmark rows with candidate scores.",
                ),
                _metric_row(
                    "candidate_score_auprc",
                    _average_precision(labels, scores),
                    scored.height,
                    "Average precision over benchmark rows with candidate scores.",
                ),
            ],
        )
    return pl.DataFrame(rows, schema_overrides=BENCHMARK_SUMMARY_SCHEMA)
