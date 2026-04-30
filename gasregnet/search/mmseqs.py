"""MMseqs2 clustering helpers."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import polars as pl

from gasregnet.errors import MissingInputError

MMSEQS_CLUSTER_COLUMNS = ["cluster_id", "sequence_id"]


def _mmseqs_binary() -> str:
    binary = shutil.which("mmseqs")
    if binary is None:
        raise MissingInputError("mmseqs binary was not found on PATH")
    return binary


def parse_mmseqs_cluster_tsv(tsv: Path) -> pl.DataFrame:
    """Parse MMseqs cluster TSV into one row per sequence assignment."""

    if not tsv.exists():
        raise MissingInputError(f"MMseqs cluster output does not exist: {tsv}")

    clusters = pl.read_csv(
        tsv,
        separator="\t",
        has_header=False,
        new_columns=MMSEQS_CLUSTER_COLUMNS,
        schema_overrides={
            "cluster_id": pl.Utf8,
            "sequence_id": pl.Utf8,
        },
    )
    return clusters.with_columns(
        (pl.col("cluster_id") == pl.col("sequence_id")).alias("is_representative"),
    ).select("sequence_id", "cluster_id", "is_representative")


def cluster_sequences(
    input_faa: Path,
    out_dir: Path,
    *,
    min_seq_id: float = 0.5,
    coverage: float = 0.8,
    threads: int = 8,
) -> pl.DataFrame:
    """Run MMseqs easy-cluster and return parsed cluster assignments."""

    if not input_faa.exists():
        raise MissingInputError(f"input FASTA does not exist: {input_faa}")

    binary = _mmseqs_binary()
    out_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = out_dir / "clusters"
    tmp_dir = out_dir / "tmp"
    command = [
        binary,
        "easy-cluster",
        str(input_faa),
        str(output_prefix),
        str(tmp_dir),
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(coverage),
        "--cov-mode",
        "0",
        "--threads",
        str(threads),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return parse_mmseqs_cluster_tsv(out_dir / "clusters_cluster.tsv")
