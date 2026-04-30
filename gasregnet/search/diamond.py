"""DIAMOND search wrapper."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import polars as pl

from gasregnet.errors import ExternalToolError, MissingInputError

DIAMOND_COLUMNS = [
    "query_id",
    "subject_id",
    "percent_identity",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovhsp",
    "scovhsp",
]
DIAMOND_OUTFMT = [
    "6",
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovhsp",
    "scovhsp",
]


def _diamond_binary() -> str:
    binary = shutil.which("diamond")
    if binary is None:
        raise MissingInputError("diamond binary was not found on PATH")
    return binary


def run_diamond(
    query_faa: Path,
    db_dmnd: Path,
    out_tsv: Path,
    *,
    evalue: float = 1e-10,
    coverage: float = 0.5,
    identity: float = 0.3,
    threads: int = 8,
    sensitivity: str = "very-sensitive",
) -> Path:
    """Run DIAMOND blastp and write a tabular output file."""

    if not query_faa.exists():
        raise MissingInputError(f"query FASTA does not exist: {query_faa}")
    if not db_dmnd.exists():
        raise MissingInputError(f"DIAMOND database does not exist: {db_dmnd}")

    binary = _diamond_binary()
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    command = [
        binary,
        "blastp",
        "--query",
        str(query_faa),
        "--db",
        str(db_dmnd),
        "--out",
        str(out_tsv),
        "--outfmt",
        *DIAMOND_OUTFMT,
        "--evalue",
        str(evalue),
        "--query-cover",
        str(coverage * 100.0),
        "--id",
        str(identity * 100.0),
        "--threads",
        str(threads),
        f"--{sensitivity}",
    ]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as error:
        stdout = (error.stdout or "").strip()
        stderr = (error.stderr or "").strip()
        context = [
            f"DIAMOND failed with exit code {error.returncode}",
            f"command: {' '.join(command)}",
        ]
        if stderr:
            context.append(f"stderr: {stderr}")
        if stdout:
            context.append(f"stdout: {stdout}")
        raise ExternalToolError("\n".join(context)) from error
    return out_tsv


def parse_diamond_output(tsv: Path) -> pl.DataFrame:
    """Parse DIAMOND outfmt 6 output into a typed Polars DataFrame."""

    if not tsv.exists():
        raise MissingInputError(f"DIAMOND output does not exist: {tsv}")

    return pl.read_csv(
        tsv,
        separator="\t",
        has_header=False,
        new_columns=DIAMOND_COLUMNS,
        schema_overrides={
            "query_id": pl.Utf8,
            "subject_id": pl.Utf8,
            "percent_identity": pl.Float64,
            "length": pl.Int64,
            "mismatch": pl.Int64,
            "gapopen": pl.Int64,
            "qstart": pl.Int64,
            "qend": pl.Int64,
            "sstart": pl.Int64,
            "send": pl.Int64,
            "evalue": pl.Float64,
            "bitscore": pl.Float64,
            "qcovhsp": pl.Float64,
            "scovhsp": pl.Float64,
        },
    )
