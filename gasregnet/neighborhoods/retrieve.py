"""Neighborhood retrieval helpers."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.errors import ConfigError
from gasregnet.io.sqlite_efi import read_efi_sqlite


def retrieve_from_efi_sqlite(
    sqlite_path: Path,
    analytes: list[str],
    *,
    cluster_filter: list[int] | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Retrieve loci and genes for one or more analytes from an EFI-GNT SQLite."""

    if not analytes:
        raise ConfigError("at least one analyte is required")

    loci_frames: list[pl.DataFrame] = []
    gene_frames: list[pl.DataFrame] = []
    for analyte in analytes:
        loci, genes = read_efi_sqlite(
            sqlite_path,
            analyte,
            cluster_filter=cluster_filter,
        )
        loci_frames.append(loci)
        gene_frames.append(genes)

    return (
        pl.concat(loci_frames, how="vertical").rechunk(),
        pl.concat(gene_frames, how="vertical").rechunk(),
    )
