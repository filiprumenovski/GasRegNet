"""Shard enumeration for partitioned corpus stores."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

import polars as pl

SHARDS_SCHEMA: dict[str, Any] = {
    "shard_id": pl.Utf8,
    "strategy": pl.Utf8,
    "phylum": pl.Utf8,
    "dataset_names": pl.List(pl.Utf8),
    "n_datasets": pl.Int64,
    "n_proteins_estimate": pl.Int64,
}


class CorpusStoreLike(Protocol):
    """Minimal interface shared with CorpusStoreHandle-like objects."""

    store_root: Path


@dataclass(frozen=True)
class Shard:
    """Resolved corpus shard."""

    shard_id: str
    strategy: str
    phylum: str
    dataset_names: tuple[str, ...]
    n_proteins_estimate: int

    @property
    def n_datasets(self) -> int:
        """Number of datasets assigned to the shard."""

        return len(self.dataset_names)


@dataclass(frozen=True)
class ShardStrategy:
    """Shard enumeration options."""

    name: str = "by_phylum"
    n_genomes_per_shard: int = 1000
    max_shards: int | None = None


DEFAULT_STRATEGY = ShardStrategy()


def _store_root(store: Path | CorpusStoreLike) -> Path:
    if isinstance(store, Path):
        return store
    return Path(store.store_root)


def _slug(value: str) -> str:
    normalized = re.sub(r"[^A-Za-z0-9]+", "-", value.strip().lower()).strip("-")
    return normalized or "unknown"


def _datasets_path(store_root: Path) -> Path:
    return store_root / "datasets.parquet"


def _read_datasets(store_root: Path) -> pl.DataFrame:
    path = _datasets_path(store_root)
    if not path.exists():
        raise FileNotFoundError(f"corpus store datasets table is missing: {path}")
    frame = pl.read_parquet(path)
    required = {"dataset_name", "phylum"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"datasets.parquet is missing required columns: {missing}")
    if "n_proteins" not in frame.columns:
        frame = frame.with_columns(pl.lit(0).alias("n_proteins"))
    return (
        frame.select(
            pl.col("dataset_name").cast(pl.Utf8),
            pl.when(pl.col("phylum").is_null() | (pl.col("phylum").cast(pl.Utf8) == ""))
            .then(pl.lit("unknown"))
            .otherwise(pl.col("phylum").cast(pl.Utf8))
            .alias("phylum"),
            pl.col("n_proteins").cast(pl.Int64, strict=False).fill_null(0),
        )
        .unique(subset=["dataset_name"], keep="first")
        .sort(["phylum", "dataset_name"])
    )


def _limit_shards(shards: list[Shard], max_shards: int | None) -> list[Shard]:
    if max_shards is None:
        return shards
    if max_shards < 1:
        raise ValueError("max_shards must be positive when provided")
    return shards[:max_shards]


def _by_phylum(frame: pl.DataFrame) -> list[Shard]:
    shards: list[Shard] = []
    for index, phylum in enumerate(sorted(frame["phylum"].unique().to_list()), start=1):
        group = frame.filter(pl.col("phylum") == phylum).sort("dataset_name")
        datasets = tuple(str(name) for name in group["dataset_name"].to_list())
        proteins = int(group["n_proteins"].sum())
        shards.append(
            Shard(
                shard_id=f"phylum-{index:04d}-{_slug(str(phylum))}",
                strategy="by_phylum",
                phylum=str(phylum),
                dataset_names=datasets,
                n_proteins_estimate=proteins,
            ),
        )
    return shards


def _by_n_genomes(frame: pl.DataFrame, n_genomes_per_shard: int) -> list[Shard]:
    if n_genomes_per_shard < 1:
        raise ValueError("n_genomes_per_shard must be positive")
    rows = list(frame.sort(["phylum", "dataset_name"]).iter_rows(named=True))
    shards: list[Shard] = []
    for start in range(0, len(rows), n_genomes_per_shard):
        chunk = rows[start : start + n_genomes_per_shard]
        phyla = sorted({str(row["phylum"]) for row in chunk})
        datasets = tuple(str(row["dataset_name"]) for row in chunk)
        proteins = sum(int(row["n_proteins"] or 0) for row in chunk)
        shards.append(
            Shard(
                shard_id=f"genomes-{len(shards) + 1:04d}",
                strategy="by_n_genomes",
                phylum=phyla[0] if len(phyla) == 1 else "mixed",
                dataset_names=datasets,
                n_proteins_estimate=proteins,
            ),
        )
    return shards


def enumerate_shards(
    store: Path | CorpusStoreLike,
    strategy: ShardStrategy | str = DEFAULT_STRATEGY,
    *,
    n_genomes_per_shard: int | None = None,
    max_shards: int | None = None,
) -> list[Shard]:
    """Enumerate deterministic corpus shards from ``datasets.parquet``."""

    resolved = (
        strategy
        if isinstance(strategy, ShardStrategy)
        else ShardStrategy(
            name=strategy,
            n_genomes_per_shard=n_genomes_per_shard or 1000,
            max_shards=max_shards,
        )
    )
    if n_genomes_per_shard is not None:
        resolved = ShardStrategy(
            name=resolved.name,
            n_genomes_per_shard=n_genomes_per_shard,
            max_shards=resolved.max_shards,
        )
    if max_shards is not None:
        resolved = ShardStrategy(
            name=resolved.name,
            n_genomes_per_shard=resolved.n_genomes_per_shard,
            max_shards=max_shards,
        )

    frame = _read_datasets(_store_root(store))
    if resolved.name == "by_phylum":
        shards = _by_phylum(frame)
    elif resolved.name == "by_n_genomes":
        shards = _by_n_genomes(frame, resolved.n_genomes_per_shard)
    else:
        raise ValueError(f"unsupported shard strategy: {resolved.name}")
    return _limit_shards(shards, resolved.max_shards)


def shards_to_frame(shards: list[Shard]) -> pl.DataFrame:
    """Convert shard objects to the persisted Parquet schema."""

    rows = [
        {
            "shard_id": shard.shard_id,
            "strategy": shard.strategy,
            "phylum": shard.phylum,
            "dataset_names": list(shard.dataset_names),
            "n_datasets": shard.n_datasets,
            "n_proteins_estimate": shard.n_proteins_estimate,
        }
        for shard in shards
    ]
    return pl.DataFrame(rows, schema=SHARDS_SCHEMA)


def write_shards(
    store: Path | CorpusStoreLike,
    out: Path | None = None,
    strategy: ShardStrategy | str = DEFAULT_STRATEGY,
    *,
    n_genomes_per_shard: int | None = None,
    max_shards: int | None = None,
) -> Path:
    """Enumerate and persist shards to ``shards.parquet``."""

    store_root = _store_root(store)
    output = out or (store_root / "shards.parquet")
    shards = enumerate_shards(
        store_root,
        strategy,
        n_genomes_per_shard=n_genomes_per_shard,
        max_shards=max_shards,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    shards_to_frame(shards).write_parquet(output)
    return output


def read_shards(path: Path) -> list[Shard]:
    """Read a persisted shard manifest."""

    frame = pl.read_parquet(path)
    missing = sorted(set(SHARDS_SCHEMA).difference(frame.columns))
    if missing:
        raise ValueError(f"shard manifest is missing required columns: {missing}")
    return [
        Shard(
            shard_id=str(row["shard_id"]),
            strategy=str(row["strategy"]),
            phylum=str(row["phylum"]),
            dataset_names=tuple(str(name) for name in row["dataset_names"]),
            n_proteins_estimate=int(row["n_proteins_estimate"] or 0),
        )
        for row in frame.sort("shard_id").iter_rows(named=True)
    ]


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--store", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument(
        "--strategy",
        choices=["by_phylum", "by_n_genomes"],
        default="by_phylum",
    )
    parser.add_argument("--n-genomes-per-shard", type=int, default=1000)
    parser.add_argument("--max-shards", type=int)
    args = parser.parse_args()
    print(
        write_shards(
            args.store,
            args.out,
            args.strategy,
            n_genomes_per_shard=args.n_genomes_per_shard,
            max_shards=args.max_shards,
        ),
    )


if __name__ == "__main__":
    _main()
