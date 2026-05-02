from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.datasets.sharding import (
    ShardStrategy,
    enumerate_shards,
    read_shards,
    write_shards,
)


def _write_datasets(store: Path) -> None:
    store.mkdir(parents=True)
    pl.DataFrame(
        {
            "dataset_name": [
                "delta_2",
                "alpha_2",
                "alpha_1",
                "beta_1",
                "unknown_1",
            ],
            "phylum": [
                "Bacillota",
                "Pseudomonadota",
                "Pseudomonadota",
                "Actinomycetota",
                "",
            ],
            "n_proteins": ["20", "12", "10", "5", None],
        },
    ).write_parquet(store / "datasets.parquet")


def test_by_phylum_enumerates_deterministic_shards(tmp_path: Path) -> None:
    store = tmp_path / "store"
    _write_datasets(store)

    shards = enumerate_shards(store, "by_phylum")

    assert [shard.shard_id for shard in shards] == [
        "phylum-0001-actinomycetota",
        "phylum-0002-bacillota",
        "phylum-0003-pseudomonadota",
        "phylum-0004-unknown",
    ]
    assert shards[2].dataset_names == ("alpha_1", "alpha_2")
    assert shards[2].n_proteins_estimate == 22


def test_by_n_genomes_assigns_each_dataset_once(tmp_path: Path) -> None:
    store = tmp_path / "store"
    _write_datasets(store)

    shards = enumerate_shards(
        store,
        ShardStrategy(name="by_n_genomes", n_genomes_per_shard=2),
    )

    assert [shard.shard_id for shard in shards] == [
        "genomes-0001",
        "genomes-0002",
        "genomes-0003",
    ]
    assert all(shard.n_datasets <= 2 for shard in shards)
    assigned = [name for shard in shards for name in shard.dataset_names]
    assert assigned == [
        "beta_1",
        "delta_2",
        "alpha_1",
        "alpha_2",
        "unknown_1",
    ]
    assert sorted(assigned) == ["alpha_1", "alpha_2", "beta_1", "delta_2", "unknown_1"]


def test_write_shards_persists_manifest(tmp_path: Path) -> None:
    store = tmp_path / "store"
    _write_datasets(store)

    output = write_shards(store, strategy="by_n_genomes", n_genomes_per_shard=3)

    assert output == store / "shards.parquet"
    assert [shard.shard_id for shard in read_shards(output)] == [
        "genomes-0001",
        "genomes-0002",
    ]
