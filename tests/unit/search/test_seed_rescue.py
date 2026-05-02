from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl
import pytest

from gasregnet.config import AnalyteConfig, load_config
from gasregnet.errors import MissingInputError
from gasregnet.schemas import AnchorHitsSchema, validate
from gasregnet.search.seed_rescue import seed_rescue_for_shard


@dataclass(frozen=True)
class Shard:
    shard_id: str
    dataset_names: list[str]


class Store:
    def __init__(self, proteins: pl.DataFrame) -> None:
        self.proteins = proteins

    def fetch_proteins_for_shard(self, _shard: Shard) -> pl.DataFrame:
        return self.proteins


def _cyd_control_analyte() -> AnalyteConfig:
    config = load_config(Path("configs"))
    return next(
        analyte for analyte in config.analytes if analyte.analyte == "cyd_control"
    )


def _proteins() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "dataset_name": ["mini", "mini"],
            "protein_accession": ["NP_TRUE.1", "NP_DECOY.1"],
            "locus_tag": ["b0001", "b0002"],
            "gene": ["cydA", "foo"],
            "product": ["cytochrome bd oxidase subunit I", "hypothetical protein"],
            "sequence": ["MLDIV", "AAAAA"],
        },
    )


def test_seed_rescue_for_shard_runs_diamond_and_returns_anchor_hits(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seeds_dir = tmp_path / "diamond_seeds"
    seeds_dir.mkdir()
    (seeds_dir / "cyd_control__cydA.dmnd").write_text("", encoding="utf-8")
    calls: list[dict[str, object]] = []

    def fake_run_diamond(
        query_faa: Path,
        db_dmnd: Path,
        out_tsv: Path,
        *,
        evalue: float,
        coverage: float,
        identity: float,
        threads: int,
        sensitivity: str,
    ) -> Path:
        calls.append(
            {
                "query_faa": query_faa,
                "db_dmnd": db_dmnd,
                "evalue": evalue,
                "coverage": coverage,
                "identity": identity,
                "threads": threads,
                "sensitivity": sensitivity,
            },
        )
        assert query_faa.read_text(encoding="utf-8").startswith(">NP_TRUE.1\nMLDIV")
        out_tsv.write_text(
            "NP_TRUE.1\tseed_cydA\t62.0\t100\t0\t0\t1\t100\t1\t100\t1e-20\t80.0\t"
            "76.0\t90.0\n"
            "NP_DECOY.1\tseed_cydA\t20.0\t100\t0\t0\t1\t100\t1\t100\t1e-4\t10.0\t"
            "90.0\t90.0\n",
            encoding="utf-8",
        )
        return out_tsv

    monkeypatch.setattr(
        "gasregnet.search.seed_rescue.diamond.run_diamond",
        fake_run_diamond,
    )

    hits = seed_rescue_for_shard(
        shard=Shard("s1", ["mini"]),
        store=Store(_proteins()),
        seeds_diamond_dir=seeds_dir,
        analytes=[_cyd_control_analyte()],
        identity_threshold=0.30,
        coverage_threshold=0.50,
        e_value=1e-10,
        threads=4,
        work_dir=tmp_path / "work",
    )

    validate(hits, AnchorHitsSchema)
    assert calls == [
        {
            "query_faa": tmp_path / "work" / "shard_proteins.faa",
            "db_dmnd": seeds_dir / "cyd_control__cydA.dmnd",
            "evalue": 1e-10,
            "coverage": 0.50,
            "identity": 0.30,
            "threads": 4,
            "sensitivity": "ultra-sensitive",
        },
    ]
    assert hits.select(
        [
            "dataset_name",
            "analyte",
            "anchor_family",
            "protein_accession",
            "locus_tag",
            "gene",
            "evidence_type",
        ],
    ).row(0) == (
        "mini",
        "cyd_control",
        "cydA",
        "NP_TRUE.1",
        "b0001",
        "cydA",
        "seed_back_confirmed",
    )
    assert hits["identity"].item() == pytest.approx(0.62)
    assert hits["coverage"].item() == pytest.approx(0.76)


def test_seed_rescue_for_shard_keeps_best_hit_per_query(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seeds_dir = tmp_path / "diamond_seeds"
    seeds_dir.mkdir()
    (seeds_dir / "cyd_control__cydA.dmnd").write_text("", encoding="utf-8")

    def fake_run_diamond(
        _query_faa: Path,
        _db_dmnd: Path,
        out_tsv: Path,
        **_kwargs: object,
    ) -> Path:
        out_tsv.write_text(
            "NP_TRUE.1\tweak_seed\t50.0\t100\t0\t0\t1\t100\t1\t100\t1e-10\t40.0\t"
            "80.0\t90.0\n"
            "NP_TRUE.1\tbest_seed\t50.0\t100\t0\t0\t1\t100\t1\t100\t1e-30\t90.0\t"
            "80.0\t90.0\n",
            encoding="utf-8",
        )
        return out_tsv

    monkeypatch.setattr(
        "gasregnet.search.seed_rescue.diamond.run_diamond",
        fake_run_diamond,
    )

    hits = seed_rescue_for_shard(
        shard=Shard("s1", ["mini"]),
        store=Store(_proteins().head(1)),
        seeds_diamond_dir=seeds_dir,
        analytes=[_cyd_control_analyte()],
        work_dir=tmp_path / "work",
    )

    assert hits.height == 1
    assert hits["bitscore"].item() == pytest.approx(90.0)
    assert hits["e_value"].item() == pytest.approx(1e-30)


def test_seed_rescue_for_shard_skips_when_diamond_binary_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seeds_dir = tmp_path / "diamond_seeds"
    seeds_dir.mkdir()
    (seeds_dir / "cyd_control__cydA.dmnd").write_text("", encoding="utf-8")

    def fake_run_diamond(*_args: object, **_kwargs: object) -> Path:
        raise MissingInputError("diamond binary was not found on PATH")

    monkeypatch.setattr(
        "gasregnet.search.seed_rescue.diamond.run_diamond",
        fake_run_diamond,
    )

    hits = seed_rescue_for_shard(
        shard=Shard("s1", ["mini"]),
        store=Store(_proteins()),
        seeds_diamond_dir=seeds_dir,
        analytes=[_cyd_control_analyte()],
        work_dir=tmp_path / "work",
    )

    assert hits.is_empty()
    validate(hits, AnchorHitsSchema)
