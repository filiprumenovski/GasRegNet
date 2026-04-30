from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.datasets.refseq import index_refseq_corpus
from scripts.run_corpus_discovery import run_corpus_discovery


def test_corpus_discovery_profile_mode_recovers_ecoli_cyd_anchors(
    tmp_path: Path,
) -> None:
    index_refseq_corpus(Path("configs/refseq_catalogs.yaml"), root=Path("."))

    out_dir = tmp_path / "corpus"
    marker = run_corpus_discovery(
        out_dir=out_dir,
        corpus_config_path=Path("configs/corpus_discovery.yaml"),
        config_path=Path("configs/headline.yaml"),
        root=Path("."),
    )

    anchor_hits = pl.read_parquet(out_dir / "intermediate" / "anchor_hits.parquet")
    recovered = set(anchor_hits["gene"].to_list())
    assert marker == out_dir / "README.txt"
    assert {"cydA", "cydB", "cydX"}.issubset(recovered)
    assert set(anchor_hits["evidence_type"].to_list()) == {"hmmer"}
    assert (out_dir / "intermediate" / "enrichment_robustness.parquet").exists()
