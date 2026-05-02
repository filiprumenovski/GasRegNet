from __future__ import annotations

import os
import time
from pathlib import Path

import pytest
import yaml  # type: ignore[import-untyped]

from gasregnet.annotation.ncbi_taxonomy import build_taxonomy_db, lookup_lineage
from gasregnet.datasets.refseq import (
    detect_refseq_anchor_hits,
    extract_refseq_neighborhoods,
    index_refseq_corpus,
)
from scripts.generate_scaling_corpus import generate_scaling_corpus
from scripts.generate_scaling_taxdump import generate_scaling_taxdump

pytestmark = [
    pytest.mark.scaling,
    pytest.mark.skipif(
        os.environ.get("GASREGNET_ACCEPTANCE") != "1",
        reason="scaling acceptance probe is opt-in; set GASREGNET_ACCEPTANCE=1",
    ),
]


def _budget() -> dict[str, int | float]:
    path = Path(__file__).with_name("budgets.yaml")
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    return payload["small"]


def test_small_scaling_probe_stays_within_budgets(tmp_path: Path) -> None:
    budget = _budget()
    started = time.perf_counter()

    paths = generate_scaling_corpus(
        tmp_path / "corpus",
        datasets=int(budget["datasets"]),
        genes_per_dataset=int(budget["genes_per_dataset"]),
        window_genes=int(budget["window_genes"]),
    )
    taxdump = generate_scaling_taxdump(
        tmp_path / "taxdump",
        datasets=int(budget["datasets"]),
    )
    taxonomy_db = build_taxonomy_db(taxdump, tmp_path / "taxonomy.duckdb")
    index_refseq_corpus(paths.manifest, root=paths.root, taxonomy_db=taxonomy_db)

    anchor_hits = detect_refseq_anchor_hits(
        paths.manifest,
        paths.scan_config,
        root=paths.root,
        mode="smoke",
    )
    loci, genes = extract_refseq_neighborhoods(
        anchor_hits,
        paths.manifest,
        root=paths.root,
        window_genes=int(budget["window_genes"]),
    )
    elapsed = time.perf_counter() - started

    assert lookup_lineage(taxonomy_db, 900_000)["species"] == "Scaling species 0"
    assert int(budget["min_anchor_hits"]) <= anchor_hits.height <= int(
        budget["max_anchor_hits"],
    )
    assert int(budget["min_loci"]) <= loci.height <= int(budget["max_loci"])
    assert int(budget["min_genes"]) <= genes.height <= int(budget["max_genes"])
    assert elapsed <= float(budget["max_seconds"])
    assert set(anchor_hits["analyte"].to_list()) == {"CO", "CN"}
    assert not loci.filter(loci["is_boundary_truncated"]).is_empty()
