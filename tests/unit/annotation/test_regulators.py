from __future__ import annotations

import polars as pl

from gasregnet.annotation.regulators import classify_regulators
from gasregnet.config import load_config
from tests.unit.test_schemas import genes_frame


def _classify_with_pfams(pfam_ids: list[str], *, is_anchor: bool = False) -> str:
    config = load_config("configs")
    genes = genes_frame().with_columns(
        pl.Series("pfam_ids", [pfam_ids], dtype=pl.List(pl.Utf8)),
        pl.lit(is_anchor).alias("is_anchor"),
        pl.lit(0 if is_anchor else -1).cast(pl.Int32).alias("relative_index"),
        pl.lit("unknown").alias("functional_class"),
    )

    classified = classify_regulators(genes, config.regulator_families)
    return str(classified["regulator_class"].item())


def test_classify_regulators_matches_configured_one_component_family() -> None:
    assert _classify_with_pfams(["PF01047"]) == "one_component"


def test_classify_regulators_lets_specific_hth_beat_generic_receiver() -> None:
    assert _classify_with_pfams(["PF00072", "PF01047"]) == "one_component"


def test_classify_regulators_prioritizes_histidine_kinase() -> None:
    assert _classify_with_pfams(["PF00512"]) == "two_component_hk"


def test_classify_regulators_prioritizes_sigma() -> None:
    assert _classify_with_pfams(["PF04545"]) == "sigma"


def test_classify_regulators_preserves_anchor_functional_class() -> None:
    config = load_config("configs")
    genes = genes_frame().with_columns(
        pl.Series("pfam_ids", [["PF01047"]], dtype=pl.List(pl.Utf8)),
        pl.lit(True).alias("is_anchor"),
    )

    classified = classify_regulators(genes, config.regulator_families)

    assert classified["functional_class"].item() == "anchor"
    assert classified["is_regulator_candidate"].item() is True


def test_classify_regulators_marks_non_matches_as_unknown() -> None:
    config = load_config("configs")
    genes = genes_frame().with_columns(
        pl.Series("pfam_ids", [["PF99999"]], dtype=pl.List(pl.Utf8)),
        pl.lit(False).alias("is_anchor"),
        pl.lit(-1).cast(pl.Int32).alias("relative_index"),
        pl.lit("regulator").alias("functional_class"),
    )

    classified = classify_regulators(genes, config.regulator_families)

    assert classified["functional_class"].item() == "unknown"
    assert classified["regulator_class"].item() == "none"
    assert classified["is_regulator_candidate"].item() is False


def test_classify_regulators_uses_text_fallback_for_unannotated_refseq_rows() -> None:
    config = load_config("configs")
    genes = genes_frame().with_columns(
        pl.Series("pfam_ids", [[]], dtype=pl.List(pl.Utf8)),
        pl.lit("DNA-binding transcriptional repressor YbgA").alias(
            "product_description",
        ),
        pl.lit(False).alias("is_anchor"),
        pl.lit(-3).cast(pl.Int32).alias("relative_index"),
        pl.lit("unknown").alias("functional_class"),
    )

    classified = classify_regulators(genes, config.regulator_families)

    assert classified["functional_class"].item() == "regulator"
    assert classified["regulator_class"].item() == "one_component"
    assert classified["is_regulator_candidate"].item() is True
