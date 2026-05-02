from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.benchmark import (
    evaluate_benchmark,
    load_benchmark_csv,
    summarize_benchmark_recovery,
)


def test_evaluate_benchmark_reports_positive_and_negative_recovery(
    tmp_path: Path,
) -> None:
    benchmark = tmp_path / "benchmark.csv"
    benchmark.write_text(
        "benchmark_id,analyte,protein_name,uniprot_accession,organism,taxon_id,"
        "anchor_family,expected_regulator_class,expected_sensory_domains,"
        "sensing_evidence_class,pmid,notes\n"
        "cn1,CN,cydA,,Escherichia coli,511145,cydA,unknown,[],inferred,[],hit\n"
        "neg1,negative_control,lacZ,,Escherichia coli,511145,,none,[],none,[],miss\n",
        encoding="utf-8",
    )
    anchor_hits = pl.DataFrame(
        {
            "dataset_name": ["ecoli_k12_mg1655"],
            "analyte": ["CN"],
            "anchor_family": ["cydA"],
            "protein_accession": ["NP_415261.2"],
            "locus_tag": ["b0733"],
            "gene": ["cydA"],
            "product": ["cytochrome bd-I ubiquinol oxidase subunit I"],
            "bitscore": [None],
            "e_value": [None],
            "identity": [None],
            "coverage": [None],
            "evidence_type": ["term_scan"],
        },
        schema_overrides={
            "bitscore": pl.Float64,
            "e_value": pl.Float64,
            "identity": pl.Float64,
            "coverage": pl.Float64,
        },
    )

    recovery = evaluate_benchmark(benchmark, anchor_hits)

    assert recovery["benchmark_id"].to_list() == ["cn1", "neg1"]
    assert recovery["hit"].to_list() == [True, True]
    assert recovery["is_negative_control"].to_list() == [False, True]


def test_summarize_benchmark_recovery_reports_validation_metrics(
    tmp_path: Path,
) -> None:
    benchmark = tmp_path / "benchmark.csv"
    benchmark.write_text(
        "benchmark_id,analyte,protein_name,uniprot_accession,organism,taxon_id,"
        "anchor_family,expected_regulator_class,expected_sensory_domains,"
        "sensing_evidence_class,pmid,notes,first_publication,verify_pmid\n"
        "pos1,CO,regA,,Org,1,coxL,one_component,[],direct,123,hit,paper,true\n"
        "neg1,negative_control,regB,,Org,1,,none,[],none,[],miss,,false\n",
        encoding="utf-8",
    )
    anchor_hits = pl.DataFrame(
        {
            "dataset_name": ["org"],
            "analyte": ["CO"],
            "anchor_family": ["coxL"],
            "protein_accession": ["regA"],
            "locus_tag": [""],
            "gene": ["regA"],
            "product": ["regA"],
            "bitscore": [1.0],
            "e_value": [1e-20],
            "identity": [1.0],
            "coverage": [1.0],
            "evidence_type": ["term_scan"],
        },
    )
    candidates = pl.DataFrame(
        {
            "gene_accession": ["regA", "regB"],
            "candidate_score": [10.0, 1.0],
            "regulation_logit_score": [0.9, 0.1],
        },
    )

    recovery = evaluate_benchmark(benchmark, anchor_hits, candidates)
    summary = summarize_benchmark_recovery(recovery)
    metrics = {row["metric"]: row for row in summary.iter_rows(named=True)}

    assert metrics["positive_recall"]["value"] == 1.0
    assert metrics["direct_positive_recall"]["value"] == 1.0
    assert metrics["negative_false_positive_rate"]["value"] == 0.0
    assert metrics["candidate_score_auroc"]["value"] == 1.0


def test_regulators_v2_benchmark_has_required_groups() -> None:
    benchmark = load_benchmark_csv(Path("configs/benchmarks/regulators_v2.csv"))
    counts = {
        str(row["analyte"]): int(row["len"])
        for row in benchmark.group_by("analyte").len().iter_rows(named=True)
    }

    assert counts["CO"] >= 10
    assert counts["NO"] >= 10
    assert counts["CN"] >= 10
    assert counts["O2"] >= 10
    assert counts["negative_control"] >= 10
    assert benchmark.filter(pl.col("verify_pmid")).height > 0
    accessions = benchmark.filter(pl.col("uniprot_accession").str.len_chars() > 0)
    assert accessions.height >= 28


def test_co_anchor_benchmark_curates_aerobic_and_anaerobic_codhs() -> None:
    benchmark = load_benchmark_csv(Path("data/benchmarks/anchors_v1.csv"))
    co_rows = benchmark.filter(pl.col("analyte") == "CO")

    cooS_rows = co_rows.filter(pl.col("anchor_family") == "cooS")
    coxL_rows = co_rows.filter(pl.col("anchor_family") == "coxL")

    assert cooS_rows.height >= 2
    assert coxL_rows.height >= 2
    assert cooS_rows.filter(pl.col("uniprot_accession") == "P31896").height == 1
    assert coxL_rows.filter(pl.col("uniprot_accession") == "P19919").height == 1
    assert cooS_rows.filter(pl.col("notes").str.contains("Anaerobic")).height >= 1
    assert coxL_rows.filter(pl.col("notes").str.contains("Aerobic")).height >= 1
