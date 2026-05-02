from __future__ import annotations

from pathlib import Path
from shutil import copytree

import pytest
from pydantic import ValidationError

from gasregnet.config import (
    BenchmarkConfig,
    GasRegNetConfig,
    load_config,
    resolve_and_dump,
)
from gasregnet.errors import ConfigError


def test_load_config_from_directory() -> None:
    config = load_config(Path("configs"))
    co_config = config.analytes[0]

    assert isinstance(config, GasRegNetConfig)
    assert [analyte.analyte for analyte in config.analytes] == [
        "CO",
        "NO",
        "CN",
        "cyd_control",
        "O2",
    ]
    assert config.seed == 20260429
    assert config.scoring.enrichment.case_control_ratio == (1, 3)
    assert config.benchmark.benchmark_csv is None
    assert config.benchmark.regulator_benchmark_csv == Path(
        "data/benchmarks/regulators_v2.csv",
    )
    assert config.sensory_domains[0].role == "sensor"
    assert config.sensory_domains[0].chemistry == "none"
    assert any(entry.role == "effector" for entry in config.sensory_domains)
    assert config.paired_evidence[0].rule_name == "cooA_heme_rescore"
    assert "sigma54_activator" in {
        entry.regulator_class for entry in config.regulator_families
    }
    assert co_config.expected_sensory_chemistry == [
        "heme",
        "iron_sulfur_4Fe4S",
    ]
    co_anchor_families = {
        family.name: (family.role, family.pfam_required, family.anchor_seeds)
        for family in co_config.anchor_families
    }
    assert co_anchor_families["coxL"] == (
        "primary",
        ["PF02738"],
        Path("data/seeds/co_coxL_anchor_seeds.faa"),
    )
    assert co_anchor_families["cooS"] == (
        "primary",
        ["PF03063"],
        Path("data/seeds/co_cooS_anchor_seeds.faa"),
    )
    assert config.analytes[2].expected_sensory_chemistry == [
        "heme",
        "iron_sulfur_4Fe4S",
        "redox_quinone",
        "cysteine_metal",
    ]


def test_load_config_from_headline_file() -> None:
    config = load_config(Path("configs/headline.yaml"))

    assert len(config.regulator_families) >= 10
    assert config.benchmark.positive_recall_threshold == 0.8


def test_load_config_rejects_unknown_regulator_family_class(tmp_path: Path) -> None:
    config_dir = tmp_path / "configs"
    copytree(Path("configs"), config_dir)
    (config_dir / "regulator_families.yaml").write_text(
        """
regulator_families:
  - family: invalid
    class: sigma_99
    pfam_required: ["PF00000"]
""",
        encoding="utf-8",
    )

    with pytest.raises(ConfigError, match="regulator_families.0.class"):
        load_config(config_dir)


def test_load_config_rejects_invalid_expected_sensory_chemistry(
    tmp_path: Path,
) -> None:
    config_dir = tmp_path / "configs"
    copytree(Path("configs"), config_dir)
    co_config = config_dir / "analytes" / "co.yaml"
    co_config.write_text(
        co_config.read_text(encoding="utf-8").replace(
            "iron_sulfur_4Fe4S",
            "iron_sulfur",
        ),
        encoding="utf-8",
    )

    with pytest.raises(ConfigError, match="expected_sensory_chemistry.1"):
        load_config(config_dir)


def test_benchmark_config_accepts_canonical_v2_paths_without_legacy_csv() -> None:
    config = BenchmarkConfig.model_validate(
        {
            "regulator_benchmark_csv": "data/benchmarks/regulators_v2.csv",
            "anchor_benchmark_csv": "data/benchmarks/anchors_v1.csv",
            "positive_recall_threshold": 0.8,
            "negative_false_positive_threshold": 0.1,
            "report_per_family": True,
        },
    )

    assert config.benchmark_csv is None
    assert config.regulator_benchmark_csv == Path("data/benchmarks/regulators_v2.csv")


def test_benchmark_config_requires_at_least_one_path() -> None:
    with pytest.raises(ValidationError, match="at least one benchmark path"):
        BenchmarkConfig.model_validate(
            {
                "positive_recall_threshold": 0.8,
                "negative_false_positive_threshold": 0.1,
                "report_per_family": True,
            },
        )


def test_load_config_rejects_invalid_paired_evidence() -> None:
    bad_config = Path("configs/sensory_domains.yaml").read_text(encoding="utf-8")
    assert "paired_evidence" in bad_config

    with pytest.raises(ConfigError, match="expected list field 'paired_evidence'"):
        from gasregnet.config import _load_optional_entries

        path = Path("/tmp/gasregnet_bad_sensory.yaml")
        path.write_text("sensory_domains: []\npaired_evidence: bad\n", encoding="utf-8")
        _load_optional_entries(path, "paired_evidence")


def test_resolve_and_dump_writes_yaml(tmp_path: Path) -> None:
    config = load_config(Path("configs"))

    output_path = resolve_and_dump(config, tmp_path)

    assert output_path == tmp_path / "config.resolved.yaml"
    assert "carbon monoxide" in output_path.read_text(encoding="utf-8")


def test_invalid_config_names_field_path(tmp_path: Path) -> None:
    config_dir = tmp_path / "configs"
    (config_dir / "analytes").mkdir(parents=True)
    (config_dir / "analytes" / "co.yaml").write_text(
        """
analyte: CO
display_name: "carbon monoxide"
anchor_seeds: data/seeds/co_anchor_seeds.faa
anchor_families: []
window_genes: 10
known_organisms_table: data/references/known_co_organisms.csv
expected_sensory_chemistry: []
seed: 20260429
""",
        encoding="utf-8",
    )
    (config_dir / "analytes" / "cn.yaml").write_text(
        """
analyte: CN
display_name: "cyanide-resistant respiration control"
anchor_seeds: data/seeds/cn_anchor_seeds.faa
anchor_families: []
window_genes: 10
known_organisms_table: data/references/known_cn_organisms.csv
expected_sensory_chemistry: []
seed: 20260429
""",
        encoding="utf-8",
    )
    (config_dir / "analytes" / "no.yaml").write_text(
        """
analyte: NO
display_name: "nitric oxide"
anchor_seeds: data/seeds/no_anchor_seeds.faa
anchor_families: []
window_genes: 10
known_organisms_table: data/references/known_no_organisms.csv
expected_sensory_chemistry: []
seed: 20260429
""",
        encoding="utf-8",
    )
    (config_dir / "analytes" / "o2.yaml").write_text(
        """
analyte: O2
display_name: "oxygen"
anchor_seeds: data/seeds/o2_anchor_seeds.faa
anchor_families: []
window_genes: 10
known_organisms_table: data/references/known_o2_organisms.csv
expected_sensory_chemistry: []
seed: 20260429
""",
        encoding="utf-8",
    )
    (config_dir / "regulator_families.yaml").write_text(
        "regulator_families: []\n",
        encoding="utf-8",
    )
    (config_dir / "sensory_domains.yaml").write_text(
        "sensory_domains: []\n",
        encoding="utf-8",
    )
    (config_dir / "scoring.yaml").write_text(
        """
locus_score_weights: {}
candidate_score_weights: {}
confidence_thresholds:
  high: 6.0
  medium: 3.5
enrichment:
  test: fisher
  multiple_comparison: benjamini_hochberg
  alpha: 0.05
  case_control_ratio: [1, 3]
  permutations: 10000
windows:
  strict: 5
  standard: 10
  extended: 20
  default: standard
robustness:
  windows_to_test: [5, 10, 20]
  weight_perturbation_pct: 20
""",
        encoding="utf-8",
    )
    (config_dir / "benchmarks.yaml").write_text(
        """
benchmark_csv: data/benchmarks/benchmark_v1.csv
positive_recall_threshold: 0.8
negative_false_positive_threshold: 0.1
report_per_family: true
""",
        encoding="utf-8",
    )

    with pytest.raises(ConfigError, match="confidence_thresholds.low"):
        load_config(config_dir)
