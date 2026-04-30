from __future__ import annotations

from pathlib import Path

import pytest

from gasregnet.config import GasRegNetConfig, load_config, resolve_and_dump
from gasregnet.errors import ConfigError


def test_load_config_from_directory() -> None:
    config = load_config(Path("configs"))

    assert isinstance(config, GasRegNetConfig)
    assert [analyte.analyte for analyte in config.analytes] == ["CO", "CN"]
    assert config.seed == 20260429
    assert config.scoring.enrichment.case_control_ratio == (1, 3)


def test_load_config_from_headline_file() -> None:
    config = load_config(Path("configs/headline.yaml"))

    assert len(config.regulator_families) >= 10
    assert config.benchmark.positive_recall_threshold == 0.8


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
display_name: "hydrogen cyanide"
anchor_seeds: data/seeds/cn_anchor_seeds.faa
anchor_families: []
window_genes: 10
known_organisms_table: data/references/known_cn_organisms.csv
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
