"""Configuration models and loaders."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

import yaml  # type: ignore[import-untyped]
from pydantic import BaseModel, ConfigDict, Field, ValidationError

from gasregnet.errors import ConfigError


class StrictModel(BaseModel):
    """Base model that rejects unknown config fields."""

    model_config = ConfigDict(extra="forbid")


class AnchorFamilyConfig(StrictModel):
    name: str
    pfam_required: list[str]
    role: Literal["primary", "accessory", "contextual"]
    pfam_supporting: list[str] = Field(default_factory=list)


class AnalyteConfig(StrictModel):
    analyte: Literal["CO", "CN"]
    display_name: str
    anchor_seeds: Path
    anchor_families: list[AnchorFamilyConfig]
    window_genes: int = Field(gt=0)
    known_organisms_table: Path
    expected_sensory_chemistry: list[str]
    seed: int


class EnrichmentConfig(StrictModel):
    test: Literal["fisher"]
    multiple_comparison: Literal["benjamini_hochberg"]
    alpha: float = Field(gt=0.0, le=1.0)
    case_control_ratio: tuple[int, int]
    permutations: int = Field(gt=0)


class WindowConfig(StrictModel):
    strict: int = Field(gt=0)
    standard: int = Field(gt=0)
    extended: int = Field(gt=0)
    default: Literal["strict", "standard", "extended"]


class RobustnessConfig(StrictModel):
    windows_to_test: list[int]
    weight_perturbation_pct: int = Field(ge=0)


class ConfidenceThresholdsConfig(StrictModel):
    high: float
    medium: float
    low: float


class ScoringConfig(StrictModel):
    locus_score_weights: dict[str, float]
    candidate_score_weights: dict[str, float]
    confidence_thresholds: ConfidenceThresholdsConfig
    enrichment: EnrichmentConfig
    windows: WindowConfig
    robustness: RobustnessConfig


class RegulatorFamilyEntry(StrictModel):
    family: str
    regulator_class: str = Field(alias="class")
    pfam_required: list[str]
    pfam_optional: list[str] = Field(default_factory=list)


class SensoryDomainEntry(StrictModel):
    domain: str
    pfam_id: str
    chemistry: Literal[
        "heme",
        "iron_sulfur",
        "cysteine_rich",
        "flavin",
        "redox",
        "metal",
        "cnmp",
        "generic_signal",
    ]


class BenchmarkConfig(StrictModel):
    benchmark_csv: Path
    positive_recall_threshold: float = Field(ge=0.0, le=1.0)
    negative_false_positive_threshold: float = Field(ge=0.0, le=1.0)
    report_per_family: bool


class GasRegNetConfig(StrictModel):
    analytes: list[AnalyteConfig]
    regulator_families: list[RegulatorFamilyEntry]
    sensory_domains: list[SensoryDomainEntry]
    scoring: ScoringConfig
    benchmark: BenchmarkConfig
    window: Literal["strict", "standard", "extended"] = "standard"
    seed: int = 20260429


def _read_yaml(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    except OSError as exc:
        raise ConfigError(f"could not read config file {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError(f"config file {path} must contain a YAML mapping")
    return data


def _load_entries(path: Path, key: str) -> list[dict[str, Any]]:
    data = _read_yaml(path)
    entries = data.get(key)
    if not isinstance(entries, list):
        raise ConfigError(f"{path}: expected list field {key!r}")
    return entries


def _compose_from_dir(config_dir: Path) -> dict[str, Any]:
    return {
        "analytes": [
            _read_yaml(config_dir / "analytes" / "co.yaml"),
            _read_yaml(config_dir / "analytes" / "cn.yaml"),
        ],
        "regulator_families": _load_entries(
            config_dir / "regulator_families.yaml",
            "regulator_families",
        ),
        "sensory_domains": _load_entries(
            config_dir / "sensory_domains.yaml",
            "sensory_domains",
        ),
        "scoring": _read_yaml(config_dir / "scoring.yaml"),
        "benchmark": _read_yaml(config_dir / "benchmarks.yaml"),
        "window": "standard",
        "seed": 20260429,
    }


def _resolve_reference(base_dir: Path, value: object) -> Path:
    if not isinstance(value, str):
        raise ConfigError(f"expected config path string, got {type(value).__name__}")
    path = Path(value)
    if not path.is_absolute() and not path.exists():
        path = base_dir / path
    return path


def _compose_from_headline(path: Path) -> dict[str, Any]:
    headline = _read_yaml(path)
    base_dir = path.parent
    analyte_refs = headline.get("analytes")
    if not isinstance(analyte_refs, list):
        raise ConfigError(f"{path}: expected list field 'analytes'")
    return {
        "analytes": [
            _read_yaml(_resolve_reference(base_dir, analyte_ref))
            for analyte_ref in analyte_refs
        ],
        "regulator_families": _load_entries(
            _resolve_reference(base_dir, headline.get("regulator_families")),
            "regulator_families",
        ),
        "sensory_domains": _load_entries(
            _resolve_reference(base_dir, headline.get("sensory_domains")),
            "sensory_domains",
        ),
        "scoring": _read_yaml(_resolve_reference(base_dir, headline.get("scoring"))),
        "benchmark": _read_yaml(
            _resolve_reference(base_dir, headline.get("benchmark")),
        ),
        "window": headline.get("window", "standard"),
        "seed": headline.get("seed", 20260429),
    }


def _format_validation_error(exc: ValidationError) -> str:
    errors = []
    for error in exc.errors():
        location = ".".join(str(part) for part in error["loc"])
        errors.append(f"{location}: {error['msg']}")
    return "; ".join(errors)


def load_config(config_dir: Path) -> GasRegNetConfig:
    """Load a composed GasRegNet config from a directory or headline YAML file."""

    path = Path(config_dir)
    try:
        raw = _compose_from_dir(path) if path.is_dir() else _compose_from_headline(path)
        return GasRegNetConfig.model_validate(raw)
    except ValidationError as exc:
        raise ConfigError(_format_validation_error(exc)) from exc


def resolve_and_dump(config: GasRegNetConfig, out_dir: Path) -> Path:
    """Write the resolved config YAML into a run output directory."""

    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / "config.resolved.yaml"
    data = config.model_dump(mode="json", by_alias=True)
    with output_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False)
    return output_path
