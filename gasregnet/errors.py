"""Typed exceptions used across GasRegNet."""

from __future__ import annotations


class GasRegNetError(Exception):
    """Base class for GasRegNet errors."""


class ConfigError(GasRegNetError):
    """Raised when configuration loading or validation fails."""


class SchemaError(GasRegNetError):
    """Raised when a table fails schema validation."""


class MissingInputError(GasRegNetError):
    """Raised when a required input file or binary is missing."""


class BenchmarkMissError(GasRegNetError):
    """Raised when a benchmark protein is absent from a corpus."""


class ScoringError(GasRegNetError):
    """Raised when score computation fails."""


class EnrichmentError(GasRegNetError):
    """Raised when enrichment computation fails."""
