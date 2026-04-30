"""GFF3 parsing helpers."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import TextIO
from urllib.parse import unquote

import polars as pl

from gasregnet.errors import MissingInputError, SchemaError

GFF_COLUMNS = [
    "seqid",
    "source",
    "feature_type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_attributes(raw_attributes: str) -> dict[str, str]:
    """Parse a GFF3 attributes field."""

    attributes: dict[str, str] = {}
    if raw_attributes == ".":
        return attributes
    for item in raw_attributes.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise SchemaError(f"invalid GFF3 attribute item: {item!r}")
        key, value = item.split("=", 1)
        attributes[unquote(key)] = unquote(value)
    return attributes


def read_gff3(path: Path) -> pl.DataFrame:
    """Read a GFF3 file into a Polars DataFrame and reject non-GFF3 files."""

    if not path.exists():
        raise MissingInputError(f"GFF file does not exist: {path}")

    rows: list[dict[str, object]] = []
    saw_version = False
    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if line_number == 1 and line != "##gff-version 3":
                raise SchemaError(
                    "only GFF3 files with '##gff-version 3' are supported",
                )
            if line.startswith("##gff-version"):
                saw_version = True
                if line != "##gff-version 3":
                    raise SchemaError("only GFF3 files are supported")
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise SchemaError(f"GFF line {line_number} does not have 9 fields")
            attributes = parse_attributes(fields[8])
            rows.append(
                {
                    "seqid": fields[0],
                    "source": fields[1],
                    "feature_type": fields[2],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "score": None if fields[5] == "." else float(fields[5]),
                    "strand": fields[6],
                    "phase": None if fields[7] == "." else int(fields[7]),
                    "attributes": attributes,
                },
            )

    if not saw_version:
        raise SchemaError("missing GFF version directive")
    return pl.DataFrame(
        rows,
        schema_overrides={
            "start": pl.Int64,
            "end": pl.Int64,
            "score": pl.Float64,
            "phase": pl.Int64,
        },
    )
