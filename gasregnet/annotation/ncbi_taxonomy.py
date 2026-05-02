"""Offline NCBI taxdump lineage resolver."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import NamedTuple

import duckdb

LINEAGE_RANKS = (
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)


class TaxonNode(NamedTuple):
    """Minimal taxdump node fields needed for lineage resolution."""

    taxon_id: int
    parent_taxon_id: int
    rank: str


class LineageRow(NamedTuple):
    """Normalized ancestor row for a taxon."""

    taxon_id: int
    ancestor_taxon_id: int
    ancestor_rank: str
    ancestor_name: str
    distance: int


def _parse_dmp_fields(line: str) -> list[str]:
    raw_fields = line.rstrip("\n").removesuffix("\t|").split("\t|\t")
    return [field.strip() for field in raw_fields]


def _parse_nodes(path: Path) -> dict[int, TaxonNode]:
    nodes: dict[int, TaxonNode] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = _parse_dmp_fields(line)
            if len(fields) < 3:
                msg = f"{path}:{line_number} has fewer than 3 taxdump fields"
                raise ValueError(msg)
            taxon_id = int(fields[0])
            nodes[taxon_id] = TaxonNode(
                taxon_id=taxon_id,
                parent_taxon_id=int(fields[1]),
                rank=fields[2],
            )
    return nodes


def _parse_scientific_names(path: Path) -> dict[int, str]:
    names: dict[int, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = _parse_dmp_fields(line)
            if len(fields) < 4:
                msg = f"{path}:{line_number} has fewer than 4 taxdump fields"
                raise ValueError(msg)
            if fields[3] == "scientific name":
                names[int(fields[0])] = fields[1]
    return names


def _lineage_rows(
    nodes: dict[int, TaxonNode],
    names: dict[int, str],
) -> Iterator[LineageRow]:
    for taxon_id in sorted(nodes):
        current_id = taxon_id
        distance = 0
        seen: set[int] = set()
        while current_id in nodes:
            if current_id in seen:
                msg = f"cycle detected while resolving taxon {taxon_id}"
                raise ValueError(msg)
            seen.add(current_id)

            node = nodes[current_id]
            yield LineageRow(
                taxon_id=taxon_id,
                ancestor_taxon_id=node.taxon_id,
                ancestor_rank=node.rank,
                ancestor_name=names.get(node.taxon_id, ""),
                distance=distance,
            )
            if node.parent_taxon_id == node.taxon_id:
                break
            current_id = node.parent_taxon_id
            distance += 1


def _quote_identifier(identifier: str) -> str:
    escaped = identifier.replace('"', '""')
    return f'"{escaped}"'


def _create_lineage_pivot(connection: duckdb.DuckDBPyConnection) -> None:
    rank_columns = [
        (
            "max(ancestor_name) FILTER "
            f"(WHERE ancestor_rank = '{rank}') AS {_quote_identifier(rank)}"
        )
        for rank in LINEAGE_RANKS
    ]
    connection.execute(
        f"""
        CREATE TABLE lineage_pivot AS
        SELECT
            taxon_id,
            max(ancestor_name) FILTER (WHERE distance = 0) AS scientific_name,
            {", ".join(rank_columns)}
        FROM lineage
        GROUP BY taxon_id
        ORDER BY taxon_id
        """,
    )


def build_taxonomy_db(taxdump_dir: Path, out_db: Path) -> Path:
    """Build a DuckDB taxonomy lineage database from NCBI taxdump files.

    The input directory must contain uncompressed ``nodes.dmp`` and ``names.dmp``.
    The output database contains:

    - ``lineage``: one row per taxon/ancestor pair.
    - ``lineage_pivot``: one row per taxon with common lineage ranks as columns.
    """

    nodes_path = taxdump_dir / "nodes.dmp"
    names_path = taxdump_dir / "names.dmp"
    missing = [path.name for path in (nodes_path, names_path) if not path.exists()]
    if missing:
        missing_text = ", ".join(missing)
        msg = f"taxdump directory is missing required file(s): {missing_text}"
        raise FileNotFoundError(msg)

    nodes = _parse_nodes(nodes_path)
    names = _parse_scientific_names(names_path)

    out_db.parent.mkdir(parents=True, exist_ok=True)
    if out_db.exists():
        out_db.unlink()

    rows = list(_lineage_rows(nodes, names))
    with duckdb.connect(str(out_db)) as connection:
        connection.execute(
            """
            CREATE TABLE lineage (
                taxon_id INTEGER NOT NULL,
                ancestor_taxon_id INTEGER NOT NULL,
                ancestor_rank VARCHAR NOT NULL,
                ancestor_name VARCHAR NOT NULL,
                distance INTEGER NOT NULL
            )
            """,
        )
        connection.executemany(
            """
            INSERT INTO lineage
                (taxon_id, ancestor_taxon_id, ancestor_rank, ancestor_name, distance)
            VALUES (?, ?, ?, ?, ?)
            """,
            rows,
        )
        _create_lineage_pivot(connection)
        connection.execute("CREATE INDEX lineage_taxon_id_idx ON lineage(taxon_id)")
        connection.execute(
            "CREATE INDEX lineage_pivot_taxon_id_idx ON lineage_pivot(taxon_id)",
        )
    return out_db


def lookup_lineage(db: Path, taxon_id: int) -> dict[str, str]:
    """Return populated lineage rank names for ``taxon_id`` from a built database."""

    columns = ("scientific_name", *LINEAGE_RANKS)
    select_columns = ", ".join(_quote_identifier(column) for column in columns)
    with duckdb.connect(str(db), read_only=True) as connection:
        result = connection.execute(
            f"""
            SELECT {select_columns}
            FROM lineage_pivot
            WHERE taxon_id = ?
            """,
            [taxon_id],
        ).fetchone()
    if result is None:
        return {}
    return {
        column: value
        for column, value in zip(columns, result, strict=True)
        if isinstance(value, str) and value
    }
