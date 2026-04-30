"""Build DuckDB reference catalogs from RefSeq FASTA and GFF3 assets."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import duckdb
import polars as pl
import yaml  # type: ignore[import-untyped]

from gasregnet.io.fasta import read_fasta
from gasregnet.io.gff import read_gff3

PROTEIN_SCHEMA = {
    "protein_accession": pl.Utf8,
    "description": pl.Utf8,
    "sequence": pl.Utf8,
    "length_aa": pl.Int64,
}
FEATURE_SCHEMA = {
    "seqid": pl.Utf8,
    "source": pl.Utf8,
    "feature_type": pl.Utf8,
    "start_nt": pl.Int64,
    "end_nt": pl.Int64,
    "strand": pl.Utf8,
    "phase": pl.Int64,
    "feature_id": pl.Utf8,
    "parent_id": pl.Utf8,
    "locus_tag": pl.Utf8,
    "protein_accession": pl.Utf8,
    "gene": pl.Utf8,
    "product": pl.Utf8,
}
METADATA_SCHEMA = {
    "key": pl.Utf8,
    "value": pl.Utf8,
}
CATALOG_SCHEMA = {
    "dataset_name": pl.Utf8,
    "protein_faa": pl.Utf8,
    "gff": pl.Utf8,
    "out_db": pl.Utf8,
}


def _feature_value(attributes: dict[str, str], *keys: str) -> str | None:
    for key in keys:
        value = attributes.get(key)
        if value:
            return value
    return None


def _proteins_frame(protein_faa: Path) -> pl.DataFrame:
    rows = [
        {
            "protein_accession": accession,
            "description": description,
            "sequence": sequence,
            "length_aa": len(sequence),
        }
        for accession, description, sequence in read_fasta(protein_faa)
    ]
    return pl.DataFrame(rows, schema=PROTEIN_SCHEMA)


def _features_frame(gff: Path) -> pl.DataFrame:
    gff_frame = read_gff3(gff)
    if gff_frame.is_empty():
        return pl.DataFrame(schema=FEATURE_SCHEMA)
    rows: list[dict[str, object]] = []
    for row in gff_frame.iter_rows(named=True):
        attributes = row["attributes"]
        if not isinstance(attributes, dict):
            attributes = {}
        rows.append(
            {
                "seqid": row["seqid"],
                "source": row["source"],
                "feature_type": row["feature_type"],
                "start_nt": row["start"],
                "end_nt": row["end"],
                "strand": row["strand"],
                "phase": row["phase"],
                "feature_id": _feature_value(attributes, "ID"),
                "parent_id": _feature_value(attributes, "Parent"),
                "locus_tag": _feature_value(attributes, "locus_tag", "Name"),
                "protein_accession": _feature_value(
                    attributes,
                    "protein_id",
                    "protein_accession",
                ),
                "gene": _feature_value(attributes, "gene", "gene_synonym"),
                "product": _feature_value(attributes, "product", "Name"),
            },
        )
    return pl.DataFrame(rows, schema=FEATURE_SCHEMA)


def _metadata_frame(
    *,
    dataset_name: str,
    protein_faa: Path,
    gff: Path,
    proteins: pl.DataFrame,
    features: pl.DataFrame,
) -> pl.DataFrame:
    rows = [
        {"key": "dataset_name", "value": dataset_name},
        {"key": "protein_faa", "value": str(protein_faa)},
        {"key": "gff", "value": str(gff)},
        {"key": "n_proteins", "value": str(proteins.height)},
        {"key": "n_features", "value": str(features.height)},
        {
            "key": "created_at",
            "value": datetime.now(UTC).isoformat(timespec="seconds"),
        },
    ]
    return pl.DataFrame(rows, schema=METADATA_SCHEMA)


def index_refseq_dataset(
    *,
    protein_faa: Path,
    gff: Path,
    out_db: Path,
    dataset_name: str,
) -> Path:
    """Index RefSeq protein and annotation assets into a DuckDB database."""

    proteins = _proteins_frame(protein_faa)
    features = _features_frame(gff)
    metadata = _metadata_frame(
        dataset_name=dataset_name,
        protein_faa=protein_faa,
        gff=gff,
        proteins=proteins,
        features=features,
    )

    out_db.parent.mkdir(parents=True, exist_ok=True)
    if out_db.exists():
        out_db.unlink()
    with duckdb.connect(str(out_db)) as connection:
        connection.register("proteins_frame", proteins)
        connection.register("features_frame", features)
        connection.register("metadata_frame", metadata)
        connection.execute("create table proteins as select * from proteins_frame")
        connection.execute("create table features as select * from features_frame")
        connection.execute("create table metadata as select * from metadata_frame")
        connection.execute(
            "create index proteins_accession_idx on proteins(protein_accession)",
        )
        connection.execute(
            "create index features_protein_idx on features(protein_accession)",
        )
        connection.execute("create index features_locus_idx on features(locus_tag)")
    return out_db


def read_refseq_catalog_manifest(manifest: Path, *, root: Path) -> pl.DataFrame:
    """Read a RefSeq catalog manifest into a normalized table."""

    payload = yaml.safe_load(manifest.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("catalogs"), list):
        raise ValueError(f"catalog manifest must contain a catalogs list: {manifest}")
    rows: list[dict[str, str]] = []
    for raw in payload["catalogs"]:
        if not isinstance(raw, dict):
            raise ValueError("catalog entries must be mappings")
        dataset_name = raw.get("dataset_name")
        protein_faa = raw.get("protein_faa")
        gff = raw.get("gff")
        out_db = raw.get("out_db")
        if not isinstance(dataset_name, str) or not dataset_name:
            raise ValueError("catalog entry is missing string field: dataset_name")
        if not isinstance(protein_faa, str) or not protein_faa:
            raise ValueError(f"catalog {dataset_name} is missing protein_faa")
        if not isinstance(gff, str) or not gff:
            raise ValueError(f"catalog {dataset_name} is missing gff")
        if not isinstance(out_db, str) or not out_db:
            raise ValueError(f"catalog {dataset_name} is missing out_db")
        rows.append(
            {
                "dataset_name": dataset_name,
                "protein_faa": str(root / protein_faa),
                "gff": str(root / gff),
                "out_db": str(root / out_db),
            },
        )
    return pl.DataFrame(rows, schema=CATALOG_SCHEMA)


def index_refseq_corpus(manifest: Path, *, root: Path) -> list[Path]:
    """Index every RefSeq catalog declared in a corpus manifest."""

    catalogs = read_refseq_catalog_manifest(manifest, root=root)
    outputs: list[Path] = []
    for row in catalogs.iter_rows(named=True):
        outputs.append(
            index_refseq_dataset(
                protein_faa=Path(str(row["protein_faa"])),
                gff=Path(str(row["gff"])),
                out_db=Path(str(row["out_db"])),
                dataset_name=str(row["dataset_name"]),
            ),
        )
    return outputs


def _metadata_map(connection: duckdb.DuckDBPyConnection) -> dict[str, str]:
    rows = connection.execute("select key, value from metadata").fetchall()
    return {str(key): str(value) for key, value in rows}


def _scalar(connection: duckdb.DuckDBPyConnection, query: str) -> Any:
    row = connection.execute(query).fetchone()
    if row is None:
        raise ValueError(f"query returned no rows: {query}")
    return row[0]


def summarize_refseq_catalog(db: Path) -> dict[str, Any]:
    """Return summary statistics for one RefSeq DuckDB catalog."""

    with duckdb.connect(str(db), read_only=True) as connection:
        metadata = _metadata_map(connection)
        n_cds = _scalar(
            connection,
            "select count(*) from features where feature_type = 'CDS'",
        )
        n_linked_features = _scalar(
            connection,
            "select count(*) from features where protein_accession is not null",
        )
        mean_length = _scalar(
            connection,
            "select avg(length_aa) from proteins",
        )
    return {
        "dataset_name": metadata["dataset_name"],
        "db": str(db),
        "n_proteins": int(metadata["n_proteins"]),
        "n_features": int(metadata["n_features"]),
        "n_cds": int(n_cds),
        "n_features_with_protein": int(n_linked_features),
        "mean_protein_length_aa": float(mean_length or 0.0),
    }


def summarize_refseq_corpus(manifest: Path, *, root: Path) -> pl.DataFrame:
    """Summarize every RefSeq catalog declared in a corpus manifest."""

    catalogs = read_refseq_catalog_manifest(manifest, root=root)
    rows = [
        summarize_refseq_catalog(Path(str(row["out_db"])))
        for row in catalogs.iter_rows(named=True)
    ]
    return pl.DataFrame(
        rows,
        schema={
            "dataset_name": pl.Utf8,
            "db": pl.Utf8,
            "n_proteins": pl.Int64,
            "n_features": pl.Int64,
            "n_cds": pl.Int64,
            "n_features_with_protein": pl.Int64,
            "mean_protein_length_aa": pl.Float64,
        },
    )


def query_refseq_catalog(db: Path, query: str, *, limit: int = 20) -> pl.DataFrame:
    """Search a DuckDB reference catalog by accession, locus tag, gene, or product."""

    pattern = f"%{query}%"
    with duckdb.connect(str(db), read_only=True) as connection:
        return connection.execute(
            """
            select
                proteins.protein_accession,
                proteins.length_aa,
                coalesce(features.locus_tag, '') as locus_tag,
                coalesce(features.gene, '') as gene,
                coalesce(features.product, proteins.description) as product,
                coalesce(features.seqid, '') as seqid,
                features.start_nt,
                features.end_nt,
                coalesce(features.strand, '') as strand
            from proteins
            left join features using (protein_accession)
            where
                proteins.protein_accession ilike ?
                or proteins.description ilike ?
                or features.locus_tag ilike ?
                or features.gene ilike ?
                or features.product ilike ?
            order by proteins.protein_accession, features.start_nt nulls last
            limit ?
            """,
            [pattern, pattern, pattern, pattern, pattern, limit],
        ).pl()


def query_refseq_corpus(
    manifest: Path,
    query: str,
    *,
    root: Path,
    limit_per_catalog: int = 20,
) -> pl.DataFrame:
    """Search every RefSeq catalog declared in a corpus manifest."""

    catalogs = read_refseq_catalog_manifest(manifest, root=root)
    frames: list[pl.DataFrame] = []
    for row in catalogs.iter_rows(named=True):
        dataset_name = str(row["dataset_name"])
        frame = query_refseq_catalog(
            Path(str(row["out_db"])),
            query,
            limit=limit_per_catalog,
        )
        if not frame.is_empty():
            frames.append(
                frame.with_columns(pl.lit(dataset_name).alias("dataset_name")),
            )
    if not frames:
        return pl.DataFrame(
            schema={
                "protein_accession": pl.Utf8,
                "length_aa": pl.Int64,
                "locus_tag": pl.Utf8,
                "gene": pl.Utf8,
                "product": pl.Utf8,
                "seqid": pl.Utf8,
                "start_nt": pl.Int64,
                "end_nt": pl.Int64,
                "strand": pl.Utf8,
                "dataset_name": pl.Utf8,
            },
        )
    return pl.concat(frames, how="vertical")
