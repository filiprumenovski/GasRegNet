"""Build DuckDB reference catalogs from RefSeq FASTA and GFF3 assets."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import duckdb
import polars as pl

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
