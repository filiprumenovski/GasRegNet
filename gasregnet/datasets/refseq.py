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
from gasregnet.schemas import AnchorHitsSchema, GenesSchema, LociSchema, validate

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
SCAN_TARGET_SCHEMA = {
    "analyte": pl.Utf8,
    "term": pl.Utf8,
}
SCAN_RESULT_SCHEMA = {
    "dataset_name": pl.Utf8,
    "analyte": pl.Utf8,
    "term": pl.Utf8,
    "protein_accession": pl.Utf8,
    "length_aa": pl.Int64,
    "locus_tag": pl.Utf8,
    "gene": pl.Utf8,
    "product": pl.Utf8,
    "seqid": pl.Utf8,
    "start_nt": pl.Int64,
    "end_nt": pl.Int64,
    "strand": pl.Utf8,
}
ANCHOR_HITS_SCHEMA = {
    "dataset_name": pl.Utf8,
    "analyte": pl.Utf8,
    "anchor_family": pl.Utf8,
    "protein_accession": pl.Utf8,
    "locus_tag": pl.Utf8,
    "gene": pl.Utf8,
    "product": pl.Utf8,
    "bitscore": pl.Float64,
    "e_value": pl.Float64,
    "identity": pl.Float64,
    "coverage": pl.Float64,
    "evidence_type": pl.Utf8,
}
LOCI_OUTPUT_SCHEMA: dict[str, Any] = {
    "locus_id": pl.Utf8,
    "analyte": pl.Utf8,
    "anchor_accession": pl.Utf8,
    "anchor_family": pl.Utf8,
    "organism": pl.Utf8,
    "taxon_id": pl.Int64,
    "cluster_id": pl.Int32,
    "contig_id": pl.Utf8,
    "window_size": pl.Int32,
    "is_boundary_truncated": pl.Boolean,
    "marker_genes_present": pl.List(pl.Utf8),
    "accessory_genes_present": pl.List(pl.Utf8),
    "locus_score": pl.Float64,
    "locus_confidence": pl.Utf8,
    "taxonomic_context_score": pl.Float64,
    "operon_integrity_score": pl.Float64,
    "created_at": pl.Datetime("us"),
}
GENES_OUTPUT_SCHEMA: dict[str, Any] = {
    "locus_id": pl.Utf8,
    "gene_accession": pl.Utf8,
    "relative_index": pl.Int32,
    "relative_start": pl.Int64,
    "relative_stop": pl.Int64,
    "strand": pl.Utf8,
    "product_description": pl.Utf8,
    "pfam_ids": pl.List(pl.Utf8),
    "pfam_descriptions": pl.List(pl.Utf8),
    "interpro_ids": pl.List(pl.Utf8),
    "interpro_descriptions": pl.List(pl.Utf8),
    "functional_class": pl.Utf8,
    "regulator_class": pl.Utf8,
    "sensory_domains": pl.List(pl.Utf8),
    "is_anchor": pl.Boolean,
    "is_regulator_candidate": pl.Boolean,
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


def read_refseq_scan_config(path: Path) -> pl.DataFrame:
    """Read analyte search terms for RefSeq catalog scans."""

    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("targets"), list):
        raise ValueError(f"scan config must contain a targets list: {path}")
    rows: list[dict[str, str]] = []
    for raw in payload["targets"]:
        if not isinstance(raw, dict):
            raise ValueError("scan target entries must be mappings")
        analyte = raw.get("analyte")
        terms = raw.get("terms")
        if not isinstance(analyte, str) or not analyte:
            raise ValueError("scan target is missing string field: analyte")
        if not isinstance(terms, list) or not all(
            isinstance(term, str) for term in terms
        ):
            raise ValueError(f"scan target {analyte} is missing string list: terms")
        rows.extend({"analyte": analyte, "term": term} for term in terms)
    return pl.DataFrame(rows, schema=SCAN_TARGET_SCHEMA)


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


def _empty_scan_results() -> pl.DataFrame:
    return pl.DataFrame(schema=SCAN_RESULT_SCHEMA)


def scan_refseq_catalog(
    db: Path,
    *,
    dataset_name: str,
    targets: pl.DataFrame,
) -> pl.DataFrame:
    """Scan a RefSeq catalog for configured analyte anchor terms."""

    frames: list[pl.DataFrame] = []
    with duckdb.connect(str(db), read_only=True) as connection:
        for row in targets.iter_rows(named=True):
            term = str(row["term"])
            pattern = f"%{term}%"
            frame = connection.execute(
                """
                select distinct
                    ? as dataset_name,
                    ? as analyte,
                    ? as term,
                    coalesce(features.protein_accession, '') as protein_accession,
                    coalesce(proteins.length_aa, 0) as length_aa,
                    coalesce(features.locus_tag, '') as locus_tag,
                    coalesce(features.gene, '') as gene,
                    coalesce(features.product, '') as product,
                    coalesce(features.seqid, '') as seqid,
                    features.start_nt,
                    features.end_nt,
                    coalesce(features.strand, '') as strand
                from features
                left join proteins using (protein_accession)
                where
                    features.protein_accession is not null
                    and (
                        features.gene ilike ?
                        or features.locus_tag ilike ?
                        or features.product ilike ?
                        or features.protein_accession ilike ?
                    )
                order by protein_accession, start_nt nulls last
                """,
                [
                    dataset_name,
                    row["analyte"],
                    term,
                    pattern,
                    pattern,
                    pattern,
                    pattern,
                ],
            ).pl()
            if not frame.is_empty():
                frames.append(frame)
    if not frames:
        return _empty_scan_results()
    return pl.concat(frames, how="vertical").unique(
        subset=["dataset_name", "analyte", "protein_accession", "locus_tag", "term"],
    )


def scan_refseq_corpus(
    manifest: Path,
    scan_config: Path,
    *,
    root: Path,
) -> pl.DataFrame:
    """Scan every RefSeq catalog for configured analyte anchor terms."""

    catalogs = read_refseq_catalog_manifest(manifest, root=root)
    targets = read_refseq_scan_config(scan_config)
    frames: list[pl.DataFrame] = []
    for row in catalogs.iter_rows(named=True):
        frame = scan_refseq_catalog(
            Path(str(row["out_db"])),
            dataset_name=str(row["dataset_name"]),
            targets=targets,
        )
        if not frame.is_empty():
            frames.append(frame)
    if not frames:
        return _empty_scan_results()
    return pl.concat(frames, how="vertical")


def _empty_anchor_hits() -> pl.DataFrame:
    return validate(pl.DataFrame(schema=ANCHOR_HITS_SCHEMA), AnchorHitsSchema)


def normalize_scan_anchor_hits(scan_results: pl.DataFrame) -> pl.DataFrame:
    """Normalize term-scan RefSeq hits into the public anchor-hit schema."""

    if scan_results.is_empty():
        return _empty_anchor_hits()
    anchor_hits = scan_results.select(
        [
            "dataset_name",
            "analyte",
            pl.col("term").alias("anchor_family"),
            "protein_accession",
            "locus_tag",
            "gene",
            "product",
            pl.lit(None, dtype=pl.Float64).alias("bitscore"),
            pl.lit(None, dtype=pl.Float64).alias("e_value"),
            pl.lit(None, dtype=pl.Float64).alias("identity"),
            pl.lit(None, dtype=pl.Float64).alias("coverage"),
            pl.lit("term_scan").alias("evidence_type"),
        ],
    ).unique(
        subset=[
            "dataset_name",
            "analyte",
            "anchor_family",
            "protein_accession",
            "locus_tag",
        ],
    )
    return validate(anchor_hits, AnchorHitsSchema)


def detect_refseq_anchor_hits(
    manifest: Path,
    scan_config: Path,
    *,
    root: Path,
    mode: str = "smoke",
    config: Path | None = None,
    profile_dir: Path = Path("data/profiles"),
    bitscore_threshold: float | None = None,
    e_value_threshold: float = 1e-20,
) -> pl.DataFrame:
    """Detect RefSeq anchors, using term scanning as the implemented smoke mode."""

    if mode == "smoke":
        return normalize_scan_anchor_hits(
            scan_refseq_corpus(manifest, scan_config, root=root),
        )
    if mode in {"profile", "hmmer"}:
        from gasregnet.search.anchors import detect_anchors_profile

        return detect_anchors_profile(
            manifest,
            config=config or Path("configs"),
            root=root,
            profile_dir=profile_dir,
            bitscore_threshold=bitscore_threshold,
            e_value_threshold=e_value_threshold,
        )
    if mode != "smoke":
        msg = (
            f"anchor detection mode {mode!r} is not implemented yet; "
            "use mode='profile' for HMMER profile discovery or mode='smoke' "
            "for term-scan anchor discovery"
        )
        raise NotImplementedError(msg)
    return _empty_anchor_hits()


def _catalog_paths(catalogs: pl.DataFrame) -> dict[str, Path]:
    return {
        str(row["dataset_name"]): Path(str(row["out_db"]))
        for row in catalogs.iter_rows(named=True)
    }


def _read_cds_features(db: Path, seqid: str) -> pl.DataFrame:
    with duckdb.connect(str(db), read_only=True) as connection:
        return connection.execute(
            """
            select
                seqid,
                source,
                feature_type,
                start_nt,
                end_nt,
                strand,
                phase,
                coalesce(feature_id, '') as feature_id,
                coalesce(parent_id, '') as parent_id,
                coalesce(locus_tag, '') as locus_tag,
                coalesce(protein_accession, '') as protein_accession,
                coalesce(gene, '') as gene,
                coalesce(product, '') as product
            from features
            where feature_type = 'CDS'
              and seqid = ?
              and protein_accession is not null
            order by start_nt, end_nt, protein_accession
            """,
            [seqid],
        ).pl()


def _anchor_feature_row(db: Path, anchor: dict[str, object]) -> dict[str, object]:
    protein_accession = str(anchor["protein_accession"])
    locus_tag = str(anchor["locus_tag"])
    gene = str(anchor["gene"])
    with duckdb.connect(str(db), read_only=True) as connection:
        row = connection.execute(
            """
            select
                seqid,
                source,
                feature_type,
                start_nt,
                end_nt,
                strand,
                phase,
                coalesce(feature_id, '') as feature_id,
                coalesce(parent_id, '') as parent_id,
                coalesce(locus_tag, '') as locus_tag,
                coalesce(protein_accession, '') as protein_accession,
                coalesce(gene, '') as gene,
                coalesce(product, '') as product
            from features
            where feature_type = 'CDS'
              and protein_accession = ?
            order by start_nt, end_nt
            limit 1
            """,
            [protein_accession],
        ).fetchone()
        if row is None and locus_tag:
            row = connection.execute(
                """
                select
                    seqid, source, feature_type, start_nt, end_nt, strand, phase,
                    coalesce(feature_id, '') as feature_id,
                    coalesce(parent_id, '') as parent_id,
                    coalesce(locus_tag, '') as locus_tag,
                    coalesce(protein_accession, '') as protein_accession,
                    coalesce(gene, '') as gene,
                    coalesce(product, '') as product
                from features
                where feature_type = 'CDS'
                  and locus_tag = ?
                order by start_nt, end_nt
                limit 1
                """,
                [locus_tag],
            ).fetchone()
        if row is None and gene:
            row = connection.execute(
                """
                select
                    seqid, source, feature_type, start_nt, end_nt, strand, phase,
                    coalesce(feature_id, '') as feature_id,
                    coalesce(parent_id, '') as parent_id,
                    coalesce(locus_tag, '') as locus_tag,
                    coalesce(protein_accession, '') as protein_accession,
                    coalesce(gene, '') as gene,
                    coalesce(product, '') as product
                from features
                where feature_type = 'CDS'
                  and gene = ?
                order by start_nt, end_nt
                limit 1
                """,
                [gene],
            ).fetchone()
    if row is None:
        raise ValueError(f"anchor feature not found: {anchor}")
    columns = list(FEATURE_SCHEMA.keys())
    return dict(zip(columns, row, strict=True))


def _feature_key(row: dict[str, object]) -> tuple[str, str, str]:
    return (
        str(row["protein_accession"]),
        str(row["locus_tag"]),
        str(row["gene"]),
    )


def _safe_identifier(value: str) -> str:
    return "".join(character if character.isalnum() else "_" for character in value)


def _empty_loci_genes() -> tuple[pl.DataFrame, pl.DataFrame]:
    return (
        validate(pl.DataFrame(schema=LOCI_OUTPUT_SCHEMA), LociSchema),
        validate(pl.DataFrame(schema=GENES_OUTPUT_SCHEMA), GenesSchema),
    )


def extract_refseq_neighborhoods(
    anchor_hits: pl.DataFrame,
    manifest: Path,
    *,
    root: Path,
    window_genes: int = 10,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Extract ±N CDS neighborhoods around RefSeq anchor hits."""

    anchor_hits = validate(anchor_hits, AnchorHitsSchema)
    if anchor_hits.is_empty():
        return _empty_loci_genes()

    catalogs = _catalog_paths(read_refseq_catalog_manifest(manifest, root=root))
    loci_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    created_at = datetime.now(UTC).replace(tzinfo=None)

    for cluster_id, anchor in enumerate(anchor_hits.iter_rows(named=True), start=1):
        dataset_name = str(anchor["dataset_name"])
        db = catalogs.get(dataset_name)
        if db is None:
            raise ValueError(f"no RefSeq catalog for dataset: {dataset_name}")
        anchor_feature = _anchor_feature_row(db, anchor)
        cds = _read_cds_features(db, str(anchor_feature["seqid"]))
        anchor_key = _feature_key(anchor_feature)
        anchor_index = None
        cds_rows = list(cds.iter_rows(named=True))
        for index, feature in enumerate(cds_rows):
            if _feature_key(feature) == anchor_key:
                anchor_index = index
                break
        if anchor_index is None:
            raise ValueError(f"anchor feature not found in contig CDS order: {anchor}")

        start_index = max(0, anchor_index - window_genes)
        end_index = min(len(cds_rows) - 1, anchor_index + window_genes)
        anchor_accession = str(anchor_feature["protein_accession"]) or str(
            anchor["protein_accession"],
        )
        anchor_family = str(anchor["anchor_family"])
        locus_id = (
            f"{str(anchor['analyte']).lower()}__{dataset_name}__"
            f"{_safe_identifier(anchor_family)}__"
            f"{anchor_accession or str(anchor_feature['locus_tag'])}"
        )
        loci_rows.append(
            {
                "locus_id": locus_id,
                "analyte": str(anchor["analyte"]),
                "anchor_accession": anchor_accession,
                "anchor_family": anchor_family,
                "organism": dataset_name,
                "taxon_id": 0,
                "cluster_id": cluster_id,
                "contig_id": str(anchor_feature["seqid"]),
                "window_size": end_index - start_index + 1,
                "is_boundary_truncated": start_index == 0
                or end_index == len(cds_rows) - 1,
                "marker_genes_present": [anchor_family],
                "accessory_genes_present": [],
                "locus_score": 0.0,
                "locus_confidence": "low",
                "taxonomic_context_score": 0.0,
                "operon_integrity_score": 0.0,
                "created_at": created_at,
            },
        )
        anchor_start = int(str(anchor_feature["start_nt"]))
        for feature_index in range(start_index, end_index + 1):
            feature = cds_rows[feature_index]
            relative_index = feature_index - anchor_index
            gene_accession = (
                str(feature["protein_accession"])
                or str(feature["feature_id"])
                or str(feature["locus_tag"])
            )
            gene_rows.append(
                {
                    "locus_id": locus_id,
                    "gene_accession": gene_accession,
                    "relative_index": relative_index,
                    "relative_start": int(str(feature["start_nt"])) - anchor_start,
                    "relative_stop": int(str(feature["end_nt"])) - anchor_start,
                    "strand": str(feature["strand"])
                    if feature["strand"] in ["+", "-"]
                    else "+",
                    "product_description": str(feature["product"]),
                    "pfam_ids": [],
                    "pfam_descriptions": [],
                    "interpro_ids": [],
                    "interpro_descriptions": [],
                    "functional_class": "anchor"
                    if relative_index == 0
                    else "unknown",
                    "regulator_class": "none",
                    "sensory_domains": [],
                    "is_anchor": relative_index == 0,
                    "is_regulator_candidate": False,
                },
            )

    loci = pl.DataFrame(
        loci_rows,
        schema_overrides=LOCI_OUTPUT_SCHEMA,
    )
    genes = pl.DataFrame(
        gene_rows,
        schema_overrides=GENES_OUTPUT_SCHEMA,
    )
    return validate(loci, LociSchema), validate(genes, GenesSchema)
