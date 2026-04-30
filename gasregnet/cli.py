"""Command-line interface for GasRegNet."""

import subprocess
from pathlib import Path

import click
import polars as pl

from gasregnet.annotation.domains import annotate_domains
from gasregnet.annotation.regulators import classify_regulators
from gasregnet.annotation.roles import (
    assign_sensor_roles,
    build_sensor_regulator_pairs,
)
from gasregnet.archetypes.cluster import cluster_archetypes
from gasregnet.assets import Downloader, fetch_assets
from gasregnet.benchmark import (
    evaluate_benchmark,
    load_benchmark_csv,
    write_benchmark_csv,
)
from gasregnet.config import load_config, resolve_and_dump
from gasregnet.datasets.refseq import (
    detect_refseq_anchor_hits,
    extract_refseq_neighborhoods,
    index_refseq_corpus,
    index_refseq_dataset,
    query_refseq_catalog,
    query_refseq_corpus,
    scan_refseq_corpus,
    summarize_refseq_corpus,
)
from gasregnet.io.parquet import read_parquet, write_parquet
from gasregnet.io.sqlite_efi import read_efi_sqlite
from gasregnet.logging import configure_logging
from gasregnet.manifest import build_manifest, write_manifest
from gasregnet.paths import ensure_out_dir
from gasregnet.reports.captions import build_result_led_captions, write_caption_files
from gasregnet.reports.figures import (
    figure_1_workflow_and_recovery,
    figure_2_locus_landscape,
    figure_3_archetype_atlas,
    figure_4_chemistry_partition,
    figure_5_candidate_ranking,
    figure_6_structure_hypotheses,
)
from gasregnet.reports.tables import write_publication_tables
from gasregnet.schemas import (
    AnchorHitsSchema,
    ArchetypesSchema,
    EnrichmentResultsSchema,
    EnrichmentRobustnessSchema,
    GenesSchema,
    LociSchema,
    RegulatorCandidatesSchema,
    SensorRegulatorPairsSchema,
)
from gasregnet.scoring.candidates import score_candidates
from gasregnet.scoring.cooccurrence import assign_phylogenetic_profile_scores
from gasregnet.scoring.enrichment import (
    run_enrichment,
    run_enrichment_robustness,
    run_stratified_enrichment,
)
from gasregnet.scoring.loci import score_loci
from gasregnet.scoring.posterior import assign_operon_regulation_posteriors
from gasregnet.search.diamond import parse_diamond_output, run_diamond
from gasregnet.simulation.synthetic_truth import simulate_synthetic_truth_corpus


@click.group(help="GasRegNet comparative genomics workflows.")
def app() -> None:
    """Run GasRegNet commands."""


def _config_hash_paths(config: Path) -> dict[str, Path]:
    if config.is_file():
        return {"config": config}
    return {
        str(path.relative_to(config)): path
        for path in sorted(config.rglob("*.yaml"))
        if path.is_file()
    }


def _write_metadata(
    *,
    out_dir: Path,
    config: Path | None,
    command: str,
    seed: int,
) -> None:
    ensure_out_dir(out_dir)
    config_paths = _config_hash_paths(config) if config is not None else {}
    if config is not None:
        resolve_and_dump(load_config(config), out_dir)
    manifest = build_manifest(
        seed=seed,
        command=command,
        config_paths=config_paths,
    )
    write_manifest(manifest, out_dir)


def _intermediate_dir(path: Path) -> Path:
    nested = path / "intermediate"
    return nested if nested.exists() else path


def _benchmark_from_candidates(candidates_path: Path) -> Path | None:
    candidate_table = (
        candidates_path.parent.parent / "tables" / "T1_benchmark_recovery.csv"
    )
    return candidate_table if candidate_table.exists() else None


def _fallback_benchmark(candidates_path: Path) -> Path:
    candidates = read_parquet(candidates_path, RegulatorCandidatesSchema)
    output = candidates_path.parent / "benchmark_recovery.csv"
    header = (
        "benchmark_id,analyte,protein_name,organism,hit,rank,"
        "regulation_posterior,candidate_score\n"
    )
    if candidates.is_empty():
        output.write_text(header, encoding="utf-8")
        return output
    sort_column = (
        "regulation_posterior"
        if "regulation_posterior" in candidates.columns
        else "candidate_score"
    )
    top = candidates.sort(sort_column, descending=True).row(0, named=True)
    posterior = top.get("regulation_posterior")
    output.write_text(
        header + f"derived_top_candidate,{top['analyte']},{top['gene_accession']},"
        f"{top['organism']},true,1,{'' if posterior is None else posterior},"
        f"{top['candidate_score']}\n",
        encoding="utf-8",
    )
    return output


@app.command("fetch-assets", help="Fetch external assets declared in a manifest.")
@click.option(
    "--manifest",
    default=Path("configs/assets.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="YAML asset manifest.",
)
@click.option(
    "--root",
    default=Path("."),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Repository root for relative output paths.",
)
@click.option(
    "--downloader",
    type=click.Choice(["auto", "aria2", "wget", "urllib"]),
    default="auto",
    show_default=True,
    help="Download backend. auto prefers aria2c, then wget.",
)
@click.option("--force", is_flag=True, help="Refetch existing assets.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def fetch_assets_command(
    manifest: Path,
    root: Path,
    downloader: Downloader,
    force: bool,
    verbose: bool,
) -> None:
    """Fetch external assets declared in a manifest."""

    configure_logging(verbose=verbose)
    written = fetch_assets(
        manifest,
        root=root,
        downloader=downloader,
        force=force,
    )
    for path in written:
        click.echo(path)


@app.command("build-profiles", help="Build HMM profiles from anchor seed FASTAs.")
@click.option(
    "--config",
    default=Path("configs"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="GasRegNet config directory or headline config.",
)
@click.option(
    "--out-dir",
    default=Path("data/profiles"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Profile output directory.",
)
@click.option(
    "--manifest-out",
    default=Path("data/profiles/profiles.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Profile manifest YAML output.",
)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def build_profiles_command(
    config: Path,
    out_dir: Path,
    manifest_out: Path,
    verbose: bool,
) -> None:
    """Build HMM profiles from configured anchor seed FASTAs."""

    configure_logging(verbose=verbose)
    from scripts.build_profiles import build_profiles

    manifest = build_profiles(
        config=config,
        out_dir=out_dir,
        manifest_out=manifest_out,
    )
    click.echo(f"wrote {manifest.height} profiles to {out_dir}")
    click.echo(manifest_out)


@app.command("index-refseq", help="Index RefSeq FASTA/GFF assets into DuckDB.")
@click.option(
    "--protein-faa",
    required=True,
    type=click.Path(path_type=Path),
    help="Protein FASTA, optionally gzip-compressed.",
)
@click.option(
    "--gff",
    required=True,
    type=click.Path(path_type=Path),
    help="GFF3 annotation, optionally gzip-compressed.",
)
@click.option(
    "--out-db",
    required=True,
    type=click.Path(path_type=Path),
    help="DuckDB database to create.",
)
@click.option(
    "--dataset-name",
    required=True,
    help="Stable dataset name stored in metadata.",
)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def index_refseq_command(
    protein_faa: Path,
    gff: Path,
    out_db: Path,
    dataset_name: str,
    verbose: bool,
) -> None:
    """Index RefSeq protein and annotation assets into DuckDB."""

    configure_logging(verbose=verbose)
    output = index_refseq_dataset(
        protein_faa=protein_faa,
        gff=gff,
        out_db=out_db,
        dataset_name=dataset_name,
    )
    click.echo(output)


@app.command("index-refseq-corpus", help="Index all RefSeq catalogs in a manifest.")
@click.option(
    "--manifest",
    default=Path("configs/refseq_catalogs.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="RefSeq catalog YAML manifest.",
)
@click.option(
    "--root",
    default=Path("."),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Repository root for relative manifest paths.",
)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def index_refseq_corpus_command(
    manifest: Path,
    root: Path,
    verbose: bool,
) -> None:
    """Index all RefSeq catalogs declared in a manifest."""

    configure_logging(verbose=verbose)
    for output in index_refseq_corpus(manifest, root=root):
        click.echo(output)


@app.command("query-refseq", help="Search a DuckDB RefSeq reference catalog.")
@click.option(
    "--db",
    required=True,
    type=click.Path(path_type=Path),
    help="DuckDB reference catalog.",
)
@click.option("--query", "query_text", required=True, help="Search text.")
@click.option("--limit", default=20, show_default=True, help="Maximum result rows.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def query_refseq_command(
    db: Path,
    query_text: str,
    limit: int,
    verbose: bool,
) -> None:
    """Search a DuckDB RefSeq reference catalog."""

    configure_logging(verbose=verbose)
    frame = query_refseq_catalog(db, query_text, limit=limit)
    click.echo(frame.write_csv())


@app.command("query-refseq-corpus", help="Search all RefSeq catalogs in a manifest.")
@click.option(
    "--manifest",
    default=Path("configs/refseq_catalogs.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="RefSeq catalog YAML manifest.",
)
@click.option(
    "--root",
    default=Path("."),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Repository root for relative manifest paths.",
)
@click.option("--query", "query_text", required=True, help="Search text.")
@click.option(
    "--limit-per-catalog",
    default=20,
    show_default=True,
    help="Maximum result rows from each catalog.",
)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def query_refseq_corpus_command(
    manifest: Path,
    root: Path,
    query_text: str,
    limit_per_catalog: int,
    verbose: bool,
) -> None:
    """Search all RefSeq catalogs declared in a manifest."""

    configure_logging(verbose=verbose)
    frame = query_refseq_corpus(
        manifest,
        query_text,
        root=root,
        limit_per_catalog=limit_per_catalog,
    )
    click.echo(frame.write_csv())


@app.command("summarize-refseq-corpus", help="Summarize RefSeq catalogs.")
@click.option(
    "--manifest",
    default=Path("configs/refseq_catalogs.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="RefSeq catalog YAML manifest.",
)
@click.option(
    "--root",
    default=Path("."),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Repository root for relative manifest paths.",
)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def summarize_refseq_corpus_command(
    manifest: Path,
    root: Path,
    verbose: bool,
) -> None:
    """Summarize all RefSeq catalogs declared in a manifest."""

    configure_logging(verbose=verbose)
    frame = summarize_refseq_corpus(manifest, root=root)
    click.echo(frame.write_csv())


@app.command("scan-refseq-corpus", help="Scan RefSeq catalogs for gas anchor terms.")
@click.option(
    "--manifest",
    default=Path("configs/refseq_catalogs.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="RefSeq catalog YAML manifest.",
)
@click.option(
    "--scan-config",
    default=Path("configs/refseq_scan.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Analyte scan target YAML.",
)
@click.option(
    "--root",
    default=Path("."),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Repository root for relative manifest paths.",
)
@click.option(
    "--out",
    type=click.Path(path_type=Path),
    help="Optional CSV output path. Writes to stdout when omitted.",
)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def scan_refseq_corpus_command(
    manifest: Path,
    scan_config: Path,
    root: Path,
    out: Path | None,
    verbose: bool,
) -> None:
    """Scan all RefSeq catalogs for configured analyte anchor terms."""

    configure_logging(verbose=verbose)
    frame = scan_refseq_corpus(manifest, scan_config, root=root)
    if out is None:
        click.echo(frame.write_csv())
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.write_csv(out)
    click.echo(out)


@app.command("detect-anchors", help="Detect corpus anchor hits.")
@click.option(
    "--manifest",
    default=Path("configs/refseq_catalogs.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="RefSeq catalog YAML manifest.",
)
@click.option(
    "--scan-config",
    default=Path("configs/refseq_scan.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Analyte scan target YAML for smoke mode.",
)
@click.option(
    "--mode",
    default="smoke",
    show_default=True,
    type=click.Choice(["smoke", "profile", "diamond", "hmmer"]),
    help="Anchor detection mode.",
)
@click.option(
    "--config",
    default=Path("configs"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="GasRegNet config directory or headline config.",
)
@click.option(
    "--profile-dir",
    default=Path("data/profiles"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Directory containing anchor-family HMM profiles.",
)
@click.option(
    "--bitscore-threshold",
    type=float,
    help="Optional minimum HMMER bitscore.",
)
@click.option(
    "--e-value-threshold",
    default=1e-20,
    show_default=True,
    type=float,
    help="Maximum HMMER E-value.",
)
@click.option(
    "--root",
    default=Path("."),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Repository root for relative manifest paths.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Anchor hits Parquet output path.",
)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def detect_anchors_command(
    manifest: Path,
    scan_config: Path,
    mode: str,
    config: Path,
    profile_dir: Path,
    bitscore_threshold: float | None,
    e_value_threshold: float,
    root: Path,
    out: Path,
    verbose: bool,
) -> None:
    """Detect anchors and write normalized anchor-hit Parquet."""

    configure_logging(verbose=verbose)
    try:
        anchor_hits = detect_refseq_anchor_hits(
            manifest,
            scan_config,
            root=root,
            mode=mode,
            config=config,
            profile_dir=profile_dir,
            bitscore_threshold=bitscore_threshold,
            e_value_threshold=e_value_threshold,
        )
    except NotImplementedError as exc:
        raise click.ClickException(str(exc)) from exc
    write_parquet(anchor_hits, out, AnchorHitsSchema)
    click.echo(out)


@app.command("detect-anchors-profile", help="Detect anchors with HMM profiles.")
@click.option(
    "--manifest",
    default=Path("configs/refseq_catalogs.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="RefSeq catalog YAML manifest.",
)
@click.option(
    "--config",
    default=Path("configs"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="GasRegNet config directory or headline config.",
)
@click.option(
    "--profile-dir",
    default=Path("data/profiles"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Directory containing anchor-family HMM profiles.",
)
@click.option(
    "--bitscore-threshold",
    type=float,
    help="Optional minimum HMMER bitscore.",
)
@click.option(
    "--e-value-threshold",
    default=1e-20,
    show_default=True,
    type=float,
    help="Maximum HMMER E-value.",
)
@click.option(
    "--root",
    default=Path("."),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Repository root for relative manifest paths.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Anchor hits Parquet output path.",
)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def detect_anchors_profile_command(
    manifest: Path,
    config: Path,
    profile_dir: Path,
    bitscore_threshold: float | None,
    e_value_threshold: float,
    root: Path,
    out: Path,
    verbose: bool,
) -> None:
    """Detect anchors with profile HMMs and write normalized anchor hits."""

    configure_logging(verbose=verbose)
    anchor_hits = detect_refseq_anchor_hits(
        manifest,
        Path("configs/refseq_scan.yaml"),
        root=root,
        mode="profile",
        config=config,
        profile_dir=profile_dir,
        bitscore_threshold=bitscore_threshold,
        e_value_threshold=e_value_threshold,
    )
    write_parquet(anchor_hits, out, AnchorHitsSchema)
    click.echo(out)


@app.command("extract-neighborhoods", help="Extract RefSeq neighborhoods.")
@click.option(
    "--anchor-hits",
    required=True,
    type=click.Path(path_type=Path),
    help="Anchor hits Parquet from detect-anchors.",
)
@click.option(
    "--manifest",
    default=Path("configs/refseq_catalogs.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="RefSeq catalog YAML manifest.",
)
@click.option(
    "--root",
    default=Path("."),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Repository root for relative manifest paths.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory for loci.parquet and genes.parquet.",
)
@click.option(
    "--window-genes",
    default=10,
    show_default=True,
    help="Number of CDS features to include on each side of an anchor.",
)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def extract_neighborhoods_command(
    anchor_hits: Path,
    manifest: Path,
    root: Path,
    out: Path,
    window_genes: int,
    verbose: bool,
) -> None:
    """Extract canonical loci and genes from RefSeq catalogs."""

    configure_logging(verbose=verbose)
    hits = read_parquet(anchor_hits, AnchorHitsSchema)
    loci, genes = extract_refseq_neighborhoods(
        hits,
        manifest,
        root=root,
        window_genes=window_genes,
    )
    write_parquet(loci, out / "intermediate" / "loci.parquet", LociSchema)
    write_parquet(genes, out / "intermediate" / "genes.parquet", GenesSchema)
    click.echo(out)


@app.command("evaluate-benchmark", help="Evaluate benchmark recovery.")
@click.option(
    "--benchmark",
    default=Path("data/benchmarks/regulators_v2.csv"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Benchmark CSV.",
)
@click.option(
    "--anchor-hits",
    required=True,
    type=click.Path(path_type=Path),
    help="Anchor hits Parquet.",
)
@click.option(
    "--candidates",
    type=click.Path(path_type=Path),
    help="Optional scored candidates Parquet.",
)
@click.option(
    "--out",
    type=click.Path(path_type=Path),
    help="Optional recovery CSV output. Writes to stdout when omitted.",
)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def evaluate_benchmark_command(
    benchmark: Path,
    anchor_hits: Path,
    candidates: Path | None,
    out: Path | None,
    verbose: bool,
) -> None:
    """Evaluate benchmark recovery from anchors and optional candidates."""

    configure_logging(verbose=verbose)
    hits = read_parquet(anchor_hits, AnchorHitsSchema)
    candidate_frame = (
        read_parquet(candidates, RegulatorCandidatesSchema)
        if candidates is not None
        else None
    )
    recovery = evaluate_benchmark(benchmark, hits, candidate_frame)
    if out is None:
        click.echo(recovery.write_csv())
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    recovery.write_csv(out)
    click.echo(out)


@app.command("validate-config", help="Validate a composed GasRegNet configuration.")
@click.option(
    "--config",
    required=True,
    type=click.Path(path_type=Path),
    help="Config directory or YAML.",
)
@click.option(
    "--out",
    type=click.Path(path_type=Path),
    help="Optional run directory for resolved metadata.",
)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def validate_config_command(
    config: Path,
    out: Path | None,
    seed: int,
    verbose: bool,
) -> None:
    """Validate a composed GasRegNet configuration."""

    configure_logging(verbose=verbose)
    load_config(config)
    if out is not None:
        _write_metadata(
            out_dir=out,
            config=config,
            command="validate-config",
            seed=seed,
        )
    click.echo("config valid")


@app.command(
    "run-sqlite",
    help="Import an EFI-GNT SQLite export into validated Parquet tables.",
)
@click.option(
    "--sqlite",
    required=True,
    type=click.Path(path_type=Path),
    help="EFI-GNT SQLite export.",
)
@click.option("--analytes", multiple=True, required=True, help="Analytes to import.")
@click.option(
    "--config",
    required=True,
    type=click.Path(path_type=Path),
    help="Config directory or YAML.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Run output directory.",
)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def run_sqlite_command(
    sqlite: Path,
    analytes: tuple[str, ...],
    config: Path,
    out: Path,
    seed: int,
    verbose: bool,
) -> None:
    """Import an EFI-GNT SQLite export into validated Parquet tables."""

    configure_logging(verbose=verbose)
    _write_metadata(out_dir=out, config=config, command="run-sqlite", seed=seed)
    for analyte in analytes:
        loci, genes = read_efi_sqlite(sqlite, analyte)
        prefix = analyte.lower()
        write_parquet(loci, out / "intermediate" / f"{prefix}_loci.parquet", LociSchema)
        write_parquet(
            genes,
            out / "intermediate" / f"{prefix}_genes.parquet",
            GenesSchema,
        )
    click.echo(f"imported {len(analytes)} analyte(s)")


@app.command("diamond-search", help="Run DIAMOND and write parsed hits as Parquet.")
@click.option(
    "--query",
    required=True,
    type=click.Path(path_type=Path),
    help="Query FASTA.",
)
@click.option(
    "--db",
    required=True,
    type=click.Path(path_type=Path),
    help="DIAMOND database.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Output Parquet path.",
)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def diamond_search_command(
    query: Path,
    db: Path,
    out: Path,
    seed: int,
    verbose: bool,
) -> None:
    """Run DIAMOND and write parsed hits as Parquet."""

    del seed
    configure_logging(verbose=verbose)
    tsv = out.with_suffix(".tsv")
    run_diamond(query, db, tsv)
    hits = parse_diamond_output(tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    hits.write_parquet(out)
    click.echo(f"wrote {out}")


@app.command("build-benchmark", help="Create the curated benchmark CSV skeleton.")
@click.option(
    "--version",
    type=click.Choice(["v1", "v2"]),
    default="v1",
    show_default=True,
    help="Benchmark version to build.",
)
@click.option(
    "--source",
    type=click.Path(path_type=Path),
    help="Source benchmark CSV for versioned benchmark builds.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Benchmark CSV path.",
)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def build_benchmark_command(
    version: str,
    source: Path | None,
    out: Path,
    seed: int,
    verbose: bool,
) -> None:
    """Create a benchmark CSV."""

    del seed
    configure_logging(verbose=verbose)
    out.parent.mkdir(parents=True, exist_ok=True)
    if version == "v2":
        source_path = source or Path("configs/benchmarks/regulators_v2.csv")
        benchmark = load_benchmark_csv(source_path)
        write_benchmark_csv(benchmark, out)
        click.echo(f"wrote {out}")
        return
    out.write_text(
        "benchmark_id,analyte,protein_name,uniprot_accession,organism,taxon_id,"
        "anchor_family,expected_regulator_class,expected_sensory_domains,"
        "sensing_evidence_class,pmid,notes\n",
        encoding="utf-8",
    )
    click.echo(f"wrote {out}")


@app.command("annotate", help="Annotate retrieved neighborhoods.")
@click.option(
    "--neighborhoods",
    required=True,
    type=click.Path(path_type=Path),
    help="Directory containing loci.parquet and genes.parquet.",
)
@click.option(
    "--config",
    required=True,
    type=click.Path(path_type=Path),
    help="Config directory or YAML.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory.",
)
@click.option(
    "--pfam-table",
    type=click.Path(path_type=Path),
    help="Optional Pfam annotation CSV.",
)
@click.option(
    "--interpro-table",
    type=click.Path(path_type=Path),
    help="Optional InterPro annotation CSV.",
)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def annotate_command(
    neighborhoods: Path,
    config: Path,
    out: Path,
    pfam_table: Path | None,
    interpro_table: Path | None,
    seed: int,
    verbose: bool,
) -> None:
    """Annotate retrieved neighborhoods."""

    configure_logging(verbose=verbose)
    ensure_out_dir(out)
    input_dir = _intermediate_dir(neighborhoods)
    cfg = load_config(config)
    loci = read_parquet(input_dir / "loci.parquet", LociSchema)
    genes = read_parquet(input_dir / "genes.parquet", GenesSchema)
    if pfam_table is not None or interpro_table is not None:
        pfam = (
            pl.read_csv(pfam_table)
            if pfam_table is not None
            else pl.DataFrame(
                schema={
                    "gene_accession": pl.Utf8,
                    "pfam_id": pl.Utf8,
                    "pfam_description": pl.Utf8,
                },
            )
        )
        interpro = (
            pl.read_csv(interpro_table)
            if interpro_table is not None
            else pl.DataFrame(
                schema={
                    "gene_accession": pl.Utf8,
                    "interpro_id": pl.Utf8,
                    "interpro_description": pl.Utf8,
                },
            )
        )
        genes = annotate_domains(genes, pfam, interpro)
    genes = classify_regulators(genes, cfg.regulator_families)
    write_parquet(loci, out / "intermediate" / "loci.parquet", LociSchema)
    write_parquet(genes, out / "intermediate" / "genes.parquet", GenesSchema)
    resolve_and_dump(cfg, out)
    write_manifest(
        build_manifest(
            seed=seed,
            command="annotate",
            config_paths=_config_hash_paths(config),
            input_paths={
                "loci": input_dir / "loci.parquet",
                "genes": input_dir / "genes.parquet",
            },
        ),
        out,
    )
    click.echo(f"wrote annotated outputs to {out}")


@app.command("assign-roles", help="Assign V2 anchor, regulator, and sensor roles.")
@click.option(
    "--neighborhoods",
    required=True,
    type=click.Path(path_type=Path),
    help="Directory containing loci.parquet and genes.parquet.",
)
@click.option(
    "--config",
    required=True,
    type=click.Path(path_type=Path),
    help="Config directory or YAML.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory.",
)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def assign_roles_command(
    neighborhoods: Path,
    config: Path,
    out: Path,
    seed: int,
    verbose: bool,
) -> None:
    """Assign V2 gene roles and sensor-regulator pairs."""

    configure_logging(verbose=verbose)
    ensure_out_dir(out)
    input_dir = _intermediate_dir(neighborhoods)
    cfg = load_config(config)
    loci = read_parquet(input_dir / "loci.parquet", LociSchema)
    genes = read_parquet(input_dir / "genes.parquet", GenesSchema)
    assigned = assign_sensor_roles(
        genes,
        regulator_families=cfg.regulator_families,
        sensory_domain_catalog=cfg.sensory_domains,
        paired_evidence_rules=cfg.paired_evidence,
    )
    pairs = build_sensor_regulator_pairs(
        assigned,
        loci,
        sensory_domain_catalog=cfg.sensory_domains,
        paired_evidence_rules=cfg.paired_evidence,
    )
    write_parquet(loci, out / "intermediate" / "loci.parquet", LociSchema)
    write_parquet(assigned, out / "intermediate" / "genes.parquet", GenesSchema)
    write_parquet(
        pairs,
        out / "intermediate" / "sensor_regulator_pairs.parquet",
        SensorRegulatorPairsSchema,
    )
    resolve_and_dump(cfg, out)
    write_manifest(
        build_manifest(
            seed=seed,
            command="assign-roles",
            config_paths=_config_hash_paths(config),
            input_paths={
                "loci": input_dir / "loci.parquet",
                "genes": input_dir / "genes.parquet",
            },
        ),
        out,
    )
    click.echo(f"wrote role assignments to {out}")


@app.command("score", help="Score annotated neighborhoods and candidate regulators.")
@click.option(
    "--neighborhoods",
    required=True,
    type=click.Path(path_type=Path),
    help="Directory containing loci.parquet and genes.parquet.",
)
@click.option(
    "--config",
    required=True,
    type=click.Path(path_type=Path),
    help="Config directory or YAML.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory.",
)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def score_command(
    neighborhoods: Path,
    config: Path,
    out: Path,
    seed: int,
    verbose: bool,
) -> None:
    """Score annotated neighborhoods and candidate regulators."""

    configure_logging(verbose=verbose)
    ensure_out_dir(out)
    input_dir = _intermediate_dir(neighborhoods)
    cfg = load_config(config)
    loci = read_parquet(input_dir / "loci.parquet", LociSchema)
    genes = read_parquet(input_dir / "genes.parquet", GenesSchema)
    scored_loci = score_loci(loci, cfg.scoring)
    candidates = score_candidates(scored_loci, genes, cfg.scoring)
    candidates = assign_phylogenetic_profile_scores(
        candidates,
        scored_loci,
        scoring=cfg.scoring,
    )
    candidates = assign_operon_regulation_posteriors(candidates)
    locus_columns = [
        column for column in LociSchema.columns if column in scored_loci.columns
    ]
    write_parquet(
        scored_loci.select(locus_columns),
        out / "intermediate" / "loci.parquet",
        LociSchema,
    )
    write_parquet(genes, out / "intermediate" / "genes.parquet", GenesSchema)
    write_parquet(
        candidates,
        out / "intermediate" / "candidates.parquet",
        RegulatorCandidatesSchema,
    )
    resolve_and_dump(cfg, out)
    write_manifest(
        build_manifest(
            seed=seed,
            command="score",
            config_paths=_config_hash_paths(config),
            input_paths={
                "loci": input_dir / "loci.parquet",
                "genes": input_dir / "genes.parquet",
            },
        ),
        out,
    )


@app.command(
    "simulate-synthetic-truth",
    help="Generate synthetic known-truth genomes for calibration tests.",
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option("--n-genomes", default=24, show_default=True, type=int)
@click.option("--positive-fraction", default=0.5, show_default=True, type=float)
@click.option("--phylogenetic-clumping", default=0.7, show_default=True, type=float)
@click.option("--annotation-noise", default=0.1, show_default=True, type=float)
@click.option("--seed", default=20260430, show_default=True, type=int)
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def simulate_synthetic_truth_command(
    out: Path,
    n_genomes: int,
    positive_fraction: float,
    phylogenetic_clumping: float,
    annotation_noise: float,
    seed: int,
    verbose: bool,
) -> None:
    """Write loci, genes, and ground truth frames for calibration runs."""

    configure_logging(verbose=verbose)
    ensure_out_dir(out)
    corpus = simulate_synthetic_truth_corpus(
        n_genomes=n_genomes,
        positive_fraction=positive_fraction,
        phylogenetic_clumping=phylogenetic_clumping,
        annotation_noise=annotation_noise,
        seed=seed,
    )
    write_parquet(corpus.loci, out / "intermediate" / "loci.parquet", LociSchema)
    write_parquet(corpus.genes, out / "intermediate" / "genes.parquet", GenesSchema)
    ground_truth_path = out / "intermediate" / "synthetic_ground_truth.csv"
    ground_truth_path.parent.mkdir(parents=True, exist_ok=True)
    corpus.ground_truth.write_csv(ground_truth_path)
    write_manifest(
        build_manifest(
            seed=seed,
            command="simulate-synthetic-truth",
            config_paths={},
            input_paths={},
        ),
        out,
    )
    click.echo(f"wrote synthetic truth corpus to {out}")


@app.command("enrich", help="Run matched-control enrichment tests.")
@click.option(
    "--scored",
    required=True,
    type=click.Path(path_type=Path),
    help="Directory containing loci.parquet and genes.parquet.",
)
@click.option(
    "--config",
    required=True,
    type=click.Path(path_type=Path),
    help="Config directory or YAML.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory.",
)
@click.option("--case-analyte", default="CO", show_default=True)
@click.option("--control-analyte", default="CN", show_default=True)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def enrich_command(
    scored: Path,
    config: Path,
    out: Path,
    case_analyte: str,
    control_analyte: str,
    seed: int,
    verbose: bool,
) -> None:
    """Run matched-control enrichment tests."""

    configure_logging(verbose=verbose)
    ensure_out_dir(out)
    input_dir = _intermediate_dir(scored)
    cfg = load_config(config)
    loci = read_parquet(input_dir / "loci.parquet", LociSchema)
    genes = read_parquet(input_dir / "genes.parquet", GenesSchema)
    case_loci = loci.filter(pl.col("analyte") == case_analyte)["locus_id"].to_list()
    control_loci = loci.filter(pl.col("analyte") == control_analyte)[
        "locus_id"
    ].to_list()
    case_genes = genes.filter(pl.col("locus_id").is_in(case_loci))
    control_genes = genes.filter(pl.col("locus_id").is_in(control_loci))
    if cfg.scoring.enrichment.test == "cmh":
        enrichment = run_stratified_enrichment(
            case_genes,
            control_genes,
            analyte=case_analyte,
            case_definition=f"{case_analyte} scored loci",
            control_definition=f"{control_analyte} scored loci",
            stratum_column=cfg.scoring.enrichment.stratum_column,
            alpha=cfg.scoring.enrichment.alpha,
            deduplication_policy=cfg.scoring.enrichment.strict_policy,
        )
    else:
        enrichment = run_enrichment(
            case_genes,
            control_genes,
            analyte=case_analyte,
            case_definition=f"{case_analyte} scored loci",
            control_definition=f"{control_analyte} scored loci",
            alpha=cfg.scoring.enrichment.alpha,
        )
    write_parquet(
        enrichment,
        out / "intermediate" / "enrichment.parquet",
        EnrichmentResultsSchema,
    )
    enrichment_robustness = run_enrichment_robustness(
        case_genes,
        control_genes,
        analyte=case_analyte,
        case_definition=f"{case_analyte} scored loci",
        control_definition=f"{control_analyte} scored loci",
        stratum_column=cfg.scoring.enrichment.stratum_column,
        alpha=cfg.scoring.enrichment.alpha,
    )
    write_parquet(
        enrichment_robustness,
        out / "intermediate" / "enrichment_robustness.parquet",
        EnrichmentRobustnessSchema,
    )
    resolve_and_dump(cfg, out)
    write_manifest(
        build_manifest(
            seed=seed,
            command="enrich",
            config_paths=_config_hash_paths(config),
            input_paths={
                "loci": input_dir / "loci.parquet",
                "genes": input_dir / "genes.parquet",
            },
        ),
        out,
    )
    click.echo(f"wrote enrichment outputs to {out}")


@app.command("archetypes", help="Cluster neighborhoods into recurrent archetypes.")
@click.option(
    "--scored",
    required=True,
    type=click.Path(path_type=Path),
    help="Directory containing loci.parquet and candidates.parquet.",
)
@click.option(
    "--config",
    required=True,
    type=click.Path(path_type=Path),
    help="Config directory or YAML.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory.",
)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def archetypes_command(
    scored: Path,
    config: Path,
    out: Path,
    seed: int,
    verbose: bool,
) -> None:
    """Cluster scored neighborhoods into recurrent archetypes."""

    configure_logging(verbose=verbose)
    ensure_out_dir(out)
    input_dir = _intermediate_dir(scored)
    cfg = load_config(config)
    loci = read_parquet(input_dir / "loci.parquet", LociSchema)
    candidates = read_parquet(
        input_dir / "candidates.parquet",
        RegulatorCandidatesSchema,
    )
    archetypes = cluster_archetypes(loci, candidates)
    write_parquet(
        archetypes,
        out / "intermediate" / "archetypes.parquet",
        ArchetypesSchema,
    )
    resolve_and_dump(cfg, out)
    write_manifest(
        build_manifest(
            seed=seed,
            command="archetypes",
            config_paths=_config_hash_paths(config),
            input_paths={
                "loci": input_dir / "loci.parquet",
                "candidates": input_dir / "candidates.parquet",
            },
        ),
        out,
    )
    click.echo(f"wrote archetype outputs to {out}")


@app.command("report", help="Build publication tables, figures, and captions.")
@click.option(
    "--results",
    required=True,
    type=click.Path(path_type=Path),
    help="Run directory containing intermediate Parquet files.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Report output directory.",
)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def report_command(
    results: Path,
    out: Path,
    seed: int,
    verbose: bool,
) -> None:
    """Build publication tables, figures, and captions."""

    configure_logging(verbose=verbose)
    ensure_out_dir(out)
    intermediate = results / "intermediate"
    loci = read_parquet(intermediate / "loci.parquet", LociSchema)
    candidates = read_parquet(
        intermediate / "candidates.parquet",
        RegulatorCandidatesSchema,
    )
    enrichment = read_parquet(
        intermediate / "enrichment.parquet",
        EnrichmentResultsSchema,
    )
    archetypes = read_parquet(intermediate / "archetypes.parquet", ArchetypesSchema)
    benchmark_path = _benchmark_from_candidates(intermediate / "candidates.parquet")
    if benchmark_path is None:
        benchmark_path = _fallback_benchmark(intermediate / "candidates.parquet")
    benchmark_recovery = pl.read_csv(benchmark_path)

    write_publication_tables(
        benchmark_recovery=benchmark_recovery,
        candidates=candidates,
        enrichment=enrichment,
        archetypes=archetypes,
        out_dir=out / "tables",
    )
    figures_dir = out / "figures"
    figure_1_workflow_and_recovery(benchmark_recovery, figures_dir)
    figure_2_locus_landscape(loci, figures_dir)
    figure_3_archetype_atlas(archetypes, figures_dir)
    figure_4_chemistry_partition(enrichment, figures_dir)
    figure_5_candidate_ranking(candidates, figures_dir)
    figure_6_structure_hypotheses(candidates, results / "structures", figures_dir)
    captions = build_result_led_captions(
        benchmark_results=benchmark_recovery,
        loci=loci,
        archetypes=archetypes,
        enrichment=enrichment,
        candidates=candidates,
        top_candidates=candidates,
    )
    write_caption_files(captions, out / "captions")
    manifest = build_manifest(
        seed=seed,
        command="report",
        input_paths={
            "loci": intermediate / "loci.parquet",
            "candidates": intermediate / "candidates.parquet",
            "enrichment": intermediate / "enrichment.parquet",
            "archetypes": intermediate / "archetypes.parquet",
        },
    )
    write_manifest(manifest, out)
    click.echo(f"wrote report artifacts to {out}")


@app.command("repro", help="Run the headline reproducibility workflow.")
@click.option(
    "--config",
    default=Path("configs/headline.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Headline config YAML.",
)
@click.option(
    "--out",
    default=Path("results/repro"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Workflow output directory.",
)
@click.option(
    "--sqlite",
    default=Path("tests/fixtures/mini_efi.sqlite"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="EFI-GNT SQLite input for the local workflow.",
)
@click.option("--cores", default=1, show_default=True, help="Snakemake cores.")
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def repro_command(
    config: Path,
    out: Path,
    sqlite: Path,
    cores: int,
    seed: int,
    verbose: bool,
) -> None:
    """Run the headline reproducibility workflow."""

    del seed
    configure_logging(verbose=verbose)
    subprocess.run(
        [
            "uv",
            "run",
            "snakemake",
            "-s",
            "workflows/sqlite_mode.smk",
            "--cores",
            str(cores),
            "--config",
            f"out_dir={out}",
            f"config_path={config}",
            f"sqlite={sqlite}",
        ],
        check=True,
    )
    click.echo(f"repro artifacts written to {out}")


@app.command("corpus-discovery", help="Run the RefSeq corpus discovery workflow.")
@click.option(
    "--corpus-config",
    default=Path("configs/corpus_discovery.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Corpus discovery config YAML.",
)
@click.option(
    "--config",
    default=Path("configs/headline.yaml"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Headline config YAML.",
)
@click.option(
    "--out",
    default=Path("results/corpus"),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Workflow output directory.",
)
@click.option(
    "--root",
    default=Path("."),
    show_default=True,
    type=click.Path(path_type=Path),
    help="Repository root for relative manifest paths.",
)
@click.option("--cores", default=1, show_default=True, help="Snakemake cores.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def corpus_discovery_command(
    corpus_config: Path,
    config: Path,
    out: Path,
    root: Path,
    cores: int,
    verbose: bool,
) -> None:
    """Run the RefSeq corpus discovery workflow."""

    configure_logging(verbose=verbose)
    subprocess.run(
        [
            "uv",
            "run",
            "snakemake",
            "-s",
            "workflows/corpus_discovery.smk",
            "--cores",
            str(cores),
            "--config",
            f"out_dir={out}",
            f"corpus_config={corpus_config}",
            f"config_path={config}",
            f"root={root}",
        ],
        check=True,
    )
    click.echo(f"corpus discovery artifacts written to {out}")
