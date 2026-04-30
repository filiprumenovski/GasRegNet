"""Command-line interface for GasRegNet."""

import subprocess
from pathlib import Path

import click
import polars as pl

from gasregnet.annotation.domains import annotate_domains
from gasregnet.annotation.regulators import classify_regulators
from gasregnet.archetypes.cluster import cluster_archetypes
from gasregnet.assets import Downloader, fetch_assets
from gasregnet.config import load_config, resolve_and_dump
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
    ArchetypesSchema,
    EnrichmentResultsSchema,
    GenesSchema,
    LociSchema,
    RegulatorCandidatesSchema,
)
from gasregnet.scoring.candidates import score_candidates
from gasregnet.scoring.enrichment import run_enrichment
from gasregnet.scoring.loci import score_loci
from gasregnet.search.diamond import parse_diamond_output, run_diamond


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
    if candidates.is_empty():
        output.write_text(
            "benchmark_id,analyte,protein_name,organism,hit,rank,candidate_score\n",
            encoding="utf-8",
        )
        return output
    top = candidates.sort("candidate_score", descending=True).row(0, named=True)
    output.write_text(
        "benchmark_id,analyte,protein_name,organism,hit,rank,candidate_score\n"
        f"derived_top_candidate,{top['analyte']},{top['gene_accession']},"
        f"{top['organism']},true,1,{top['candidate_score']}\n",
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
    "--out",
    required=True,
    type=click.Path(path_type=Path),
    help="Benchmark CSV path.",
)
@click.option("--seed", default=20260429, show_default=True, help="Random seed.")
@click.option("--verbose", is_flag=True, help="Enable debug logs.")
def build_benchmark_command(
    out: Path,
    seed: int,
    verbose: bool,
) -> None:
    """Create the curated benchmark CSV skeleton."""

    del seed
    configure_logging(verbose=verbose)
    out.parent.mkdir(parents=True, exist_ok=True)
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
    write_parquet(
        scored_loci.select(list(LociSchema.columns.keys())),
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
    click.echo(f"wrote scored outputs to {out}")


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
    enrichment = run_enrichment(
        genes.filter(pl.col("locus_id").is_in(case_loci)),
        genes.filter(pl.col("locus_id").is_in(control_loci)),
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
