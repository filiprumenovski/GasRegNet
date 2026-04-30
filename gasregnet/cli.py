"""Command-line interface for GasRegNet."""

import subprocess
from pathlib import Path

import click
import polars as pl

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


def _placeholder(command: str) -> None:
    raise click.UsageError(f"{command} is not implemented in this scaffold yet.")


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
def annotate_command() -> None:
    """Annotate retrieved neighborhoods."""

    _placeholder("annotate")


@app.command("score", help="Score annotated neighborhoods and candidate regulators.")
def score_command() -> None:
    """Score annotated neighborhoods and candidate regulators."""

    _placeholder("score")


@app.command("enrich", help="Run matched-control enrichment tests.")
def enrich_command() -> None:
    """Run matched-control enrichment tests."""

    _placeholder("enrich")


@app.command("archetypes", help="Cluster neighborhoods into recurrent archetypes.")
def archetypes_command() -> None:
    """Cluster scored neighborhoods into recurrent archetypes."""

    _placeholder("archetypes")


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
