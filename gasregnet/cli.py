"""Command-line interface for GasRegNet."""

from pathlib import Path

import click

from gasregnet.config import load_config, resolve_and_dump
from gasregnet.io.parquet import write_parquet
from gasregnet.io.sqlite_efi import read_efi_sqlite
from gasregnet.logging import configure_logging
from gasregnet.manifest import build_manifest, write_manifest
from gasregnet.paths import ensure_out_dir
from gasregnet.schemas import GenesSchema, LociSchema
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
def report_command() -> None:
    """Build publication tables, figures, and captions."""

    _placeholder("report")


@app.command("repro", help="Run the headline reproducibility workflow.")
def repro_command() -> None:
    """Run the headline reproducibility workflow."""

    _placeholder("repro")
