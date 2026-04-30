"""Command-line interface for GasRegNet."""

from __future__ import annotations

import typer

app = typer.Typer(help="GasRegNet comparative genomics workflows.")


@app.callback()
def main() -> None:
    """Run GasRegNet commands."""
