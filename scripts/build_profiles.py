"""Build pyhmmer HMM profiles from configured anchor seed FASTA files."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import polars as pl
import yaml  # type: ignore[import-untyped]
from pyhmmer import easel, plan7

from gasregnet.config import AnalyteConfig, load_config
from gasregnet.hashing import file_sha256
from gasregnet.io.fasta import read_fasta
from gasregnet.paths import ensure_out_dir

PROFILE_MANIFEST_SCHEMA: dict[str, Any] = {
    "analyte": pl.Utf8,
    "anchor_family": pl.Utf8,
    "profile": pl.Utf8,
    "seed_faa": pl.Utf8,
    "seed_accession": pl.Utf8,
    "selection_rule": pl.Utf8,
    "sha256": pl.Utf8,
}


def _gene_symbol(description: str) -> str:
    match = re.search(r"\bGN=([A-Za-z0-9_.-]+)", description)
    return match.group(1) if match else ""


def _select_seed(
    analyte: AnalyteConfig,
    anchor_family: str,
) -> tuple[str, str, str, str]:
    records = list(read_fasta(analyte.anchor_seeds))
    family_lower = anchor_family.lower()
    for accession, description, sequence in records:
        if _gene_symbol(description).lower() == family_lower:
            return accession, description, sequence, "gene_symbol"
    for accession, description, sequence in records:
        searchable = f"{accession} {description}".lower()
        if family_lower in searchable:
            return accession, description, sequence, "description_match"
    accession, description, sequence = records[0]
    return accession, description, sequence, "analyte_seed_fallback"


def _build_single_sequence_hmm(
    *,
    name: str,
    sequence: str,
    out_hmm: Path,
) -> None:
    alphabet = easel.Alphabet.amino()
    text_sequence = easel.TextSequence(
        name=name.encode("utf-8"),
        sequence=sequence.replace("*", ""),
    )
    digital_sequence = text_sequence.digitize(alphabet)
    builder = plan7.Builder(alphabet)
    background = plan7.Background(alphabet)
    hmm, _, _ = builder.build(digital_sequence, background)
    hmm.name = name.encode("utf-8")
    out_hmm.parent.mkdir(parents=True, exist_ok=True)
    with out_hmm.open("wb") as handle:
        hmm.write(handle)


def build_profiles(
    *,
    config: Path,
    out_dir: Path,
    manifest_out: Path,
) -> pl.DataFrame:
    """Build one HMM profile per configured anchor family."""

    loaded = load_config(config)
    ensure_out_dir(out_dir)
    rows: list[dict[str, object]] = []
    for analyte in loaded.analytes:
        for family in analyte.anchor_families:
            accession, _description, sequence, selection_rule = _select_seed(
                analyte,
                family.name,
            )
            out_hmm = out_dir / f"{family.name}.hmm"
            _build_single_sequence_hmm(
                name=family.name,
                sequence=sequence,
                out_hmm=out_hmm,
            )
            rows.append(
                {
                    "analyte": analyte.analyte,
                    "anchor_family": family.name,
                    "profile": str(out_hmm),
                    "seed_faa": str(analyte.anchor_seeds),
                    "seed_accession": accession,
                    "selection_rule": selection_rule,
                    "sha256": file_sha256(out_hmm),
                },
            )

    manifest = pl.DataFrame(rows, schema=PROFILE_MANIFEST_SCHEMA)
    manifest_out.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_csv(manifest_out.with_suffix(".csv"))
    manifest_out.write_text(
        yaml.safe_dump({"profiles": rows}, sort_keys=False),
        encoding="utf-8",
    )
    return manifest


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs"))
    parser.add_argument("--out-dir", type=Path, default=Path("data/profiles"))
    parser.add_argument(
        "--manifest-out",
        type=Path,
        default=Path("data/profiles/profiles.yaml"),
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    build_profiles(
        config=args.config,
        out_dir=args.out_dir,
        manifest_out=args.manifest_out,
    )
    print(args.manifest_out)


if __name__ == "__main__":
    main()
