"""Build pyhmmer HMM profiles from configured anchor seed FASTA files."""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path
from typing import Any, cast

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
    "pfam_id": pl.Utf8,
    "profile": pl.Utf8,
    "seed_faa": pl.Utf8,
    "seed_accession": pl.Utf8,
    "selection_rule": pl.Utf8,
    "alignment": pl.Utf8,
    "nseq": pl.Int64,
    "sha256": pl.Utf8,
}
PFAM_FULL_ALIGNMENT_URL = (
    "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/"
    "{pfam_id}?annotation=alignment:full&download"
)


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


def _download_pfam_alignment(pfam_id: str, out_path: Path) -> Path:
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(
        PFAM_FULL_ALIGNMENT_URL.format(pfam_id=pfam_id),
        out_path,
    )
    return out_path


def _alignment_sequences(
    alignment_gz: Path,
    *,
    max_sequences: int,
) -> list[tuple[str, str]]:
    alphabet = easel.Alphabet.amino()
    with gzip.open(alignment_gz, "rb") as handle:
        with easel.MSAFile(
            cast(Any, handle),
            format="stockholm",
            digital=False,
            alphabet=alphabet,
        ) as msa_file:
            msa = cast(easel.TextMSA | None, msa_file.read())
    if msa is None:
        raise ValueError(f"no alignment records in {alignment_gz}")
    sequences: list[tuple[str, str]] = []
    step = max(1, len(msa.sequences) // max_sequences)
    for index in range(0, len(msa.sequences), step):
        sequence = msa.sequences[index]
        clean = re.sub(r"[^A-Za-z]", "", str(sequence.sequence)).upper()
        if clean:
            sequences.append((sequence.name.decode("utf-8"), clean))
        if len(sequences) >= max_sequences:
            break
    return sequences


def _write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")


def _build_alignment_hmm(
    *,
    name: str,
    pfam_id: str,
    out_hmm: Path,
    alignment_dir: Path,
    max_sequences: int,
) -> tuple[Path, int]:
    raw_alignment = _download_pfam_alignment(
        pfam_id,
        alignment_dir / f"{pfam_id}.full.stockholm.gz",
    )
    records = _alignment_sequences(raw_alignment, max_sequences=max_sequences)
    if len(records) < 50:
        raise ValueError(f"{pfam_id} alignment has {len(records)} sequences; need >=50")
    with tempfile.TemporaryDirectory() as tmp:
        input_fasta = Path(tmp) / f"{pfam_id}.unaligned.faa"
        aligned_fasta = Path(tmp) / f"{pfam_id}.mafft.faa"
        _write_fasta(records, input_fasta)
        with aligned_fasta.open("w", encoding="utf-8") as handle:
            subprocess.run(
                ["mafft", "--auto", "--thread", "1", str(input_fasta)],
                check=True,
                stdout=handle,
                stderr=subprocess.DEVNULL,
            )
        alphabet = easel.Alphabet.amino()
        with easel.MSAFile(
            str(aligned_fasta),
            format="afa",
            digital=True,
            alphabet=alphabet,
        ) as msa_file:
            msa = cast(easel.DigitalMSA | None, msa_file.read())
        if msa is None:
            raise ValueError(f"MAFFT produced no alignment for {pfam_id}")
        msa.name = name.encode("utf-8")
        builder = plan7.Builder(alphabet)
        background = plan7.Background(alphabet)
        hmm, _, _ = builder.build_msa(msa, background)
        hmm.name = name.encode("utf-8")
        out_hmm.parent.mkdir(parents=True, exist_ok=True)
        with out_hmm.open("wb") as handle:
            hmm.write(handle)
    return raw_alignment, len(records)


def build_profiles(
    *,
    config: Path,
    out_dir: Path,
    manifest_out: Path,
    source: str = "pfam",
    max_sequences: int = 500,
) -> pl.DataFrame:
    """Build one HMM profile per configured anchor family."""

    loaded = load_config(config)
    ensure_out_dir(out_dir)
    rows: list[dict[str, object]] = []
    for analyte in loaded.analytes:
        for family in analyte.anchor_families:
            pfam_id = family.pfam_required[0]
            accession, _description, sequence, selection_rule = _select_seed(
                analyte,
                family.name,
            )
            out_hmm = out_dir / f"{family.name}.hmm"
            alignment = ""
            nseq = 1
            if source == "pfam":
                if shutil.which("mafft") is None:
                    raise RuntimeError("MAFFT is required for Pfam profile builds")
                alignment_path, nseq = _build_alignment_hmm(
                    name=family.name,
                    pfam_id=pfam_id,
                    out_hmm=out_hmm,
                    alignment_dir=out_dir / "alignments",
                    max_sequences=max_sequences,
                )
                alignment = str(alignment_path)
                selection_rule = "pfam_full_mafft"
            elif source == "seed":
                _build_single_sequence_hmm(
                    name=family.name,
                    sequence=sequence,
                    out_hmm=out_hmm,
                )
                selection_rule = f"{selection_rule}:single_seed"
            else:
                raise ValueError(f"unknown profile source: {source}")
            rows.append(
                {
                    "analyte": analyte.analyte,
                    "anchor_family": family.name,
                    "pfam_id": pfam_id,
                    "profile": str(out_hmm),
                    "seed_faa": str(analyte.anchor_seeds),
                    "seed_accession": accession,
                    "selection_rule": selection_rule,
                    "alignment": alignment,
                    "nseq": nseq,
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
    parser.add_argument("--source", choices=["pfam", "seed"], default="pfam")
    parser.add_argument("--max-sequences", type=int, default=500)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    build_profiles(
        config=args.config,
        out_dir=args.out_dir,
        manifest_out=args.manifest_out,
        source=args.source,
        max_sequences=args.max_sequences,
    )
    print(args.manifest_out)


if __name__ == "__main__":
    main()
