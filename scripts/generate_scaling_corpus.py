"""Generate a deterministic synthetic RefSeq-style corpus for scaling probes."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import yaml  # type: ignore[import-untyped]

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
ANCHOR_PLAN = {
    "CO": ("coxL", "carbon monoxide dehydrogenase large subunit"),
    "CN": ("cydA", "cytochrome bd ubiquinol oxidase subunit I"),
}


@dataclass(frozen=True)
class CorpusPaths:
    """Paths emitted by the scaling corpus generator."""

    root: Path
    manifest: Path
    scan_config: Path
    corpus_config: Path
    benchmark: Path
    metadata: Path


def _sequence(seed: int, length: int = 80) -> str:
    return "".join(
        AMINO_ACIDS[(seed + offset) % len(AMINO_ACIDS)]
        for offset in range(length)
    )


def _taxon(dataset_index: int) -> dict[str, str | int]:
    taxon_id = 900_000 + dataset_index
    return {
        "taxon_id": taxon_id,
        "superkingdom": "Bacteria",
        "phylum": "Pseudomonadota" if dataset_index % 2 == 0 else "Bacillota",
        "class": "Alphaproteobacteria" if dataset_index % 2 == 0 else "Bacilli",
        "order": f"ScalingOrder{dataset_index}",
        "family": f"ScalingFamily{dataset_index}",
        "genus": f"ScalingGenus{dataset_index}",
        "species": f"Scaling species {dataset_index}",
    }


def _gene_for_index(dataset_index: int, gene_index: int) -> tuple[str, str, str]:
    if gene_index % 14 == 4:
        gene, product = ANCHOR_PLAN["CO"]
    elif gene_index % 14 == 10:
        gene, product = ANCHOR_PLAN["CN"]
    elif gene_index % 7 == 1:
        gene = f"reg{gene_index}"
        product = "PAS domain transcriptional regulator"
    elif gene_index % 7 == 5:
        gene = f"hk{gene_index}"
        product = "sensor histidine kinase"
    else:
        gene = f"orf{gene_index}"
        product = "hypothetical protein"
    accession = f"SC{dataset_index:03d}_{gene_index:05d}"
    return accession, gene, product


def _write_dataset(
    *,
    root: Path,
    dataset_index: int,
    genes_per_dataset: int,
) -> dict[str, object]:
    dataset_name = f"scaling_dataset_{dataset_index:03d}"
    dataset_dir = root / "data" / "scaling_corpus" / dataset_name
    dataset_dir.mkdir(parents=True, exist_ok=True)
    protein_faa = dataset_dir / "proteins.faa"
    gff = dataset_dir / "genomic.gff"
    seqid = f"{dataset_name}_contig_1"

    with protein_faa.open("w", encoding="utf-8") as faa_handle, gff.open(
        "w",
        encoding="utf-8",
    ) as gff_handle:
        gff_handle.write("##gff-version 3\n")
        for gene_index in range(genes_per_dataset):
            accession, gene, product = _gene_for_index(dataset_index, gene_index)
            locus_tag = f"SC{dataset_index:03d}_{gene_index:05d}"
            start = 1 + gene_index * 300
            end = start + 239
            strand = "+" if gene_index % 3 else "-"
            faa_handle.write(
                f">{accession} {product}\n{_sequence(dataset_index + gene_index)}\n",
            )
            attributes = {
                "ID": f"cds-{locus_tag}",
                "locus_tag": locus_tag,
                "protein_id": accession,
                "gene": gene,
                "product": product,
            }
            encoded = ";".join(
                f"{key}={value.replace(' ', '%20')}"
                for key, value in attributes.items()
            )
            gff_handle.write(
                "\t".join(
                    [
                        seqid,
                        "GasRegNet",
                        "CDS",
                        str(start),
                        str(end),
                        ".",
                        strand,
                        "0",
                        encoded,
                    ],
                )
                + "\n",
            )

    return {
        "dataset_name": dataset_name,
        "protein_faa": str(protein_faa.relative_to(root)),
        "gff": str(gff.relative_to(root)),
        "out_db": f"databases/{dataset_name}.duckdb",
        "taxonomy": _taxon(dataset_index),
    }


def _write_yaml(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def generate_scaling_corpus(
    out_dir: Path,
    *,
    datasets: int = 3,
    genes_per_dataset: int = 32,
    window_genes: int = 4,
) -> CorpusPaths:
    """Generate FASTA/GFF assets and configs for a small scaling probe."""

    out_dir.mkdir(parents=True, exist_ok=True)
    catalogs = [
        _write_dataset(
            root=out_dir,
            dataset_index=dataset_index,
            genes_per_dataset=genes_per_dataset,
        )
        for dataset_index in range(datasets)
    ]

    configs_dir = out_dir / "configs"
    manifest = configs_dir / "refseq_catalogs.yaml"
    scan_config = configs_dir / "refseq_scan.yaml"
    corpus_config = configs_dir / "corpus_discovery.yaml"
    benchmark = out_dir / "data" / "benchmarks" / "regulators_v2.csv"
    metadata = out_dir / "metadata.json"

    _write_yaml(manifest, {"catalogs": catalogs})
    _write_yaml(
        scan_config,
        {
            "targets": [
                {"analyte": "CO", "terms": ["coxL"]},
                {"analyte": "CN", "terms": ["cydA"]},
            ],
        },
    )
    _write_yaml(
        corpus_config,
        {
            "mode": "smoke",
            "out_dir": "results/scaling_corpus",
            "catalogs": str(manifest.relative_to(out_dir)),
            "scan_config": str(scan_config.relative_to(out_dir)),
            "benchmark": str(benchmark.relative_to(out_dir)),
            "window_genes": window_genes,
            "sharding": {"strategy": "all", "shards": ["all"]},
            "seed": 20260430,
        },
    )

    benchmark.parent.mkdir(parents=True, exist_ok=True)
    benchmark_columns = [
        "benchmark_id",
        "analyte",
        "protein_name",
        "uniprot_accession",
        "organism",
        "anchor_family",
        "expected_regulator_class",
        "expected_sensory_domains",
        "sensing_evidence_class",
        "pmid",
        "notes",
        "first_publication",
        "verify_pmid",
    ]
    benchmark_row = [
        "scaling_co_anchor",
        "CO",
        "coxL",
        "",
        "scaling_dataset_000",
        "coxL",
        "one_component",
        "PAS",
        "genomic_context",
        "[]",
        "synthetic acceptance anchor",
        "synthetic",
        "false",
    ]
    benchmark.write_text(
        ",".join(benchmark_columns) + "\n" + ",".join(benchmark_row) + "\n",
        encoding="utf-8",
    )
    metadata.write_text(
        json.dumps(
            {
                "datasets": datasets,
                "genes_per_dataset": genes_per_dataset,
                "window_genes": window_genes,
                "catalogs": len(catalogs),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return CorpusPaths(
        root=out_dir,
        manifest=manifest,
        scan_config=scan_config,
        corpus_config=corpus_config,
        benchmark=benchmark,
        metadata=metadata,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--datasets", type=int, default=3)
    parser.add_argument("--genes-per-dataset", type=int, default=32)
    parser.add_argument("--window-genes", type=int, default=4)
    args = parser.parse_args()
    paths = generate_scaling_corpus(
        args.out_dir,
        datasets=args.datasets,
        genes_per_dataset=args.genes_per_dataset,
        window_genes=args.window_genes,
    )
    print(paths.metadata)


if __name__ == "__main__":
    main()
