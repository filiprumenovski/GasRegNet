# Workflow profiles

GasRegNet ships two Snakemake profiles for corpus discovery runs:

```bash
uv run snakemake -s workflows/corpus_discovery.smk --profile workflows/profiles/local
uv run snakemake -s workflows/corpus_discovery.smk --profile workflows/profiles/slurm
```

The local profile is the default workstation shape. The SLURM profile records
per-rule resources for schedulers that consume Snakemake resource declarations;
site-specific executor flags can be layered on top without changing the
workflow.
