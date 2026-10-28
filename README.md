# FinemappingMR_manuscript

Analyses for the finemappingMR manuscript, run as a Snakemake pipeline.

## Setup

```bash
./setup.sh                          # conda env + MR-link-2 + renv::restore()
conda activate finemappingMR_manuscript
```

The environment pins Snakemake 7.12.1; the Snakefile aborts if any other version is used.
GitHub-only R packages (finemappingMR, cisMRcML, gsmr2, GWASBrewer, TwoSampleMR) are pinned
by commit in `renv.lock`.

Simulation inputs (`00.inputs/`, ~1 GB: LD matrices, allele frequencies, genotype reference panels)
are not tracked; place or symlink them in the repository root.

## Running

```bash
snakemake -n                              # dry run
snakemake --profile profiles/slurm        # run on SLURM
```

## Contents

| Path | Description |
|---|---|
| `Snakefile` | Workflow |
| `scripts/` | R scripts called by the workflow |
| `config/gene_ranges.tsv` | cis regions for the 5 genes |
| `profiles/slurm/` | SLURM cluster profile |
| `01.sim_resampled/` | Output: resampled-LD simulations (untracked) |

## Simulations

1. **Resampled LD** (`01.sim_resampled`): 5 genes, 5 causal exposure variants, 0 or 2 pleiotropic outcome variants, 200 replicates.
