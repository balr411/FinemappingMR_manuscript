#!/usr/bin/env bash
# One-time setup: conda environment, R packages (renv), and MR-link-2.
set -euo pipefail
cd "$(dirname "$0")"

ENV_NAME=finemappingMR_manuscript
MRLINK2_COMMIT=1ee3ff1ab59e1f851eeb4dd23abc77f13cf10aa4

# 1. conda environment (mamba if available)
SOLVER=$(command -v mamba || command -v conda)
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    "$SOLVER" env update -n "$ENV_NAME" -f environment.yml --prune
else
    "$SOLVER" env create -f environment.yml
fi

# 2. MR-link-2 at a pinned commit
if [ ! -d external/mrlink2 ]; then
    mkdir -p external
    git clone https://github.com/adriaan-vd-graaf/mrlink2.git external/mrlink2
fi
git -C external/mrlink2 checkout "$MRLINK2_COMMIT"

# 3. R packages from renv.lock
conda run -n "$ENV_NAME" --no-capture-output Rscript -e 'renv::restore(prompt = FALSE)'

echo "Done. Activate with: conda activate $ENV_NAME"
