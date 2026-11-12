#!/bin/bash
# One job script for both steps of VAR sweep 03 (number of variables):
#   - without an array index it runs setup (simulate data, compile models);
#   - as an array task it fits one (rep, mode, N, p, model) grid cell.
# Submit with slurm/submit-var-sweep-03.sh from the repository root.

#SBATCH --job-name=var-sweep-03
#SBATCH --account=cademartorigrp            # DEAC account
#SBATCH --partition=small          # DEAC partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         # one per chain
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-user=hatan25@wfu.edu
#SBATCH --output=slurm/logs/%x_%A_%a.out

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

RSCRIPT=/deac/opt/rocky9-noarch/r/4.3.3/lib64/R/bin/Rscript
export CMDSTAN=/deac/sta/classes/sta720/software/cmdstan/2.37.0

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  $RSCRIPT slurm/var-sweep-03.R setup
else
  $RSCRIPT slurm/var-sweep-03.R task "$SLURM_ARRAY_TASK_ID"
fi
