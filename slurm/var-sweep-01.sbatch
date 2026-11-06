#!/bin/bash
# One job script for both steps of the VAR lag sweep:
#   - without an array index it runs setup (simulate data, compile models);
#   - as an array task it fits one (model, p, mode) grid cell.
# Submit with slurm/submit-var-sweep-01.sh from the repository root.

#SBATCH --job-name=var-sweep-01
#SBATCH --account=cademartorigrp            # DEAC account
#SBATCH --partition=small          # DEAC partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         # one per chain
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --mail-user=hatan25@wfu.edu
#SBATCH --output=slurm/logs/%x_%A_%a.out

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module load R                     # TODO: exact module name from `module avail R`
export CMDSTAN=/deac/sta/classes/sta720/software/cmdstan/2.37.0

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  Rscript slurm/var-sweep-01.R setup
else
  Rscript slurm/var-sweep-01.R task "$SLURM_ARRAY_TASK_ID"
fi
