#!/bin/bash
# Submit the VAR lag sweep: a setup job, then the array once setup succeeds.
# Run from the repository root:  bash slurm/submit-var-sweep-01.sh
set -euo pipefail

N_TASKS=64        # 4 models x 8 lag orders x 2 modes (fit, prior); must match var-sweep-01.R
MAX_RUNNING=16    # array concurrency cap

setup_id=$(sbatch --parsable \
  --job-name=var-sweep-01-setup \
  --cpus-per-task=1 --time=00:30:00 \
  --output=slurm/logs/%x_%j.out \
  slurm/var-sweep-01.sbatch)

array_id=$(sbatch --parsable \
  --dependency=afterok:"$setup_id" \
  --array=1-"$N_TASKS"%"$MAX_RUNNING" \
  slurm/var-sweep-01.sbatch)

echo "setup: $setup_id   array: $array_id"
