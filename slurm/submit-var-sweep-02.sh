#!/bin/bash
# Submit VAR lag sweep 02: a setup job, then the array once setup succeeds.
# Run from the repository root:  bash slurm/submit-var-sweep-02.sh
set -euo pipefail

RSCRIPT=/deac/opt/rocky9-noarch/r/4.3.3/lib64/R/bin/Rscript
MAX_RUNNING=16    # array concurrency cap

# Array size comes from the grid in var-sweep-02.R, so it cannot drift out of sync.
N_TASKS=$($RSCRIPT slurm/var-sweep-02.R ntasks 2>/dev/null | tr -d '[:space:]')
if ! [[ "$N_TASKS" =~ ^[0-9]+$ ]]; then
  echo "Could not read the task count; run: $RSCRIPT slurm/var-sweep-02.R ntasks" >&2
  exit 1
fi

setup_id=$(sbatch --parsable \
  --job-name=var-sweep-02-setup \
  --cpus-per-task=1 --time=00:30:00 \
  --output=slurm/logs/%x_%j.out \
  slurm/var-sweep-02.sbatch)

array_id=$(sbatch --parsable \
  --dependency=afterok:"$setup_id" \
  --array=1-"$N_TASKS"%"$MAX_RUNNING" \
  slurm/var-sweep-02.sbatch)

echo "tasks: $N_TASKS   setup: $setup_id   array: $array_id"
