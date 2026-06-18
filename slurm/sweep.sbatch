#!/bin/bash
#SBATCH --job-name=var2-lambda-sweep
#SBATCH --partition=small
#SBATCH --account=cademartorigrp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4GB
#SBATCH --time=00-01:00:00
#SBATCH --mail-user=hatan25@wfu.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm/logs/%x-%A_%a.out
#SBATCH --error=slurm/logs/%x-%A_%a.err
# NOTE: do not hardcode --array here; pass it at submit time so it matches the
# sweep CSV, e.g.  sbatch --array=0-4 --export=ALL,SWEEP=configs/sweeps/lambda_rate.csv slurm/sweep.sbatch
# (slurm/submit.sh computes the range for you.)

set -euo pipefail

# Environment
export CMDSTAN=${CMDSTAN:-/deac/sta/classes/sta720/software/cmdstan/2.37.0}
module load python/3.14.5 
source .venv/bin/activate

SWEEP=${SWEEP:-configs/sweeps/lambda_rate.csv}
DATA=${DATA:-notebooks/simdata/Y_sim.csv}
SWEEP_NAME=$(basename "$SWEEP" .csv)

# Pick job tasks by row in SLURM array
ROW=$(awk -F, -v id="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==id {print; exit}' "$SWEEP")
if [[ -z "$ROW" ]]; then
  echo "No row with idx=$SLURM_ARRAY_TASK_ID in $SWEEP" >&2
  exit 1
fi

# Extract hyperparameters
SHAPE=$(echo "$ROW" | cut -d, -f2)
RATE=$(echo "$ROW" | cut -d, -f3)

# Tasks write to results folder, including an .out log
OUTDIR="results/${SWEEP_NAME}/idx${SLURM_ARRAY_TASK_ID}_rate${RATE}"
echo "[sbatch] task ${SLURM_ARRAY_TASK_ID}: shape=${SHAPE} rate=${RATE} -> ${OUTDIR}"

# Hand off to run.py
python run.py \
  --data "$DATA" \
  --shape "$SHAPE" \
  --rate "$RATE" \
  --outdir "$OUTDIR" \
  --cmdstan "$CMDSTAN"
