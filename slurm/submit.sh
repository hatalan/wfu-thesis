#!/bin/bash
# Generate (if needed), compile once, then submit the lambda-sweep array.
#
# Compiling once up front avoids every array task racing to build the same
# Stan executable. The array job is held until the compile job succeeds.
#
# Usage:
#   slurm/submit.sh                              # default sweep configs/sweeps/lambda_rate.csv
#   slurm/submit.sh configs/sweeps/grid.csv      # a different sweep CSV

set -euo pipefail
cd "$(dirname "$0")/.."

# Set sweep .csv and cmdstan location
SWEEP=${1:-configs/sweeps/lambda_rate.csv}
export CMDSTAN=${CMDSTAN:-/deac/sta/classes/sta720/software/cmdstan/2.37.0}

# Auto-generate sweep if missing
if [[ ! -f "$SWEEP" ]]; then
  echo "[submit] sweep $SWEEP not found; generating default with generate_sweep.py"
  python generate_sweep.py --name "$(basename "$SWEEP" .csv)"
fi


# Count the tasks
N=$(($(wc -l < "$SWEEP") - 1))   # rows minus header
if (( N < 1 )); then
  echo "[submit] no spec rows in $SWEEP" >&2
  exit 1
fi
echo "[submit] $SWEEP has $N specs -> array 0-$((N-1))"

# Compile the model
COMPILE_JID=$(sbatch --parsable \
  --job-name=var2-compile --partition=small --account=cademartorigrp \
  --nodes=1 --ntasks-per-node=1 --cpus-per-task=2 --mem=2GB --time=00-00:20:00 \
  --output=slurm/logs/%x-%j.out --error=slurm/logs/%x-%j.err \
  --export=ALL,CMDSTAN="$CMDSTAN" \
  --wrap="module load apps/python/3.14.5 && source .venv/bin/activate && python run.py --rate 1 --outdir /tmp/compile --compile-only")
echo "[submit] compile job: $COMPILE_JID"

# Submit the array, held until compile succeeds
ARRAY_JID=$(sbatch --parsable \
  --dependency=afterok:"$COMPILE_JID" \
  --array=0-$((N-1)) \
  --export=ALL,SWEEP="$SWEEP",CMDSTAN="$CMDSTAN" \
  slurm/sweep.sbatch)
echo "[submit] sweep array job: $ARRAY_JID (depends on $COMPILE_JID)"
echo "[submit] watch with: squeue -u \$USER ; results land in results/$(basename "$SWEEP" .csv)/"
