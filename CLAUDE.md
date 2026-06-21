# CLAUDE.md

Working notes for Claude Code on this repo. This file is **gitignored** — it's a scratchpad for progress/context, not committed project docs.

## Project

WFU thesis (STA 720, Cademartori group) on **Bayesian Vector Autoregression (VAR)** models with
**Minnesota / Litterman shrinkage priors**. Current work is methodological — fitting VAR(2) models
to *simulated* data to validate the estimators before applying to real bike-share demand data.

Intended applied datasets (not yet ingested):
- `data/raw/indego/` — Indego bike share (Philadelphia)
- `data/raw/cta/` — CTA (Chicago Transit Authority)

## Layout

- `notebooks/` — **the active work lives here**
  - `var-2.ipynb` — Python/cmdstanpy driver: simulates a stable VAR(2), fits the two Stan models, checks recovery of true params.
  - `var-2.stan` — flat-prior VAR(2) (Normal(0,10) on coeffs, LKJ(2) on corr, half-Normal on sigma).
  - `var-2-hierarchical.stan` — VAR(P) with **Minnesota prior** sampled hierarchically (`lambda` tightness ~ Gamma(2,10), `theta` cross-eq ~ Beta(2,2), `alpha_lag=1.0` fixed in transformed data).
  - `var-2-conjugate-litterman.qmd` — R/Quarto: **conjugate** Minnesota posterior (closed-form, fixed innovation sigmas from univariate AR residuals). Includes `conj_litter()` helper + a sweep over `lambda` plotting posterior coefficients.
  - `simdata/Y_sim.csv` — 200×2 simulated VAR(2) series (gitignored). `Y_sim.csv` also sits in `notebooks/` root.
- `configs/sweeps/` — sweep CSVs (e.g. `lambda_rate.csv`: `idx,lambda_shape,lambda_rate`). `configs/base.yaml` still empty.
- `models/var-2-hierarchical.stan` — **cluster/sweep copy** of the notebook hierarchical model, but the gamma hyperprior on `lambda` is parameterized as data (`lambda_shape`, `lambda_rate`). `shape=2, rate=10` reproduces the notebook's `gamma(2,10)`. (Two copies exist on purpose — notebook one is interactive; consolidate later if desired.)
- `run.py` — single-spec runner: loads data, fits the model for one `(shape, rate)`, writes `summary.csv`, `lambda_theta.csv`, `metadata.json` (divergences, max R-hat, timing) to `--outdir`. `--compile-only` to precompile. Reads `CMDSTAN` env (falls back to DEAC path).
- `generate_sweep.py` — writes a sweep CSV as the Cartesian product of `--shapes` × `--rates`. Default: shape=2, rates {2.5,5,10,20,40} → prior means {0.8,0.4,0.2,0.1,0.05}.
- `slurm/sweep.sbatch` — array job; each task picks its CSV row by `$SLURM_ARRAY_TASK_ID`, runs `run.py`. Do **not** hardcode `--array` (passed at submit time). 4 cpus / 4 chains, 4GB, 1h.
- `slurm/submit.sh` — convenience: generates sweep if missing, submits a one-shot compile job, then the array `--dependency=afterok` on it (avoids tasks racing to compile the same Stan exe).
- `cleaning/`, `results/`, `data/` — scaffolding, mostly `.gitkeep` only.
- `README.md` — empty placeholder.

### Lambda sweep workflow (cluster)
```
python generate_sweep.py                 # -> configs/sweeps/lambda_rate.csv
bash slurm/submit.sh                      # compile + submit array; results/lambda_rate/idx*/
```
Each run dir gets `summary.csv`, `lambda_theta.csv` (posterior draws of the swept hyperparams), `metadata.json`.

## Environment

- Python deps in `requirements.txt`: numpy, pandas, cmdstanpy 1.3.0, matplotlib, seaborn.
- **Stan via cmdstan 2.37.0 on the DEAC cluster**: `set_cmdstan_path("/deac/sta/classes/sta720/software/cmdstan/2.37.0")` — hardcoded in the notebook; won't run Stan locally without cmdstan installed.
- R/Quarto for the `.qmd` (uses tidyverse, ggplot2).

## Simulated DGP (in var-2.ipynb)

```
c   = [2.5, 1.5]
A_1 = [[0.45,-0.35],[-0.15,0.30]]
A_2 = [[0.05,-0.05],[-0.05,0.03]]
Sigma = [[1,0.1],[0.1,1]],  T=200, seed=1   (stable: companion eigenvalues < 1)
```

## Conventions / gotchas

- Stan models index lags `1:P`; data loops run `t in (P+1):T`.
- Minnesota SD: diagonal `lambda / ell^alpha`, off-diagonal `(lambda*theta*sigma_i)/(ell^alpha * sigma_j)`.
- The `.qmd` reads `notebooks/simdata/Y_sim.csv` with a path relative to the **project root** (Quarto run from root).
- Data, SLURM logs, results, compiled Stan binaries, and `*.html` are gitignored (see `.gitignore`).

## Progress log

- 2026-06-15 — Initial repo review; created this CLAUDE.md. Repo state: simulation + 3 model
  implementations (flat Stan, hierarchical Minnesota Stan, conjugate Litterman R) all in place;
  validated on simulated data. Real bike-share data not yet ingested; run/sweep infra still stubs.
- 2026-06-15 — Built the lambda gamma-prior SLURM sweep system: parameterized
  `models/var-2-hierarchical.stan` (gamma `(shape,rate)` as data), `run.py` runner,
  `generate_sweep.py`, `slurm/sweep.sbatch` array job, `slurm/submit.sh` (compile→array w/ dependency).
  Default sweep fixes shape=2, varies rate. py_compile + bash -n + awk row-parse verified locally;
  Stan/cmdstanpy not runnable locally (no numpy/cmdstan) — must run on DEAC.
