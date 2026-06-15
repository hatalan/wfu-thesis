#!/usr/bin/env python3
"""Fit the hierarchical Minnesota VAR for a single gamma-prior specification.

One invocation == one point in a sweep. The SLURM array (slurm/sweep.sbatch)
calls this once per (lambda_shape, lambda_rate) row of a sweep CSV.

Outputs written to --outdir:
  metadata.json       run config + diagnostics (divergences, max R-hat, timing)
  summary.csv         cmdstanpy posterior summary for all parameters
  lambda_theta.csv    posterior draws of lambda and theta (the swept hyperparams)
  *.csv (optional)    full cmdstan draws, only with --save-draws

Example:
  python run.py --rate 10 --outdir results/lambda_rate/idx2_rate10
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Default cmdstan install on the DEAC cluster (matches notebooks/var-2.ipynb).
DEFAULT_CMDSTAN = "/deac/sta/classes/sta720/software/cmdstan/2.37.0"
DEFAULT_MODEL = Path(__file__).parent / "models" / "var-2-hierarchical.stan"
DEFAULT_DATA = Path(__file__).parent / "notebooks" / "simdata" / "Y_sim.csv"


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", default=str(DEFAULT_DATA),
                   help="CSV of the (T x K) series. Default: simulated Y_sim.csv")
    p.add_argument("--has-header", action="store_true",
                   help="Set if the data CSV has a header row (Y_sim.csv does not).")
    p.add_argument("--lags", type=int, default=2, help="VAR order P (default 2)")
    p.add_argument("--shape", type=float, default=2.0,
                   help="Gamma prior shape for lambda (default 2.0, held fixed in the sweep)")
    p.add_argument("--rate", type=float, required=True,
                   help="Gamma prior rate for lambda (the swept quantity)")
    p.add_argument("--model", default=str(DEFAULT_MODEL), help="Path to the .stan file")
    p.add_argument("--outdir", required=True, help="Directory for this run's outputs")
    p.add_argument("--cmdstan", default=os.environ.get("CMDSTAN", DEFAULT_CMDSTAN),
                   help="cmdstan install dir (env CMDSTAN overrides; falls back to DEAC path)")
    p.add_argument("--chains", type=int, default=4)
    p.add_argument("--iter-warmup", type=int, default=1000)
    p.add_argument("--iter-sampling", type=int, default=1000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--save-draws", action="store_true",
                   help="Also copy the full cmdstan CSV draws into --outdir")
    p.add_argument("--compile-only", action="store_true",
                   help="Compile the Stan model and exit (run once before submitting an array)")
    return p.parse_args(argv)


def setup_cmdstan(path):
    from cmdstanpy import set_cmdstan_path
    if path and os.path.isdir(path):
        set_cmdstan_path(path)
        print(f"[run] cmdstan path: {path}", flush=True)
    else:
        print(f"[run] cmdstan path '{path}' not found; using cmdstanpy default", flush=True)


def main(argv=None):
    args = parse_args(argv)
    from cmdstanpy import CmdStanModel

    setup_cmdstan(args.cmdstan)

    if args.compile_only:
        model = CmdStanModel(stan_file=args.model)  # triggers compilation
        print(f"[run] compiled {args.model} -> {model.exe_file}", flush=True)
        return 0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Load data ---
    header = 0 if args.has_header else None
    Y = pd.read_csv(args.data, header=header).to_numpy(dtype=float)
    T, K = Y.shape
    P = args.lags
    if T <= P:
        raise ValueError(f"Need T > P; got T={T}, P={P}")

    stan_data = {
        "T": T, "K": K, "P": P, "Y": Y,
        "lambda_shape": args.shape, "lambda_rate": args.rate,
    }
    prior_mean = args.shape / args.rate

    print(f"[run] data={args.data} T={T} K={K} P={P} | "
          f"lambda ~ gamma(shape={args.shape}, rate={args.rate}) "
          f"(prior mean {prior_mean:.4g})", flush=True)

    # --- Fit ---
    model = CmdStanModel(stan_file=args.model)
    t0 = time.time()
    fit = model.sample(
        data=stan_data,
        chains=args.chains,
        parallel_chains=args.chains,
        iter_warmup=args.iter_warmup,
        iter_sampling=args.iter_sampling,
        seed=args.seed,
        show_progress=False,
    )
    elapsed = time.time() - t0

    # --- Persist results ---
    summary = fit.summary()
    summary.to_csv(outdir / "summary.csv")
    fit.draws_pd(vars=["lambda", "theta"]).to_csv(outdir / "lambda_theta.csv", index=False)
    if args.save_draws:
        fit.save_csvfiles(str(outdir))

    # --- Diagnostics (best-effort; never lose a completed fit over these) ---
    diagnostics = {}
    try:
        diagnostics["num_divergent"] = int(fit.method_variables()["divergent__"].sum())
    except Exception as e:  # noqa: BLE001
        diagnostics["num_divergent"] = None
        diagnostics["divergent_error"] = str(e)
    try:
        diagnostics["max_rhat"] = float(np.nanmax(summary["R_hat"].to_numpy()))
    except Exception as e:  # noqa: BLE001
        diagnostics["max_rhat"] = None
        diagnostics["rhat_error"] = str(e)
    try:
        diagnostics["diagnose"] = fit.diagnose()
    except Exception as e:  # noqa: BLE001
        diagnostics["diagnose"] = f"diagnose() failed: {e}"

    metadata = {
        "data": os.path.abspath(args.data),
        "model": os.path.abspath(args.model),
        "T": T, "K": K, "P": P,
        "lambda_shape": args.shape,
        "lambda_rate": args.rate,
        "lambda_prior_mean": prior_mean,
        "chains": args.chains,
        "iter_warmup": args.iter_warmup,
        "iter_sampling": args.iter_sampling,
        "seed": args.seed,
        "sampling_seconds": round(elapsed, 2),
        "lambda_post_mean": float(summary.loc["lambda", "Mean"]),
        "theta_post_mean": float(summary.loc["theta", "Mean"]),
        **diagnostics,
    }
    with open(outdir / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"[run] done in {elapsed:.1f}s | lambda post mean "
          f"{metadata['lambda_post_mean']:.4g} | divergences "
          f"{diagnostics['num_divergent']} | max R-hat {diagnostics['max_rhat']} | "
          f"-> {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
