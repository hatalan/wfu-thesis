#!/usr/bin/env python3
"""Generate a sweep CSV of gamma-prior specifications for lambda.

Each row is one SLURM array task. Columns: idx, lambda_shape, lambda_rate.
Rows are the Cartesian product of --shapes x --rates, so the default
(one shape, several rates) gives the "fix shape, vary rate" tightness sweep:
lambda ~ gamma(shape, rate) has prior mean shape/rate.

Examples:
  # default: shape=2, rates -> prior means {0.8, 0.4, 0.2, 0.1, 0.05}
  python generate_sweep.py

  # full grid over shape and rate
  python generate_sweep.py --name grid --shapes 1 2 4 --rates 5 10 20
"""
import argparse
import csv
import itertools
from pathlib import Path

DEFAULT_OUTDIR = Path(__file__).parent / "configs" / "sweeps"


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--name", default="lambda_rate", help="Sweep name -> <name>.csv")
    p.add_argument("--shapes", type=float, nargs="+", default=[2.0],
                   help="Gamma shape value(s). Default: 2.0 (held fixed).")
    p.add_argument("--rates", type=float, nargs="+", default=[2.5, 5.0, 10.0, 20.0, 40.0],
                   help="Gamma rate value(s). Default sweep.")
    p.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{args.name}.csv"

    rows = list(itertools.product(args.shapes, args.rates))
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["idx", "lambda_shape", "lambda_rate"])
        for idx, (shape, rate) in enumerate(rows):
            w.writerow([idx, shape, rate])

    n = len(rows)
    print(f"Wrote {n} specs -> {path}")
    print(f"SLURM array range: --array=0-{n - 1}")
    print("  idx  shape   rate   prior_mean(shape/rate)")
    for idx, (shape, rate) in enumerate(rows):
        print(f"  {idx:>3}  {shape:>5}  {rate:>5}   {shape / rate:.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
