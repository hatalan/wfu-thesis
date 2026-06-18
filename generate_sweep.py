#!/usr/bin/env python3

"""
Generate a sweep CSV of gamma-prior specifications for lambda.
"""

# Import packages
import argparse
import csv
import itertools
from pathlib import Path

DEFAULT_OUTDIR = Path(__file__).parent / "configs" / "sweeps"

# Defining the flags
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Generates sweep name
    p.add_argument("--name", default="lambda_rate", help="Sweep name -> <name>.csv")

    # Shape parameter of Gamma distribution (held fixed)
    p.add_argument("--shapes", type=float, nargs="+", default=[2.0],
                   help="Gamma shape value(s). Default: 2.0 (held fixed).")
    
    # Rate parameter of Gamma distribution (sweep across values)
    p.add_argument("--rates", type=float, nargs="+", default=[2.5, 5.0, 10.0, 20.0, 40.0],
                   help="Gamma rate value(s). Default sweep.")
    
    # Where to write the file (defaulting to configs/sweeps)
    p.add_argument("--outdir", default=str(DEFAULT_OUTDIR))

    # Returns a container whose attributes are parsed flags (args.name, args.shapes, sargs.rates, and arms.outdir)
    return p.parse_args(argv)


def main(argv=None):
    # Parse the arguments and generte the output path
    args = parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{args.name}.csv"

    # Generates combination rows from the shapes and rates list 
    rows = list(itertools.product(args.shapes, args.rates))

    # Write the .csv
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["idx", "lambda_shape", "lambda_rate"])
        for idx, (shape, rate) in enumerate(rows):
            w.writerow([idx, shape, rate])

    # Write the console report (how many specs were made, and the exact array range used, etc.)
    n = len(rows)
    print(f"Wrote {n} specs -> {path}")
    print(f"SLURM array range: --array=0-{n - 1}")
    print("  idx  shape   rate   prior_mean(shape/rate)")
    for idx, (shape, rate) in enumerate(rows):
        print(f"  {idx:>3}  {shape:>5}  {rate:>5}   {shape / rate:.4g}")
    return 0

# Run only when the file is executed directly, not when it's imported
if __name__ == "__main__":
    raise SystemExit(main())
