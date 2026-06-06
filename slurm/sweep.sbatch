#!/bin/bash
#SBATCH --job-name=""
#SBATCH --partition=small
#SBATCH --account=cademartorigrp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1GB
#SBATCH --mail-user=hatan25@wfu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=00-01:00:00
#SBATCH --output=hello-%j.out
#SBATCH --array1-5

echo "Hello, World!"
hostname