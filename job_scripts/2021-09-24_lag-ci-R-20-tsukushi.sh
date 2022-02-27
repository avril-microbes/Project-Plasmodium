#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=8
#SBATCH --job-name="test"
 module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
 Rscript 2021-09-24_lag-ci-R-20-tsukushi.R

