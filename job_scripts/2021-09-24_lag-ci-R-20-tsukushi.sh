#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=16
#SBATCH --job-name="test"
 module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2
 module load parallel
 Rscript 2021-09-24_lag-ci-R-20-tsukushi.R

