#!/bin/bash
#SBATCH --account=def-mideon  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48     # add this line to make sure that slurm uses multiple node4              # number of processes
#SBATCH --mem=0      # memory; default unit is megabytes
#SBATCH --time=01:00:00         # time (HH:MM:SS)

module load gcc/9.3.0
module load r/4.0.2
module load openmpi/4.0.3
export R_LIBS=~/local/R_libs/

mpirun -np 1 R CMD BATCH 2022-02-01_pso-si-I-log.R test.txt