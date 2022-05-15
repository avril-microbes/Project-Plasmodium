#!/bin/bash
#SBATCH --account=def-mideon   # replace this with your supervisors account
#SBATCH --nodes=1                # number of node MUST be 1
#SBATCH --ntasks=1        # number of processes      # memory; default unit is megabytes
#SBATCH --cpus-per-task=32
#SBATCH --time=00:45:00        # time (DD-HH:MM)
#SBATCH --mail-user=avril.wang@mail.utoronto.com # Send email updates to you or someone else
#SBATCH --mail-type=ALL          # send an email in all cases (job started, job ended, job aborted)

module load gcc/9.3.0 r/4.0.2
R -f 2022-05-12_I1-I2-cf.R