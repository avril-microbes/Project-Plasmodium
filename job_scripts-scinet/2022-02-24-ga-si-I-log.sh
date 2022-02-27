#!/bin/bash
#SBATCH --account=def-mideon  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48     # add this line to make sure that slurm uses multiple node4              # number of processes
#SBATCH --mem=0      # memory; default unit is megabytes
#SBATCH --time=06:00:00         # time (HH:MM:SS)

cd $SLURM_SUBMIT_DIR
module load gcc/8.3.0 intel/2019u4 r/4.1.2

# load R script to all nodes
${SCINET_R_ROOT}/job_scripts --no-restore 2022-02-24-ga-si-I-log.R

