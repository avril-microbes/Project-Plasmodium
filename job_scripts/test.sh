#!/bin/bash
#SBATCH --mail-user=avril.wang@mail.utoronto.ca
#SBATCH --account=def-mideon  
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10     # add this line to make sure that slurm uses multiple node4              # number of processes
#SBATCH --time=00:30:00         # time (HH:MM:SS)

module load gcc/9.3.0 r/4.0.2

# Export the nodes names. 
# If all processes are allocated on the same node, NODESLIST contains : node1 node1 node1 node1
# Cut the domain name and keep only the node name
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))
R -f test.R
