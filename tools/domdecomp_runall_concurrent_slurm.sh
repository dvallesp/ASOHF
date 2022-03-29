#!/bin/bash

#SBATCH --job-name=asohf
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --array=0-7

#SBATCH --output=asohf_%A_%a.out
#SBATCH --error=omp_%j.err

ulimit -c 10
ulimit -s unlimited
ulimit -d unlimited
export OMP_STACKSIZE=4000m
export KMP_BLOCKTIME=100000

declare -i i
declare -i j
declare -i k
declare -i ID
ID=$SLURM_ARRAY_TASK_ID
i=$ID/4 
j=$ID-4*$i
j=$j/2
k=$ID-4*$i-2*$j

cd domain_${i}_${j}_${k}
echo "domain_${i}_${j}_${k}"

srun ./asohf.x

## SEND TO THE QUEUE WITH A COMMAND LIKE:
## (this is an example for supercomputer LluisVives@UV
## sbatch --qos=thin_astro --mem 180g domdecomp_runall_concurrent_slurm.sh 

