#!/bin/bash
#SBATCH --job-name=distributions
#SBATCH --account=pi-kravtsov
#SBATCH --error=distr.err
#SBATCH --out=distr.out
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python distributions.py 6 3 > halo_6_distr.out
