#!/bin/bash
#SBATCH --job-name=load_halos
#SBATCH --account=pi-kravtsov
#SBATCH --error=load_halos.err
#SBATCH --out=load_halos.out
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python halo_testing.py 24 3 new > halo_24.out  
 

