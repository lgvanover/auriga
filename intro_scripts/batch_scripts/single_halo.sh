#!/bin/bash -l
#PBS -l select=1:system=polaris
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -q preemptable 
#PBS -A Auriga
#PBS -l filesystems=home:eagle

#cd ${PBS_O_WORKDIR}

# MPI example w/ 8 MPI ranks per node spread evenly across cores
NNODES=`wc -l < $PBS_NODEFILE`

CPU_THREAD=0
HALO=27
echo $HALO
mpiexec -n 1 --ppn 1 --cpu-bind list:$CPU_THREAD,$((CPU_THREAD+32)) python ~/arepo-snap-util/load_new_eps.py $HALO 3 > ~/arepo-snap-util/batch_scripts/outputs/halo$HALO.out
