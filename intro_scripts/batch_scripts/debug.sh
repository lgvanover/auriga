#!/bin/bash -l
#PBS -l select=1:system=polaris
#PBS -l place=scatter
#PBS -l walltime=10:00:00
#PBS -q preemptable
#PBS -A Auriga
#PBS -l filesystems=home:eagle

#cd ${PBS_O_WORKDIR}

# MPI example w/ 8 MPI ranks per node spread evenly across cores
NNODES=`wc -l < $PBS_NODEFILE`

CPU_THREAD=0
MAX_HALO=30
echo $HALO
for (( i=1; i <= $MAX_HALO; ++i))
do
    mpiexec -n 1 --ppn 1 --cpu-bind list:0,1 hostname
   
    CPU_THREAD=$((CPU_THREAD+1))
done

wait
