#!/bin/bash
#COBALT -n 1 
#COBALT -t 30 
#COBALT -A Auriga
#COBALT -q debug-cache-quad
#COBALT --attrs filesystems=home,eagle

NNODES=${COBALT_JOBSIZE}
NRANKS_PER_NODE=1
NTHREADS_PER_CORE=1
NDEPTH=64

NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))

# option  long version            (explanation)
#
# -n                              "PEs" (ranks)
# -N      --pes-per-node          ranks per node
# -d      --cpus-per-pe           hyperthreads per rank
# -cc     --cpu-binding depth
# -j                              cpus (hyperthreads) per compute unit (core)


aprun -n 1 -N 1 -d 64 -j 1 python load_new_eps.py 6 3
status=$?


echo "Exit status of aprun is: $status"
exit $status
