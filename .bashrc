module load conda
module load cray-hdf5
#FILESYSTEM=/eagle/projects/Auriga/bin/lib
#CONDA=/home/lgvanover:$CONDA
#OTHER=/home/lgvanover/arepo-snap-util/libs:$OTHER
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/eagle/projects/Auriga/bin/lib:/home/lgvanover
#PATH=$PATH:$FILESYSTEM:$CONDA:$OTHER
export LD_LIBRARY_PATH
conda activate arepo_util
#export PATH

