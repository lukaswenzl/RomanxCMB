#PBS -S /bin/bash
#PBS -N fisher
#PBS -m abe
#PBS -M ljw232@cornell.edu
#PBS -l select=1:ncpus=28:mpiprocs=14:model=bro
#PBS -l place=scatter:excl
#PBS -l walltime=00:30:00
#PBS -q devel

# Load modules
source /usr/local/lib/global.profile
module load mpi-hpe/mpt.2.21 comp-intel/2018.3.222 python3/3.7.0

# Lukas
#source /nasa/jupyter/4.4/miniconda/etc/profile.d/conda.sh
#conda activate cosmosis
# Source cosmosis config
source $HOME/cosmosis/config/setup-pleiades-cosmosis
#
# ​Cyrille
# Source cosmosis config
#source $HOME/cosmosis/config/setup-my-cosmosis

export RUN_FOLDER=$PBS_O_WORKDIR
export RUN_NAME="${PBS_O_WORKDIR##*/}"

cd $RUN_FOLDER/..
sh check_versions.sh
NCORES=$(wc -w < $PBS_NODEFILE)

export OMP_NUM_THREADS=2
export DATAFILE=6x2pt_Roman_SO_v1_2_bf26108.fits

#export MPI_LAUNCH_TIMEOUT=40
#/u/scicon/tools/bin/several_tries 
mpiexec -n $NCORES cosmosis --mpi ${RUN_FOLDER}/fisher.ini
## mpiexec -n 14 cosmosis --mpi ${RUN_FOLDER}/fisher.ini
