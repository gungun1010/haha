#!/bin/bash
#PBS -l procs=2,tpn=2,mem=46gb
#PBS -l walltime=15:00
#PBS -N result
#PBS -r n
#PBS -j oe
#PBS -q cpsc424

module load Langs/Intel/14 MPI/OpenMPI/1.6.5
pwd
cd $PBS_O_WORKDIR
pwd
cat $PBS_NODEFILE
time mpiexec -n 2 nonBlockingMPI
