#!/bin/bash
#PBS -q cpsc424

module load Langs/Intel/14 MPI/OpenMPI/1.6.5
printf "****************************************************************************\n"
printf "*                       making blockingMPI                                 *\n"
printf "****************************************************************************\n"
cd ../BlockingMPI
make clean
make

printf "****************************************************************************\n"
printf "*                        making Non blockingMPI                            *\n"
printf "****************************************************************************\n"
cd ../NonBlockingMPI
make clean
make

printf "****************************************************************************\n"
printf "*                        making balance load  non blockingMPI              *\n"
printf "****************************************************************************\n"
cd ../BalanceLoad
make clean
make

cd ../testEntry

#PBS -l procs=8,tpn=8,mem=46gb
#PBS -l walltime=15:00
#PBS -N result
#PBS -r n
#PBS -j oe

pwd
cd $PBS_O_WORKDIR
pwd
cat $PBS_NODEFILE
time mpiexec -n 8 ../BlockingMPI/blockingMPI
