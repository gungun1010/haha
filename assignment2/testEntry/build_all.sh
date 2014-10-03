#!/bin/bash
#PBS -q cpsc424
#PBS -l procs=8,tpn=8,mem=46gb
#PBS -l walltime=15:00
#PBS -N result
#PBS -r n
#PBS -j oe

module load Langs/Intel/14 MPI/OpenMPI/1.6.5
printf "****************************************************************************\n"
printf "*                       making blockingMPI                                 *\n"
printf "****************************************************************************\n"

pwd
cd $PBS_O_WORKDIR
pwd
#cat $PBS_NODEFILE
cd ../BlockingMPI
make clean
make
printf "\n"
printf "****************************************************************************\n"
printf "*                        making Non blockingMPI                            *\n"
printf "****************************************************************************\n"
pwd
cd $PBS_O_WORKDIR
pwd
#cat $PBS_NODEFILE
cd ../NonBlockingMPI
make clean
make
printf "\n"
printf "****************************************************************************\n"
printf "*                        making balance load  non blockingMPI              *\n"
printf "****************************************************************************\n"
pwd
cd $PBS_O_WORKDIR
pwd
#cat $PBS_NODEFILE
cd ../BalanceLoad
make clean
make

cd ../testEntry

printf "\n"
pwd
cd $PBS_O_WORKDIR
pwd
cat $PBS_NODEFILE
printf "\n"
printf "****************************************************************************\n"
printf "*                        running  blockingMPI                              *\n"
printf "****************************************************************************\n"
time mpiexec -n 8 ../BlockingMPI/blockingMPI
printf "\n"
printf "****************************************************************************\n"
printf "*                        running  non blockingMPI                          *\n"
printf "****************************************************************************\n"
time mpiexec -n 8 ../NonBlockingMPI/nonBlockingMPI
printf "\n"
printf "****************************************************************************\n"
printf "*                        running  load balance non blocking                *\n"
printf "****************************************************************************\n"
time mpiexec -n 8 ../BalanceLoad/nonBlockingMPI
