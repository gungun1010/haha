#!/bin/bash
#PBS -l procs=8,tpn=4,mem=46gb,walltime=30:00
#PBS -q cpsc424
#PBS -j oe

cd $PBS_O_WORKDIR

module load Langs/Intel/14

time mpiexec -n 8 task2 < ../data/actualdata1 > ./results/actualdata1_c.out
# To run other data set (2 - 4) , replace all occurance of "actualdata1" with "actualdataN", where N is the dataset number 2-4
