#!/bin/bash
#PBS -l procs=8,tpn=8,mem=46gb,walltime=30:00
#PBS -q cpsc424
#PBS -j oe

cd $PBS_O_WORKDIR
module load Langs/Intel/14
./task1 ../inputs/Square-n/Square-n.14.0.gr ../inputs/Square-n/Square-n.14.0.ss
