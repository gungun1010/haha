#!/bin/bash
#PBS -l procs=8,tpn=4,mem=46gb,walltime=30:00
#PBS -q cpsc424
#PBS -j oe

cd $PBS_O_WORKDIR

module load Langs/Intel/14

time mpiexec -n 8 parallel < ../data/actualdata4 > ./results/actualdata4_c.out
