#!/bin/bash
#PBS -l procs=8,tpn=8,mem=46gb,walltime=30:00
#PBS -q cpsc424
#PBS -j oe

cd $PBS_O_WORKDIR

module load Langs/Intel/14

./serial < data/actualdata3 > ./results/actualdata3_c.out
./serial < data/actualdata4 > ./results/actualdata4_c.out
