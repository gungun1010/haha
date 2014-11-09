#!/bin/bash
#PBS -l procs=8,tpn=8,mem=46gb,walltime=30:00
#PBS -q cpsc424
#PBS -j oe

cd $PBS_O_WORKDIR
module load Langs/Intel/14
./task2 /home/fas/hpcprog/ahs3/cpsc424/assignment4/inputs/Square-n/Square-n.19.0.gr /home/fas/hpcprog/ahs3/cpsc424/assignment4/inputs/Square-n/Square-n.19.0.ss > results/square_19_4threads.out
