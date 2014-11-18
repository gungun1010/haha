#!/bin/bash
#PBS -l procs=1,tpn=1,mem=46gb,walltime=30:00
#PBS -q cpsc424
#PBS -j oe

cd $PBS_O_WORKDIR
module load Langs/Intel/14
#./task2 /home/fas/hpcprog/ahs3/cpsc424/assignment4/inputs/Square-n/Square-n.14.0.gr /home/fas/hpcprog/ahs3/cpsc424/assignment4/inputs/Square-n/Square-n.14.0.ss
./task2 /home/fas/hpcprog/ahs3/cpsc424/assignment4/inputs/Long-n/Long-n.14.0.gr /home/fas/hpcprog/ahs3/cpsc424/assignment4/inputs/Long-n/Long-n.14.0.ss
