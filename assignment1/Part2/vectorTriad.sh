#!/bin/bash
#PBS -l procs=1,tpn=1,mem=24gb,walltime=60:00 
#PBS -q cpsc424
cd $PBS_O_WORKDIR
./vectorTriad_o0
