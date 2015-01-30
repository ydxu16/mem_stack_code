#!/bin/bash 
#parallel job using 1 nodes and 4 CPU cores and runs for 1 hours(Max).
#PBS -l nodes=1:ppn=4,walltime=24:00:00
#PBS -o out -e err 
#sends mail if the process abors, when it begins, and when it ends(abe) 
#PBS -m abe 
module load intel/11.1/64/11.1.075 
cd $PBS_O_WORKDIR 
#numprocs=`echo ${PBS_NODEFILE} | wc -l` 
./solver 0.001 400000 5 0.001 128 128 1 1 80 0.3