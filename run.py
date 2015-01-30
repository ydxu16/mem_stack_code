from subprocess import call
from numpy import arange
from sys import argv
dt = 0.001
Nstep = 400000
num_layers = 5
coupling = 0.001 
Nx = 128
Ny = 128
eta_M = 1
eta_S = 1
N_print_gap = 80
ave_psi = 0.3
st=  "#!/bin/bash \n  #parallel job using 1 nodes and 4 CPU cores and runs for 1 hours(Max). \n #PBS -l nodes=1:ppn=4,walltime=24:00:00 \n #PBS -o out -e err \n  #sends mail if the process abors, when it begins, and when it ends(abe) \n  #PBS -m abe \n module load intel/11.1/64/11.1.075 \n cd $PBS_O_WORKDIR \n  #numprocs=`echo ${PBS_NODEFILE} | wc -l` \n ./solver" + str(dt) + str(Nstep) + str(num_layers) + str(coupling) + str(Nx) + str(Ny) + str(eta_M) + str(eta_S) + str(N_print_gap) + str(ave_psi)
 
