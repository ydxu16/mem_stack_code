#!/bin/sh
#SBATCH -N 1                           # nodes=4
#SBATCH --ntasks-per-node=1            # ppn=1
#SBATCH --cpus-per-task=4
#SBATCH -t 24:00:00                     # 13 hours walltime
#SBATCH -o out -e err

#SBATCH -J job_name          # job name
#SBATCH --mail-user=yuanda@princeton.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load openmpi
module load fftw/intel-14.0/3.3.3
./solver 0.001 1200000 20
