#!/bin/bash
#SBATCH --ntasks=100             # number of processes
#SBATCH --mem-per-cpu=4000M      # memory; default unit is megabytes
#SBATCH --time=0-00:59           # time (DD-HH:MM)


module load r/3.5.0

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_LIBS=/home/pnickchi/R/x86_64-pc-linux-gnu-library/3.5
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))

R CMD BATCH --no-save 4_IBD_dCorN.R