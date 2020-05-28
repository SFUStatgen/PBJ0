#!/bin/bash
#SBATCH --ntasks=100     # number of processes; instead of few w/ long run time try requesting lots w/ short run time which the scheduler shd favors; note that Cedar has 94000 processors so asking for 100 is a drop in the bucket.
#SBATCH --mem-per-cpu=4000M      # memory; default unit is megabytes
#SBATCH --time=02:00:00           # time (HH:MM:SS) - run time per processor; changed from 2days to 2h!


module load gcc/5.4.0 # needed for SNPknock
module load r/3.5.0

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_LIBS=/home/jgraham/R/x86_64-pc-linux-gnu-library/3.5 #This will need to be changed to your compute canada account info
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))

R CMD BATCH --no-save dCorKnockoff.R
