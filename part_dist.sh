#!/bin/bash
#SBATCH --ntasks=4               # number of processes
#SBATCH --mem-per-cpu=4000M      # memory; default unit is megabytes
#SBATCH --time=0-1:30            # time (DD-HH:MM)


module load r/3.5.0

R CMD BATCH --no-save 2_part_dist.R