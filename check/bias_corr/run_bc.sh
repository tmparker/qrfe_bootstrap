#!/bin/bash
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00-08:00 # time (DD-HH:MM)
#SBATCH --array=100-199,1000-1099

module load r
Rscript sim_bc.R
