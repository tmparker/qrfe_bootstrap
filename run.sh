#!/bin/bash
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00-24:00 # time (DD-HH:MM)
#SBATCH --array=100-1099

module load r
Rscript sim.R $SLURM_ARRAY_TASK_ID
