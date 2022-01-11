#!/bin/bash
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=250M
#SBATCH --time=00-01:15 # time (DD-HH:MM)
#SBATCH --array=1-1

module load gcc r
Rscript sim.R $SLURM_ARRAY_TASK_ID
