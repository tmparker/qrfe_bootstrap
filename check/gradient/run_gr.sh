#!/bin/bash
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00-06:00 # time (DD-HH:MM)
#SBATCH --array=0-99

module load r
Rscript sim_gr.R $SLURM_ARRAY_TASK_ID
