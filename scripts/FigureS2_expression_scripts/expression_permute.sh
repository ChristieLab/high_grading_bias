#!/bin/bash
#SBATCH --job-name=permute_expression
#SBATCH --array=1-100
#SBATCH -A standby 
#SBATCH --mem=20G
#SBATCH -N 1
#SBATCH -n 9
#SBATCH -t 4:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=lee3617@purdue.edu
#SBATCH --output ./slurm_files/slurm_permute_%A_%a.out
module purge
module load r/4.3.3

Rscript expression_permute.R $SLURM_ARRAY_TASK_ID
