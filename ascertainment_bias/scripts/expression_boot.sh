#!/bin/bash
#SBATCH --job-name=boot_expression
#SBATCH --array=1-100
#SBATCH -A standby 
#SBATCH --mem=20G
#SBATCH -N 1
#SBATCH -n 9
#SBATCH -t 4:00:00
#SBATCH --output ./slurm_files/slurm_boot_pca_%A_%a.out
module purge
module load r/4.3.3

Rscript expression_boot.R $SLURM_ARRAY_TAS