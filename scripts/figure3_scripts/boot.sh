#!/bin/bash
#SBATCH --job-name=boot_PCA
#SBATCH --array=1-100
#SBATCH -A standby
#SBATCH --mem=140G
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 4:00:00
#SBATCH --output ./slurm_files/slurm_boot_pca_%A_%a.out
module purge
module load r/4.2

data=$1
out=$2
Rscript boot.R $data $out $SLURM_ARRAY_TASK_ID
