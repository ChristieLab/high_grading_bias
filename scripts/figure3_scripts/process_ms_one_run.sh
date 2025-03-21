#!/bin/bash
#SBATCH --job-name=boot_PCA
#SBATCH -A standby
#SBATCH --mem=120G
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 4:00:00
#SBATCH --array=1-100
#SBATCH --output slurm_files/slurm_boot_pca_%A_%a.out
module purge
module load r/4.2

data=$1
iter=$2
out=../results/PCA_boots/no_struct_multi/permutation

Rscript boot.R $data ${out}_${iter} ${SLURM_ARRAY_TASK_ID}
