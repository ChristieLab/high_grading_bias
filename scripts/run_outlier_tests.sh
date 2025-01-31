#!/bin/bash
#SBATCH --job-name=outliers
#SBATCH -A christ99
#SBATCH --mem=40G
#SBATCH -t 4:00:00
#SBATCH --output=../slurm_files/slurm_%x_%A_%a.out
#SBATCH --error=../slurm_files/slurm_%x_%A_%a.err

set -e

module load r/4.2

infile=../data/island_sim.RDS
outfile=../results/outliers_island

Rscript run_outlier_tests.R $infile $outfile
