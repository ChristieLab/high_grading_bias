#!/bin/bash
#SBATCH --job-name=scrm
#SBATCH -A standby
#SBATCH --mem=120G
#SBATCH -t 4:00:00
#SBATCH -A christ99
module purge
module load r/4.2

Rscript sim_eight_pops_scrm.R
