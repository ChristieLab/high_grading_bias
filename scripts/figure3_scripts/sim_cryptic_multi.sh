#!/bin/bash
#SBATCH --job-name=multi_cryptic
#SBATCH --array=1-108
#SBATCH -A standby
#SBATCH --mem=200G
#SBATCH -t 4:00:00
#SBATCH --output slurm_ms_%A_%a.out

set -e

module load r/4.2

parmfile=cryptic_parms.txt
outfile=cryptic_multi

Rscript sim_cryptic_multi.R $parmfile $SLURM_ARRAY_TASK_ID ${outfile}_${SLURM_ARRAY_TASK_ID}

sbatch boot.sh ${outfile}_${SLURM_ARRAY_TASK_ID}.RDS ${outfile}_${SLURM_ARRAY_TASK_ID}
