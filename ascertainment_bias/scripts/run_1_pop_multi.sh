#!/bin/bash
#SBATCH --job-name=multi_ms
#SBATCH --array=18,19,20,21,24,25,27,32,33,35,36,37,43,44,45,50,55,56,57,58,59,62,63,64,65,66,67,68,69,70,71,72,74,76,78,79,82,83,84,89,90,96,97,98,99,100
#SBATCH -A standby
#SBATCH --mem=40G
#SBATCH -t 4:00:00
#SBATCH --output slurm_ms_%A_%a.out

set -e

module load r/4.2

#Rscript scrm_call.R $SLURM_ARRAY_TASK_ID

#Rscript convert_and_randpop_ms.R ms_out_${SLURM_ARRAY_TASK_ID}.ms $SLURM_ARRAY_TASK_ID

#mv ms_out_${SLURM_ARRAY_TASK_ID}.ms_converted.RDS ms_out_${SLURM_ARRAY_TASK_ID}.RDS

#Rscript process_ms_one_run.R ms_out_${SLURM_ARRAY_TASK_ID}.RDS ${SLURM_ARRAY_TASK_ID}

sbatch process_ms_one_run.sh ms_out_${SLURM_ARRAY_TASK_ID}.RDS ${SLURM_ARRAY_TASK_ID}

