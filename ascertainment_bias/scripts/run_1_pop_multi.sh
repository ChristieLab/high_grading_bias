#!/bin/bash
#SBATCH --job-name=multi_ms
#SBATCH --array=2-100%9
#SBATCH -A standby
#SBATCH --mem=40G
#SBATCH -t 4:00:00
#SBATCH --output slurm_ms_%A_%a.out

set -e

module load r/4.2

Rscript scrm_call.R $SLURM_ARRAY_TASK_ID

Rscript convert_and_randpop_ms.R ms_out_${SLURM_ARRAY_TASK_ID}.ms $SLURM_ARRAY_TASK_ID

mv ms_out_${SLURM_ARRAY_TASK_ID}.ms_converted.RDS ms_out_${SLURM_ARRAY_TASK_ID}.RDS

#Rscript process_ms_one_run.R ms_out_${SLURM_ARRAY_TASK_ID}.RDS ${SLURM_ARRAY_TASK_ID}

sbatch --wait process_ms_one_run.sh ms_out_${SLURM_ARRAY_TASK_ID}.RDS ${SLURM_ARRAY_TASK_ID}
wait
