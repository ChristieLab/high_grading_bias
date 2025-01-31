#!/bin/bash
#SBATCH --job-name=multi_cryptic
#SBATCH --array=1-300%10
#SBATCH -A christ99
#SBATCH --mem=1G
#SBATCH -t 4-4:00:00
#SBATCH --output slurm_ms_%A_%a.out

# set -e

module load r/4.2

parmfile=candidate_parms.txt
outfile=candidate_cryptic_multi

#Rscript sim_cryptic_multi_w_control.R $parmfile $SLURM_ARRAY_TASK_ID ${outfile}_${SLURM_ARRAY_TASK_ID}

############################################################


if [ -f ${outfile}_${SLURM_ARRAY_TASK_ID}.RDS ]; then
    echo "submitting boots for test"
    sbatch --wait boot.sh ${outfile}_${SLURM_ARRAY_TASK_ID}.RDS ${outfile}_${SLURM_ARRAY_TASK_ID}
    wait
fi


############################################################


if [ -f ${outfile}_${SLURM_ARRAY_TASK_ID}_control.RDS ]; then
    echo "submitting boots for control"
    sbatch --wait boot.sh ${outfile}_${SLURM_ARRAY_TASK_ID}_control.RDS ${outfile}_${SLURM_ARRAY_TASK_ID}_control
    wait
fi

echo "finished"
