#!/bin/bash
#SBATCH --job-name=multi_cryptic
#SBATCH --array=103,104,105,106,107,114,115,117,122,123,124,127,130,131,134,135,137,139,143,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200%10
#SBATCH -A christ99
#SBATCH --mem=200G
#SBATCH -t 6-23:59:59
#SBATCH --output slurm_ms_%A_%a.out

# set -e

module load r/4.2

parmfile=candidate_parms.txt
outfile=candidate_cryptic_multi
outdir=cryptic_multi_permutation

Rscript sim_cryptic_multi_w_control.R $parmfile $SLURM_ARRAY_TASK_ID ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}

############################################################


if [ -f ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}.RDS ]; then
    echo "submitting boots for test"
    sbatch --wait boot.sh ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}.RDS ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}
    wait
fi


############################################################


if [ -f ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}_control.RDS ]; then
    echo "submitting boots for control"
    sbatch --wait boot.sh ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}_control.RDS ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}_control
    wait
fi

echo "finished"
