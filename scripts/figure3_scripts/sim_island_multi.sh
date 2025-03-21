#!/bin/bash
#SBATCH --job-name=multi_island
#SBATCH --array=1-100%9
#SBATCH -A christ99
#SBATCH --mem=200G
#SBATCH -t 6-23:59:59
#SBATCH --output slurm_ms_%A_%a.out

# set -e

module load r/4.2

outfile=island_multi
outdir=island_permutation

Rscript sim_island.R $SLURM_ARRAY_TASK_ID ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}

############################################################
if [ -f ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}.RDS ]; then
    echo "submitting boots for test"
    sbatch --wait boot.sh ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}.RDS ${outdir}/${outfile}_${SLURM_ARRAY_TASK_ID}
    wait
fi

echo "all boots finished"

