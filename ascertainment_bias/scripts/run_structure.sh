#!/bin/bash -l
#SBATCH -n 4
#SBATCH --array=1-10
#SBATCH --mem=20G
#SBATCH -t 12:00:00
#SBATCH -J structure

module load R/4.2

dat=monarchs_random.RDS
target_dir=all_sites_structure_res

Rscript run_structure.R struct_parmfile.txt $SLURM_ARRAY_TASK_ID $dat $target_dir
