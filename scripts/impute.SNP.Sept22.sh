#!/bin/bash
#SBATCH --job-name=impute
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=1
#SBATCH --array=1-146
 
module load R/4.1.2

file="helperfiles/callhaps.txt"

pool=`head -n $SLURM_ARRAY_TASK_ID $file | cut -f 1 | tail -n 1`
Rscript scripts/impute.SNP.Sept22.R $pool

