#!/bin/bash
#SBATCH --job-name=rawPoll
#SBATCH -A tdlong_lab
#SBATCH -p highmem          
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00  

module load R/4.2.2
Rscript scripts/impute_SNPs.R
Rscript scripts/impute2_SNPs.R

