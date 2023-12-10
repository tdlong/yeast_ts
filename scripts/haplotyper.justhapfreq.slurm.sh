#!/bin/bash
#SBATCH --job-name=haplotyper
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1

module load R/4.1.2
file=$1
folder=$2
SNPtable=$3
foundernames=$4
NSNPs=$5

pool=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f 1` 

echo "Rscript scripts/haplotyper.justhapfreq.R $pool $folder $SNPtable $foundernames $NSNPs"
Rscript --verbose scripts/haplotyper.justhapfreq.R $pool $folder $SNPtable $foundernames $NSNPs


