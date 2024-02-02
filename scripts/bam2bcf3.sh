#!/bin/bash
#SBATCH --job-name=callSNPs
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --array=1-17

module load bwa/0.7.17
module load samtools/1.10
module load bcftools/1.10.2

ref="/dfs7/adl/tdlong/RL_yeast/ref/S288c.fasta"
dir1="results/bcf3"

declare -a chrs=("chrI" "chrII" "chrIII" "chrIV" "chrIX" "chrM" "chrV" "chrVI" "chrVII" "chrVIII" "chrX" "chrXI" "chrXII" "chrXIII" "chrXIV" "chrXV" "chrXVI")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

bcftools mpileup -I -d 4000 -t $mychr -a "FORMAT/AD,FORMAT/DP" -f $ref -b helperfiles/bam.list.Jan3_23.txt | bcftools call -mv -Ob > $dir1/calls.$mychr.bcf  

