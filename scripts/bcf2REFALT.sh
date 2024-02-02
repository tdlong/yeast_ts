#!/bin/bash
#SBATCH --job-name=callSNPs
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --array=1-17

module load bwa/0.7.17
module load samtools/1.10
module load bcftools/1.10.2
module load R/4.2.2

ref="/dfs7/adl/tdlong/RL_yeast/ref/S288c.fasta"
dir1="results/bcf3"
dir2="results/process"

declare -a chrs=("chrI" "chrII" "chrIII" "chrIV" "chrIX" "chrM" "chrV" "chrVI" "chrVII" "chrVIII" "chrX" "chrXI" "chrXII" "chrXIII" "chrXIV" "chrXV" "chrXVI")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

echo -ne "CHROM\tPOS" > $dir2/RefAlt.$mychr.txt
bcftools query -l $dir1/calls.$mychr.bcf | awk '{printf("\tREF_%s\tALT_%s",$1,$1)}' >> $dir2/RefAlt.$mychr.txt
echo -ne "\n" >> $dir2/RefAlt.$mychr.txt
bcftools view -m2 -M2 -v snps -i 'QUAL>20' $dir1/calls.$mychr.bcf | bcftools query -e'GT ="./."'  -e'QUAL<60' -f'%CHROM %POS [ %AD{0} %AD{1}] [%GT]\n' | grep -v '\.' | awk 'NF-=1' >>$dir2/RefAlt.$mychr.txt
Rscript scripts/newSNP.R $mychr
