# yeast_ts
code for yeast time series paper

## Obtain raw fastq files for: 
- founders @ SRA::PRJNA552112
- BAS02 @ SRA::SRX6465400 
- evolved populations @ SRA::   (this paper) 
	
## align raw reads to reference genome (S288c)

```bash
bwa mem -t 1 -M $ref ${R1name} ${R2name} | samtools view -bS - > $name.temp.bam
samtools sort $name.temp.bam -o $name.bam
samtools index $name.bam
```
	
## possibly merge samples sequenced across two or more runs

```bash
bamtools merge -in $name.rep1.bam -in $name.rep2.bam -out $name.merged.bam
```
	
## create a table of SNP frequencies from bams by chromosome 

requires: 
- a file with a list of paths to all the bams to be considered (bam.list.txt) 
- a list of SNPs to consider => $Granny="Granny.in.allFound.txt" 
- a helper = yeast.freqtab.pl 
- a helper = yeast.counttab.pl 
	
```bash
declare -a chrs=("chrI" "chrII" "chrIII" "chrIV" "chrIX" "chrM" "chrV" "chrVI" "chrVII" "chrVIII" "chrX" "chrXI" "chrXII" "chrXIII" "chrXIV" "chrXV" "chrXVI")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

bcftools mpileup -I -d 1000 -t $mychr -a "FORMAT/AD,FORMAT/DP" -f $ref -b bam.list.txt | bcftools call -mv -Ob > calls.$mychr.bcf  
bcftools view -T $Granny calls.$mychr.bcf | bcftools query -e'QUAL<60' -f'%CHROM %POS %REF %ALT [ %AD{0} %AD{1}]\n' | perl yeast.freqtab.pl > fSNP.$mychr.txt
bcftools view -T $Granny calls.$mychr.bcf | bcftools query -e'QUAL<60' -f'%CHROM %POS %REF %ALT [ %AD{0} %AD{1}]\n' | perl yeast.counttab.pl > cSNP.$mychr.txt
```

## cleanup SNP table to include only well behaved SNPs and other QC

```bash
# names.txt are the headers for the frequency table
# "chr\tpos\tsample_1\t ... sample_n\n"
# samples are in order of "bam.list.txt"
cat names.txt > SNP.freq.txt
cat fSNP.chr*.txt >> SNP.freq.txt
```

plus a little more cleanup in R

```R
xx = read.table("SNP.freq.txt",header=TRUE)
yy = read.table("founder.file.Sept5_21.txt",header=FALSE)   # the subset of samples that are founders, 1st column matches names
founderroots = yy[,1]

founders = match(founderroots,names(xx))
Nf = length(founders)
# drop SNPs that are polymorphic in founders
seg = apply(xx[,founders],1,function(x) sum(x>0.05 & x<0.95))
xx = xx[seg==0,]
seg = apply(xx[,founders],1,function(x) sum(unlist(lapply(x,function(y) min((1-y)^2,(0-y)^2)))))
xx = xx[!is.na(seg),]
seg = apply(xx[,founders],1,function(x) sum(unlist(lapply(x,function(y) min((1-y)^2,(0-y)^2)))))
# this throws out about 2%
xx = xx[seg<0.001,]

# manually remove a few SNPs that are problematic to haplotype calling
bchrom = c("chrIII","chrIV","chrIV","chrIV")
bpos = c(159720,617177,1202187,99689)
for (bad in 1:length(bchrom)){
	xx = xx[!(xx$CHROM==bchrom[bad] & xx$POS==bpos[bad]),]
	}	
write.table(xx,"SNP.freq.txt")
```

This SNPtable is available for download (see below)

## call haplotypes 

requires: 
- list of samples to impute haplotypes (=callhaps.txt), with as many lines as size of array job 
- an output folder (=OUTPUT) 
- the frequency table above 
- list of founders (=founder.file.Sept5_21.txt) 
- window size (=50)
<!-- -->
- haplotyper.justhapfreq.slurm.sh 
- haplotyper.justhapfreq.R 
- haplotyper.limSolve.simple.code.R 
<!-- -->
- install limSolve in R (i.e., install.packages("limSolve")) 
	
```bash
sbatch --array=1-146 haplotyper.justhapfreq.slurm.sh callhaps.txt OUTPUT SNP.freq.txt founder.file.Sept5_21.txt 50

# merge data over samples
cat OUTPUT/sample1_hap_freq.txt | head -n 1 > allhaps.txt
awk FNR-1 OUTPUT/*_hap_freq.txt >> allhaps.txt
cat allhaps.txt | gzip -c > Sept5_22.allhaps.Ns50.txt.gz
```

This haplotype table is available for download (see below)

## impute SNPs from haplotypes 

requires 
- impute2_SNPs.R
- impute_SNPs.R
- impute_SNPs.sh
<!-- -->
- haplotype calls (=Sept5_22.allhaps.Ns50.txt.gz) 
- SNP table (=SNP.freq.txt) 

```bash
sbatch impute_SNPs.sh

# output = imputedSNPs.final.Jul23.RDS
# which was converted to a more portable data frame 
```

This *imputed* SNPtable is available for download (see below)

## new mutations

```bash
sbatch scripts/bam2bcf3.sh
sbatch scripts/bcf2REFALT.sh
```

## scripts to make figures

```bash
Rscript Figure2.R
Rscript Figure1_supp.R
Rscript Figure3.R
Rscript Figure4.R
Rscript make_newmutation_figure.R
```

## tarball of 3 tables generated via this pipeline here:

yeast_ts.tar.gz

