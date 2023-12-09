# yeast_ts
code for yeast time series paper

Obtain raw fastq files for: 
- founders	@ 
- BAS02 @ 
- evolved populations @    (this paper) 
	
align raw reads to reference genome (S288c)

```bash
bwa mem -t 1 -M $ref ${R1name} ${R2name} | samtools view -bS - > $name.temp.bam
samtools sort $name.temp.bam -o $name.bam
samtools index $name.bam
```
	
possibly merge samples sequenced across two or more runs

```bash
bamtools merge -in $name.rep1.bam -in $name.rep2.bam -out $name.merged.bam
```
	
create a table of SNP frequencies from bams by chromosome 

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

cleanup SNP table to include only well behaved SNPs and other QC

```bash
# names.txt are the headers for the frequency table
# "chr\tpos\tsample_1\t ... sample_n\n"
# samples are in order of "bam.list.txt"
cat names.txt > SNP.freq.txt
cat fSNP.chr*.txt >> SNP.freq.txt
```

do a little cleanup, and saves more R friendly 
in R

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

{include SNPtable on website}

call haplotypes 

requires: 
- list of samples to impute haplotypes (=callhaps.txt), with as many lines as size of array job 
- an output folder (=OUTPUT) 
- the frequency table above 
- list of founders (=founder.file.Sept5_21.txt) 
- window size (=50) 

- haplotyper.justhapfreq.slurm.sh 
- haplotyper.justhapfreq.R 
- haplotyper.limSolve.simple.code.R 
 
- install limSolve in R (i.e., install.packages("limSolve")) 
	
```bash
sbatch --array=1-146 haplotyper.justhapfreq.slurm.sh callhaps.txt OUTPUT SNP.freq.txt founder.file.Sept5_21.txt 50

# merge data over samples
cat OUTPUT/sample1_hap_freq.txt | head -n 1 > allhaps.txt
awk FNR-1 OUTPUT/*_hap_freq.txt >> allhaps.txt
cat allhaps.txt | gzip -c > Sept5_22.allhaps.Ns50.txt.gz
```

{include haplotpe calls on website}

impute SNPs from haplotypes 

requires 
- impute.SNP.Sept22.sh 
- impute.SNP.Sept22.R   # paths may be hardwired... 
 
- haplotype calls (=Sept5_22.allhaps.Ns50.txt.gz) 
- samples to call (=callhaps.txt) 
- SNP table (=SNP.freq.txt) 

```bash
sbatch scripts/impute.SNP.Sept22.sh

# results in -> impute_SNP_Sept22/
# CHROM	POS	POOL	iGeno	oGeno
# chrIV 	 18422 	 SEE01B02CD600R04 	 0.9994 	 1 
# iGeno = imputed allele frequency frequency
# oGeno = observed {raw freq} allele frequency

cat impute_OUT/$sample_1.impute.SNP.txt | head -n 1 > allSNPs.txt
awk FNR-1 impute_OUT/*.impute.SNP.txt >> allSNPs.txt
cat allSNPs.txt | gzip -c > allSNPs.Ns50impute.txt.gz
```

Then some tweaks fairly specific to these particular samples in R 

```R
xx = read.table("allSNPs.Ns50impute.txt.gz", header=TRUE)
library(tidyverse)
xx2 = xx %>%
	filter(! POOL %in% c("A1xB3_diploid_controls","BAS02","BAS02r")) %>%
	mutate(drug = substr(as.character(POOL), 9, 10)) %>%
	mutate(rep = substr(as.character(POOL), 15, 16)) %>%
	mutate(week = substr(as.character(POOL), 4, 5)) %>%
	mutate(base = "N")
	
xx3 = xx %>% filter(POOL == "BAS02") %>%	
	mutate(drug = "XX") %>%
	mutate(rep = "00") %>%
	mutate(week = "00") %>%
	mutate(base = "T")

xx4 = xx %>% filter(POOL == "BAS02r") %>%	
	mutate(drug = "XX") %>%
	mutate(rep = "00") %>%
	mutate(week = "00") %>%
	mutate(base = "R")

###  I am trying to encode the base base populations here
###  they are common to all treatments
###  there are two of them BAS02 = "T" and BAS02r = "R" 
###  or "N" = not a base population
	
write.table(rbind(xx2,xx3,xx4),"allSNPs2.Ns50impute.txt")
```

{include imputed SNPs on website}

script to make figures
