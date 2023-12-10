args = commandArgs(trailingOnly=TRUE)
pool = args[1]
myfilename = paste0("impute_SNP_Sept22/",pool,".impute.SNP.txt",sep='')
cat("CHROM\tPOS\tPOOL\tiGeno\toGeno\n",file=myfilename)
SNP = read.table("results/bcf2/SNP.freq.txt")
found = SNP[,c(1:14,16:20)]
het = apply(found[,-c(1:2)],1,function(x) sum(x*(1-x)*2))
found = found[het<0.001,]
obsGeno = cbind(SNP[,c(1,2)],SNP[,pool])
colnames(obsGeno)[3] = pool
xx = read.table("Sept5_22.allhaps.Ns50.txt.gz",header=TRUE)
xxp = xx[xx$pool == pool,]
rm(xx)
rm(SNP)
for(chrs in c("chrIV","chrXV","chrVII","chrXII","chrXVI","chrXIII","chrII","chrXIV","chrX","chrXI","chrV","chrVIII","chrIX","chrIII","chrVI","chrI")){
	xxc = xxp[xxp$chr==chrs,]
	foundc = found[found$CHROM==chrs,]
	obsGenoc = obsGeno[obsGeno$CHROM==chrs,]
	for(i in 1:nrow(foundc)){
		mpos = foundc$POS[i]
		weights = as.numeric(foundc[i,3:19])	
		closest = xxc$pos[which.min(abs(xxc$pos - mpos))]
		hapf = as.numeric(unlist(strsplit(xxc$founderfreqs[xxc$pos==closest],";")))
		cat(chrs,"\t",mpos,"\t",pool,"\t",sum(hapf * weights),"\t",obsGenoc[obsGenoc$POS==mpos,pool],"\n",file=myfilename,append=TRUE) 
		} # SNPs in chr
	} #chr


