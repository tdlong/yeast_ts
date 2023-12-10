blah = function(par, tw, NSNPs){
	ww = exp(-unlist(tw)/(2*par^2))
	ww = ww/sum(ww)
	# the idea is to choose sigma, such that the NSNPs closest to pos
	# account for 50% of the weight
	(sum(sort(ww,decreasing=TRUE)[(NSNPs+1):length(ww)])-0.50)^2
	}
    
runscan = function(poolroot, folder, SNPtable, foundernames, NSNPs){

	#  poolroot = pool to get frequencies for
	#  folder = output folder for results
	#  SNPtable = file of SNP frequencies in pool and founders.  The structure of the file matters
	#     the 1st column is a chromosome name (CHROM), and the 2nd a position (POS), followed by SNP frequencies
	#  foundernames = file with the names of the founder haplotypes (and a short name)
	#  NSNPs = number of SNPs closest to position that account for 50% of Gaussian weights
	xx = read.table(SNPtable,header=TRUE)
	tempchrs=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrM")
	chromlengths=c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066,85779)
	KONSTANT = data.frame(CHROM=tempchrs,nchr=1:17,offset=c(0,cumsum(chromlengths[1:16])))
	yy = read.table(file=foundernames,header=FALSE)
	founderroots = yy[,1]
	founders = match(founderroots,names(xx))
	pool = match(poolroot,names(xx))
	Nf = length(founders)
	foundershortnames = yy[,2]
	filename=paste(folder,"/",poolroot,"_hap_freq.txt",sep='')
	
	if (!file.exists(filename)){	
	
	stepSize = 1000
	WINSIZE = 50000
	mynames = paste(foundershortnames, collapse=";")
	
	lppp = 0
	for(ccc in 1:17){
		Maxpos = max(xx$POS[xx$CHROM == tempchrs[ccc]],na.rm=TRUE)
		Minpos = min(xx$POS[xx$CHROM == tempchrs[ccc]],na.rm=TRUE)
		ppp = seq(round((Minpos+5000)/100,0)*100,round((Maxpos-5000)/100,0)*100,stepSize)
		lppp = lppp + length(ppp)
		}
	ddd = data.frame(poolroot=rep("",lppp),chr=rep("",lppp),nchr=rep(0,lppp),pos=rep(0,lppp),mpos=rep(0,lppp),NSNPs=rep(0,lppp),foundernames=rep("",lppp), cutree=rep("",lppp),founderfreqs=rep("",lppp),stringsAsFactors=FALSE)

	i=1
	for(ccc in 1:17){
		Maxpos = max(xx$POS[xx$CHROM == tempchrs[ccc]],na.rm=TRUE)
		Minpos = min(xx$POS[xx$CHROM == tempchrs[ccc]],na.rm=TRUE)
		ppp = seq(round((Minpos+5000)/100,0)*100,round((Maxpos-5000)/100,0)*100,stepSize)
		for(pos in ppp){
			#pos=9700
			#pos=50700
			# median non-recombined block size is 14kb
			#  so step size should be much smaller than that (500bp = 1/30th)
			#  I would guess optimum window size is +/- 7kb (= 14 total) ish.  
			#  but markers are limiting...
			predictors = subset(xx, CHROM == tempchrs[ccc] & POS > pos - WINSIZE & POS < pos + WINSIZE, select=founders)
			Y = subset(xx, CHROM == tempchrs[ccc] & POS > pos - WINSIZE & POS < pos + WINSIZE, select=pool)
			tw = (subset(xx, CHROM == tempchrs[ccc] & POS > pos - WINSIZE & POS < pos + WINSIZE, select=POS)-pos)^2
			out <- optimize(blah,c(2000,50000), tw, NSNPs)
			sigma=out$minimum
			weights = exp(-tw/(2*sigma^2))

			# clearly this is very much a heuristic
			# the current code basically will not calculate haplotype frequencies if two founders are really
			# similar in some window.  Another strategy would be to calculate haplotype frequencies, but print 
			# minManhattan.  Even better would be to think about which haplotypes are not distinguishable, and print that
		
			# Note that the haplotype frequencies are constrained to sum to one and have a minimum 
			# frequency of 0.002.  This means the maximum haplotpe frequency is 96.8%.
			# There are 17 haplotypes, if 16 have a frequency of 0.2% then the 17th
			# is 96.8%

			MynumberSNP = nrow(predictors)
			Freqs = rep(NA,Nf)
			Groups = rep(NA,Nf)
			if(nrow(predictors)>=10){
				Groups = cutree(hclust(as.dist(1-cor(predictors*unlist(weights))^2)),h=0.005)
				nGroups = nlevels(as.factor(Groups))
				
				###  this is the only place that the predictors change over samples...
				###  I moved it to after I decide on the groups
				###  the new danger ... is that once I drop rows of predictors haplotypes cannot be distinguished anymore
				
				predictNotMissing = apply(predictors,1,function(x) sum(is.na(x))==0) & (!is.na(Y))
				predictors = predictors[predictNotMissing,]
				Y = Y[predictNotMissing,]
				tw = tw[predictNotMissing,]
				weights = weights[predictNotMissing]
				
				Groups2 = cutree(hclust(as.dist(1-cor(predictors*unlist(weights))^2)),h=0.005)
				nGroups2 = nlevels(as.factor(Groups2))
				if(paste(as.numeric(Groups),collapse=";") == paste(as.numeric(Groups2),collapse=";")){
					d = ncol(predictors)
					A = predictors
					B = Y
					E = t(matrix(rep(1,d)))
					F = 1
					G = diag(rep(1,d))
					H = matrix(rep(0.0003,d))
					Wa = weights
					out = lsei(A=A,B=B,E=E,F=F,G=G,H=H,Wa=Wa,verbose=TRUE)
					Freqs = out$X
					} # groups not changed
				} # nrow > 10 
				
			mpos = pos + KONSTANT$offset[KONSTANT$CHROM==tempchrs[ccc]]
			nchr = KONSTANT$nchr[KONSTANT$CHROM==tempchrs[ccc]]
			groups = paste(as.numeric(Groups),collapse=";")
			freqs = paste(round(as.numeric(Freqs),4),collapse=";")
			ddd[i,] = c(poolroot, tempchrs[ccc], nchr, pos, mpos, MynumberSNP, mynames, groups, freqs)
			i=i+1
			}  # over positions
		} # over chromosome
	write.table(ddd,filename,row.names=FALSE,sep="\t",quote=FALSE)
	} # exists 
} # function specific for pool


