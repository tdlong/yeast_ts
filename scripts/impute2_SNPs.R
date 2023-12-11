library(tidyverse)

flavor_SNP = function(ff, founderfreqs, cutree){

#	ff=BAS$ff[3]
#	founderfreqs=BAS$founderfreqs[3]
#	cutree=BAS$cutree[3]

	blah="AB1;AB2;AB3;AB4;A5;A6;A7;A8;A9;A11;A12;B5;B7;B8;B9;B11;B12"
	templabs = unlist(strsplit(blah, split=';'))
	foStates = as.numeric(unlist(strsplit(ff, split=';')))	
	foFreqs = as.numeric(unlist(strsplit(founderfreqs, split=';')))
	ctv = as.numeric(unlist(strsplit(cutree, split=';')))

	df = data.frame(templabs=templabs,foFreqs=foFreqs,foStates=foStates,ctv=ctv)
	df2 = df %>%
		group_by(ctv) %>%
		summarize(templabs=paste0(templabs, collapse=""),foFreqs=sum(foFreqs),foStates=mean(foStates))	%>%
		ungroup() %>%
		select(-ctv) %>%
		mutate(hap2 = templabs) %>%
		mutate(hap2 = case_when(
			hap2 %in% c("A11","A12","B11","A11A12","A11B11","A12B11","A11A12B11") ~ "A11",
			hap2 %in% c("AB2","A7","B12","AB2A7","AB2B12","A7B12","AB2A7B12") ~ "AB2",
			.default = as.character(hap2)
			)) %>%
		mutate(hap2 = case_when(
			!(hap2 %in% c("AB2","A11","A9","AB4","AB3","B9","B8","B5","A5","A6","AB1","B7","A8","A5B5","AB3B8")) ~ "UNR",
			.default = as.character(hap2)
			)) %>%
		select(-templabs) %>%
		group_by(hap2) %>%
		summarize(foFreqs=sum(foFreqs),foStates=mean(foStates)) %>%
		ungroup() %>%
		filter(foFreqs >0.0005) %>%
		mutate(foFreqs = foFreqs/sum(foFreqs))
	
	Nhaps=nrow(df2)
	MAC = round(min(sum(df2$foStates), Nhaps - sum(df2$foStates)),0)
	MACHAP = NA
	phase = NA
	if(MAC ==1 & round(sum(df2$foStates)) == 1){
		MACHAP = df2$hap2[which.max(df2$foStates)]
		phase = 1
		}else{
		MACHAP = df2$hap2[which.min(df2$foStates)]
		phase = -1
		}
	list(Nhaps=Nhaps,MAC=MAC,MACHAP=MACHAP,phase=phase)
	}	
	
snps_haps = readRDS("imputedSNPs.Jul23.RDS")

BAS = snps_haps %>%
	filter(poolroot == "BAS02")	%>%
# 	head(n=100) %>%
	mutate(result = pmap(list(ff,founderfreqs,cutree), flavor_SNP)) %>%
	unnest_wider(result) %>%
	mutate(biSNP = iSNP) %>%
	select(c(CHROM,POS,Nhaps,MAC,MACHAP,phase,biSNP)) 

bigsnp = snps_haps %>%
	filter(! poolroot %in% c("A1xB3_diploid_controls","BAS02r")) %>%
	select(c(CHROM,POS,pos,poolroot,iSNP)) %>%
	left_join(BAS) %>%
	mutate(dSNP = iSNP - biSNP) %>%
	select(-biSNP)
	
saveRDS(bigsnp,"imputedSNPs.final.Jul23.RDS")

