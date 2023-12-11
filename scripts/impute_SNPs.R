library(tidyverse)

#zcat Sept5_22.allhaps.Ns50.txt.gz | head
# poolroot	chr	nchr	pos	mpos	NSNPs	foundernames	cutree	founderfreqs
# A1xB3_diploid_controls	chrI	1	32800	32800	1144	AB1;AB2;AB3;AB4;A5;A6;A7;A8;A9;A11;A12;B5;B7;B8;B9;B11;B12	1;2;3;4;5;6;7;8;9;10;10;11;12;13;14;10;2	3e-04;3e-04;3e-04;0.0742;3e-04;0.0635;0.1477;3e-04;3e-04;3e-04;3e-04;3e-04;3e-04;0.2192;0.3234;3e-04;0.1688

haps = read_tsv("Sept5_22.allhaps.Ns50.txt.gz") %>%
#	filter(chr=="chrVI") %>%
	select(-c(nchr,mpos,NSNPs,foundernames)) 
	
haps_small = haps %>%
	select(c(chr,pos)) %>%
	arrange(chr,pos) %>%
	unique() %>%
	mutate(posL = pos, posR = pos) %>%
	select(-pos)	
	
snps = as_tibble(read.table("results/bcf2/SNP.freq.txt")) %>%
	select(1:20) %>%
#	filter(CHROM=="chrVI") %>%
	select(-B6) %>%
	unite(ff, -c(CHROM,POS),sep = ";")

blah="AB1;AB2;AB3;AB4;A5;A6;A7;A8;A9;A11;A12;B5;B7;B8;B9;B11;B12"
templabs = unlist(strsplit(blah, split=';'))

impute_SNP = function(ff, founderfreq){
	ffv = as.numeric(unlist(strsplit(ff, split=';')))	
	fpv = as.numeric(unlist(strsplit(founderfreq, split=';')))
	sum(ffv*fpv)
	}

flavor_SNP = function(ff, founderfreq, cutree, templabs){
	ffv = as.numeric(unlist(strsplit(ff, split=';')))	
	fpv = as.numeric(unlist(strsplit(founderfreq, split=';')))
	ctv = as.numeric(unlist(strsplit(cutree, split=';')))	
	cc=c()
	dd=c()   # names of each group
	tf=c()   # freq of each group
	obsf = c()
	Ndiffhaps = 0
	for(i in 1:17){
		if(! ctv[i] %in% cc){
			cc = append(cc, ctv[i])
			Ndiffhaps = Ndiffhaps + 1
			dd = append(dd, templabs[i])
			tf = append(tf, ffv[i])
			obsf = append(obsf, fpv[i])
			} else {
			obsf[ctv[i]] = obsf[ctv[i]] + fpv[i]
			}
		}
	list(foLabs=dd,foStates=tf, foFreqs=obsf)
	}

snps_haps = snps %>% 
	left_join(haps_small %>% select(-posR), join_by(CHROM==chr, closest(POS >= posL))) %>%
	left_join(haps_small %>% select(-posL), join_by(CHROM==chr, closest(POS <  posR))) %>%
	mutate(pos = ifelse(abs(POS-posL) < abs(POS-posR),posL, posR)) %>%
	mutate(pos = ifelse(is.na(pos),c(posL,posR)[!is.na(c(posL,posR))],pos)) %>%
	select(-c(posL,posR)) %>%
	left_join(haps, join_by(CHROM==chr,pos),multiple = "all") %>%
	mutate(iSNP = map2_dbl(ff, founderfreqs, impute_SNP)) 
# %>%
#	mutate(result = pmap(~flavor_SNP(.ff, .founderfreqs, .cutree, templabs))) %>%
#	unnest_wider(result) 
	
saveRDS(snps_haps,"imputedSNPs.Jul23.RDS")
