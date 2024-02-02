library(tidyverse)
library(data.table)
library(dtplyr)
args = commandArgs(trailingOnly=TRUE)
mychr=as.character(args[1])
# mychr="chrII"

filein=paste0("results/process/RefAlt.",mychr,".txt")
rdsfile=paste0("results/process/newSNP.",mychr,".rds")


founders=c("A11_YJM978", "A12_YJM975", "A5_DBVPG6765", "A6_BC187_b",     
	"A7_SK1", "A8_L_1374", "A9_UWOPSO3_461_4", "AB1b",            
	"AB2b", "AB3b", "AB4b", "B11_YJM981",     
	"B12_Y55", "B5_273614N", "B6_YPS606", "B7_L_1528",       
	"B8_UWOPS83_787_3", "B9_UWOPS87_2421")
	
ts = c("SEE01B02CD600R04", "SEE01B02CD600R06", "SEE01B02CP010R03", "SEE01B02CP010R04", "SEE01B02DI027R01", "SEE01B02DI027R06", "SEE01B02GA120R01", "SEE01B02GA120R05", "SEE01B02NC950R01", "SEE01B02NC950R02",
"SEE01B02YP000R03", "SEE01B02YP000R04", "SEE05B02CD600R04", "SEE05B02CD600R06", "SEE05B02CP010R03", "SEE05B02CP010R04", "SEE05B02DI027R01", "SEE05B02DI027R06", "SEE05B02GA120R01", "SEE05B02GA120R05",
"SEE05B02NC950R01", "SEE05B02NC950R02", "SEE05B02YP000R03", "SEE05B02YP000R04", "SEE09B02CD600R04", "SEE09B02CD600R06", "SEE09B02CP010R03", "SEE09B02CP010R04", "SEE09B02DI027R01",
"SEE09B02DI027R06", "SEE09B02GA120R01", "SEE09B02GA120R05", "SEE09B02NC950R01", "SEE09B02NC950R02", "SEE09B02YP000R03", "SEE09B02YP000R04", "SEE13B02CD600R04", "SEE13B02CD600R06", "SEE13B02CP010R03",
"SEE13B02CP010R04", "SEE13B02DI027R01", "SEE13B02DI027R06", "SEE13B02GA120R01", "SEE13B02GA120R05", "SEE13B02NC950R01", "SEE13B02NC950R02", "SEE13B02YP000R03", "SEE13B02YP000R04", "SEE16B02CD600R04",
"SEE16B02CD600R06", "SEE16B02CP010R03", "SEE16B02CP010R04", "SEE16B02DI027R01", "SEE16B02DI027R06", "SEE16B02GA120R01", "SEE16B02GA120R05", "SEE16B02NC950R01", "SEE16B02NC950R02", "SEE16B02YP000R03",
"SEE16B02YP000R04", "SEE21B02CD600R04", "SEE21B02CD600R06", "SEE21B02CP010R03", "SEE21B02CP010R04", "SEE21B02DI027R01", "SEE21B02DI027R06", "SEE21B02GA120R01", "SEE21B02GA120R05", "SEE21B02NC950R01",
"SEE21B02NC950R02", "SEE21B02YP000R03", "SEE21B02YP000R04")


df = lazy_dt(read.table(filein,header=TRUE))
df2 = df %>%
	pivot_longer(c(-CHROM,-POS), names_to = "lab", values_to = "count") %>%
	mutate(RefAlt = str_sub(lab,1,3)) %>%
	mutate(name = str_sub(lab,5)) %>%
	select(-lab) %>%
#	separate(lab, c("RefAlt", "name"), "_", extra = "merge") %>%
	pivot_wider(names_from = RefAlt, values_from = count) %>%
	mutate(freq = REF/(REF+ALT), N = REF+ALT) %>%
	select(-c("REF","ALT")) %>%
	mutate(name = str_remove(name, fixed(".dfs7.adl.tdlong.RL_yeast.DipX_QTL.May1.bam."))) %>%
	mutate(name = str_remove(name, fixed(".bam"))) %>%
	mutate(name = str_remove(name, fixed("results.mergebam."))) %>%
	filter((name %in% founders | name %in% ts | name == "BAS02")) %>%
	as_tibble()
	

# identify SNPs that are fixed in founders
good_SNPs = df2 %>%
	filter(name %in% founders) %>%
	group_by(CHROM,POS) %>%
	summarize(zeros=sum(N==0),allRef=sum(N!=0 & freq < 0.05),allAlt=sum(N!=0 & freq > 0.95)) %>%
	ungroup() %>%
	# non-informative SNPs are potentially new mutations
	filter(zeros==0 & ( allRef == length(founders) | allAlt == length(founders) )) %>%
	select(c(CHROM,POS))

# now subset the entire dataset for:
# 	- rare in BAS
#	- at least 5% in one of the 72 focus populations
df3 = good_SNPs %>% 
	left_join(df2, multiple = "all") %>%
	filter(! (name %in% founders)) %>%
	select(-N) %>%
	pivot_wider(names_from="name", values_from="freq") %>%
	filter(BAS02 < 0.01 | BAS02 > 0.99) %>%
	rowwise() %>%
	mutate(maxfreq=max(c_across(starts_with("SEE")))) %>%
	filter((maxfreq*(1-maxfreq)) >= 0.05)

if(nrow(df3)>=1){saveRDS(df3, file = rdsfile)}
