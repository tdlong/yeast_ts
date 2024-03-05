library(tidyverse)
df <- list.files(path="results/process", pattern = ".rds",full.names=TRUE) %>%
  map(readRDS) %>% 
  bind_rows()
# df has 26 rows...
write.table(df,"denovo.SNPs.txt")

df2 = df %>%
	select(-maxfreq) %>%
	pivot_longer(!c(CHROM,POS), names_to = "name", values_to = "freq") %>%
	filter(name != "BAS02") %>%   # add me back...
	# SEE01B02NC950R02
	mutate(gen = as.numeric(str_sub(name, 4, 5))) %>%
	mutate(drug_rep = paste0(str_sub(name,9,10),"_",str_sub(name,15,16))) %>%
	select(-name) %>%
	relocate(drug_rep, .after = POS) %>%
	relocate(gen, .after = drug_rep)

drug_rep = as.character(levels(as.factor(df2$drug_rep)))

dfbas = df %>%
	select(-maxfreq) %>%
	pivot_longer(!c(CHROM,POS), names_to = "name", values_to = "freq") %>%
	filter(name == "BAS02") %>%
	mutate(gen = 0) %>%
	select(-name) %>%
	expand(nesting(CHROM,POS,gen,freq),drug_rep) %>%
	relocate(drug_rep, .after = POS) %>%
	relocate(gen, .after = drug_rep)
	
df3 = bind_rows(df2,dfbas) %>%
	arrange(CHROM,POS,drug_rep,gen)	
	
write.table(df3,"denovo.SNPs.2.txt")

### now make some plots.
#  CHROM    POS drug_rep   gen    freq
#   <chr>  <int> <chr>    <dbl>   <dbl>
# 1 chrI  228756 CD_04        0 0.00433
# 2 chrI  228756 CD_04        1 0      
# 3 chrI  228756 CD_04        5 0      

df3 = read.table("denovo.SNPs.2.txt") %>%
	mutate(drug = str_sub(drug_rep,1,2)) %>%
	mutate(CHROM_POS = paste0(CHROM,"_",POS))
	
df4 = df3 %>%
	filter(CHROM_POS != "chrXIV_642529") %>%
	filter(freq > 0.001) %>%
	rename(week=gen)

NewMut = ggplot(df4, aes(x=week,y=log10(freq), group=drug_rep)) +
#	geom_line(aes(color = drug)) +
	geom_point(aes(color = drug)) + 
	facet_wrap(~ CHROM_POS)
	
tiff("Figure_2_supp.tiff", units="in", width=7, height=9, res=600)
NewMut
dev.off()

# almost certainly NOT a new mutation!!
# ggplot(df3 %>% filter(CHROM_POS == "chrXIV_642529"),aes(x=gen,y=freq, group=drug_rep)) +
#	geom_line(aes(color = drug)) +
#	geom_point(aes(color = drug))

