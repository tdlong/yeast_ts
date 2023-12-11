library(tidyverse)
xx=read_tsv("Sept5_22.allhaps.Ns50.txt.gz")
	
xx2 = xx %>%
	select(-c(nchr,mpos,NSNPs)) %>%
	filter(! poolroot %in% c("A1xB3_diploid_controls","BAS02r")) %>%
	separate_longer_delim(c(foundernames,cutree,founderfreqs), delim = ";") %>% 
	mutate(founderfreqs=as.numeric(founderfreqs)) %>%
	group_by(chr,pos,poolroot,cutree) %>%
	summarize(Nhaps=n(),frq=sum(founderfreqs),longhap=paste0(foundernames,collapse="")) %>%
	select(cutree)
	
temp = xx2 %>%
	ungroup() %>%
	filter(poolroot=="BAS02") %>%
	select(-poolroot) %>%
	mutate(BAS02frq = frq) %>%
	select(-frq)

xx3 = xx2 %>%
	ungroup() %>%
	filter(poolroot!="BAS02") %>%
	left_join(temp) %>%
	mutate(frqdiff = frq - BAS02frq) %>%
	select(-BAS02frq) %>%
	mutate(drug = substr(as.character(poolroot), 9, 10)) %>%
	mutate(rep = paste0(substr(as.character(poolroot), 9, 10),"_",substr(as.character(poolroot), 15, 16))) %>%
	mutate(week = substr(as.character(poolroot), 4, 5)) %>%
	mutate(week = as.numeric(week)) %>%
	filter(week %in% c(1,5,9,13,16,21))
	
write.table(xx3,"Sept5_22.allhaps.tidy.Ns50.txt")
## frqdiff = diff from BAS02, longhap is concatenated haps undistinguishable, Nhaps is Nhaps in longhap
## colnames = poolroot, chr, pos, drug, rep, week, frq, frqdiff, Nhaps, longhap	

xx3 = read.table("Sept5_22.allhaps.tidy.Ns50.txt")

library(tidyverse)
library(Polychrome)
library(patchwork)

xx3 = read.table("Sept5_22.allhaps.tidy.Ns50.txt")

# build-in color palette
Glasbey = glasbey.colors(18)
Glasbey2 = Glasbey[-c(1,5)]
Glasbey2[16] = "#808080"
# UNR = UN Resolved
names(Glasbey2) = c("AB1","AB2","AB3","AB4","A5","A6","A8","A9","A11","B5","B7","B8","B9","A5B5","AB3B8","UNR")
swatch(Glasbey2)

xx_fig2A = xx3 %>%
	filter(chr=="chrIV") %>%
	filter(drug %in% c("GA","NC")) %>%
	mutate(hap2 = longhap) %>%
	mutate(hap2 = case_when(
		hap2 %in% c("A11","A12","B11","A11A12","A11B11","A12B11","A11A12B11") ~ "A11",
		hap2 %in% c("AB2","A7","B12","AB2A7","AB2B12","A7B12","AB2A7B12") ~ "AB2",
		.default = as.character(hap2)
		)) %>%
	mutate(hap2 = case_when(
		!(hap2 %in% c("AB2","A11","A9","AB4","AB3","B9","B8","B5","A5","A6","AB1","B7","A8","A5B5","AB3B8")) ~ "UNR",
		.default = as.character(hap2)
		)) %>%
	filter(week %in% c(1,9,16,21)) %>%
	select(-longhap) %>%
	select(c(chr,pos,drug,rep,week,frq,frqdiff,hap2,Nhaps)) %>%
	group_by(chr,pos,drug,rep,week,hap2) %>%
	summarize(frq = sum(frq),frqdiff=sum(frqdiff),Nhaps=sum(Nhaps)) %>%
	ungroup() %>%
	mutate(hap3 = factor(hap2,levels=names(Glasbey2))) %>%
	arrange(chr,pos,drug,rep,week,desc(frqdiff)) %>%
	group_by(chr,pos,drug,rep,week) %>%
	slice(1:2) %>%
	ungroup()

write.table(xx_fig2A,"xx_fig2A.txt")


####  now dig in check subset are OK, and make figures.

# most INCREASED haplotype
F2A = ggplot(xx_fig2A,
	aes(x=pos/1000, y=frqdiff, colour = hap3 )) +
    geom_point(size=0.05) +
    ylim(0, 1) +
    theme_light() +
    facet_grid(rep~week) +
    labs(x = "chrIV (kb)", y = "freq change (MCH)") +
    scale_colour_manual(values=Glasbey2) +
	theme(legend.title=element_blank()) +
	guides(colour = guide_legend(override.aes = list(size=5)))


#	theme(legend.position = "bottom",
#		legend.title=element_blank(),
#		legend.text = element_text(size = 5)) +
#	guides(color = guide_legend(override.aes = list(size = 1), nrow=2, byrow=TRUE))


## Now ... make other figures using the same rules....colors...

## This is the correlation over Most changed HAPs (most changed on average over two
## replicates)
## (this will become the bottom portion of Figure 1)

xx3 = read.table("Sept5_22.allhaps.tidy.Ns50.txt")
temp = xx3 %>%
#	filter(chr=="chrI" & pos==32800 & drug=="CD" & week==1) %>%
	filter(chr != "chrM") %>%
	mutate(hap2 = longhap) %>%
	mutate(hap2 = case_when(
		hap2 %in% c("A11","A12","B11","A11A12","A11B11","A12B11","A11A12B11") ~ "A11",
		hap2 %in% c("AB2","A7","B12","AB2A7","AB2B12","A7B12","AB2A7B12") ~ "AB2",
		.default = as.character(hap2)
		)) %>%
	mutate(hap2 = case_when(
		!(hap2 %in% c("AB2","A11","A9","AB4","AB3","B9","B8","B5","A5","A6","AB1","B7","A8","A5B5","AB3B8")) ~ "UNR",
		.default = as.character(hap2)
		)) %>%
	select(-longhap) %>%
	select(c(chr,pos,drug,rep,week,frq,frqdiff,hap2,Nhaps)) %>%
	group_by(chr,pos,drug,rep,week,hap2) %>%
	summarize(frq = sum(frq),frqdiff=sum(frqdiff),Nhaps=sum(Nhaps)) %>%
	ungroup()

#temp %>%
#	filter(chr=="chrIII" & pos==116800 & drug=="CD" & week == 1)

temp2 = temp %>%
	group_by(chr, pos, drug, week, hap2) %>%
	mutate(simple_rep = as.numeric(as.factor(rep))) %>%
	mutate(simple_rep = paste0("R", simple_rep)) %>%
	ungroup() %>%
	select(-c(rep,frq,Nhaps)) %>%
	pivot_wider(names_from = simple_rep, values_from = frqdiff) %>%
	mutate(mfrqdiff = (R1 + R2)/2) %>%
	group_by(chr,pos,drug,week) %>%
	arrange(desc(mfrqdiff)) %>%
	slice(1)

CHR=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
xx_3ex = temp2 %>%
	filter(chr %in% c("chrIV", "chrIX", "chrXIII")) %>%
	group_by(chr,drug,week) %>%
	summarize(avg_diff = mean(abs(R1-R2))) %>%
	ungroup()

xx_genome = temp2 %>%
	group_by(drug,week) %>%
	summarize(avg_diff = mean(abs(R1-R2))) %>%
	mutate(chr = "genome")

xx_fig2C = xx_3ex %>%
	rbind(xx_genome) %>%
	mutate(chr2 = factor(chr, levels=c("genome","chrIV", "chrIX", "chrXIII")))

F2B = ggplot(xx_fig2C, aes(x=week, y=avg_diff, color = drug)) +
	geom_line() +
	facet_wrap(vars(chr2),nrow=1) +
    theme_light() +
    labs(x = "week", y = "avg freq change diff") +
	guides(color = guide_legend(override.aes = list(linewidth = 3 )))

library(patchwork)
patched = F2A + F2B + plot_layout(ncol = 1, heights = c(3, 1)) + plot_annotation(tag_levels = 'A')

tiff("Figure_2.tiff",width = 7.5, height = 7, units = 'in', res = 600)
patched
dev.off()

