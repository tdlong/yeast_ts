bigsnp = readRDS("imputedSNPs.final.Jul23.RDS")

library(Polychrome)
# build-in color palette
Glasbey = glasbey.colors(18)
Glasbey2 = Glasbey[-c(1,5)]
Glasbey2[16] = "#808080"
# UNR = UN Resolved
names(Glasbey2) = c("AB1","AB2","AB3","AB4","A5","A6","A8","A9","A11","B5","B7","B8","B9","A5B5","AB3B8","UNR")
swatch(Glasbey2)

xx2 = bigsnp %>%
	ungroup() %>%
	mutate(drug = substr(as.character(poolroot), 9, 10)) %>%
	mutate(rep = paste0(substr(as.character(poolroot), 9, 10),"_",substr(as.character(poolroot), 15, 16))) %>%
	mutate(week = substr(as.character(poolroot), 4, 5)) %>%
	mutate(week = as.numeric(week)) %>%
	filter(week %in% c(1,5,9,13,16,21)) %>%
	filter(MAC != 0)
	
table(xx2$MAC)
#      1       2       3       4       5       6       7 
#5863464 1289376  862200  978552 1044864  333936     288 

#
#57% of SNPs are "singletons"... where singleton mean minor allele is private to a single founder haplotype
#	& that haplotype occurs at a "measurable" frequency in BAS

xx_fig2B = xx2 %>%
	filter(drug %in% c("GA","NC")) %>%
	filter(MAC == 1) %>%
	mutate(phase_dSNP = phase*dSNP) %>%
	mutate(hap3 = factor(MACHAP,levels=names(Glasbey2))) %>%
	filter(week %in% c(1,9,16,21) & CHROM=="chrIV" & drug %in% c("GA","NC"))

table(xx_fig2B$Nhaps)	

#    6     7     8     9    10    11    12    13 
#   96   848  3440 11376 32016 51296 48288 17792 


F3A = ggplot(xx_fig2B, 
	aes(x=POS/1000, y=phase_dSNP, color = hap3)) +
    geom_point(size=0.04) +
    theme_light() +
    facet_grid(rep~week) +
    labs(x = "chrIV (kb)", y = "phased frq change") +
    scale_colour_manual(values=Glasbey2) +
    theme(legend.position = "none")
    
    
#	theme(legend.position = "bottom",
#		legend.title=element_blank(),
#		legend.text = element_text(size = 5)) +
#	guides(color = guide_legend(override.aes = list(size = 1), nrow=2, byrow=TRUE))



####   same figure straight up SNP differences
# (scp .../results/bcf2/SNP.freq.txt)

rawsnps = as_tibble(read.table("SNP.freq.txt")) %>%
	select(-c(3:20)) %>%       # founders
	pivot_longer(-c(1:2),names_to="poolroot",values_to="frq")

rawBAS = rawsnps %>%
	filter(poolroot == "BAS02")	%>%
	mutate(bfrq = frq) %>%
	select(c(CHROM,POS,bfrq)) 

rawsnp2 = rawsnps %>%
	filter(! poolroot %in% c("A1xB3_diploid_controls","BAS02r")) %>%
	left_join(rawBAS) %>%
	mutate(dfrq = frq - bfrq) %>%
	select(-bfrq) %>%
	mutate(drug = substr(as.character(poolroot), 9, 10)) %>%
	mutate(rep = paste0(substr(as.character(poolroot), 9, 10),"_",substr(as.character(poolroot), 15, 16))) %>%
	mutate(week = substr(as.character(poolroot), 4, 5)) %>%
	mutate(week = as.numeric(week)) %>%
#	filter(week %in% c(1,5,9,13,16,21)) %>%
	filter(week %in% c(1,9,16,21)) %>%
	filter(drug %in% c("GA","NC")) %>%
	filter(CHROM=="chrIV")

F3B = ggplot(rawsnp2, aes(x=POS/1000, y=abs(dfrq), alpha=abs(dfrq))) +
    geom_point(size=0.04) +
    theme_light() +
    facet_grid(rep~week) +
    labs(x = "chrIV (kb)", y = "abs(frq change)") +
	theme(legend.position = "none")

library(patchwork)

patched = F3B + F3A +
	plot_layout(ncol = 1) +
	plot_annotation(tag_levels = 'A')

tiff("Figure_4_SNP.tiff",width = 7.5, height = 10, units = 'in', res = 600)
patched
dev.off()   

