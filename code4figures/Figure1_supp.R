library(tidyverse)
library(Polychrome)
library(patchwork)

# build-in color palette
Glasbey = glasbey.colors(18)
Glasbey2 = Glasbey[-c(1,5)]
Glasbey2[16] = "#808080"
# UNR = UN Resolved
names(Glasbey2) = c("AB1","AB2","AB3","AB4","A5","A6","A8","A9","A11","B5","B7","B8","B9","A5B5","AB3B8","UNR")
swatch(Glasbey2)

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


xx_fig2D = temp2 %>%
	group_by(chr,drug,week) %>%
	summarize(avg_diff = mean(abs(R1-R2))) %>%
	ungroup() %>%
	mutate(chr2 = factor(chr, levels = CHR))

tiff("Figure_1_supp.tiff",width = 8, height = 8, units = 'in', res = 600)

ggplot(xx_fig2D, aes(x=week, y=avg_diff, color = drug)) +
	geom_line() +
	facet_wrap(vars(chr2)) +
    theme_light() +
    labs(x = "week", y = "avg freq change difference") +
	guides(color = guide_legend(override.aes = list(linewidth = 3 )))

dev.off()

