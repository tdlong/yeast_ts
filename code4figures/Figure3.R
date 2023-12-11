library(Polychrome)
# build-in color palette
Glasbey = glasbey.colors(18)
Glasbey2 = Glasbey[-c(1,5)]
Glasbey2[16] = "#808080"
# UNR = UN Resolved
names(Glasbey2) = c("AB1","AB2","AB3","AB4","A5","A6","A8","A9","A11","B5","B7","B8","B9","A5B5","AB3B8","UNR")
swatch(Glasbey2)

xx3 = read.table("Sept5_22.allhaps.tidy.Ns50.txt")
temp1_som = xx3 %>%
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
	filter(week %in% c(1,5,9,13,16,21)) %>%
	select(-longhap) %>%
	select(c(chr,pos,drug,rep,week,frq,frqdiff,hap2,Nhaps)) %>%
	group_by(chr,pos,drug,rep,week,hap2) %>%
##  This summarize is adding up freqs within a haplotype pool!
	summarize(frq = sum(frq),frqdiff=sum(frqdiff),Nhaps=sum(Nhaps)) %>%
	ungroup() %>%
	mutate(hap3 = factor(hap2,levels=names(Glasbey2)))	


#  I should add time zero
temp2_som = temp1_som %>%
	filter(week==1) %>%
	mutate(week = 0) %>%
	mutate(frq = frq - frqdiff) %>%
	mutate(frqdiff = 0)


####
myweeks = c(1,5,9,16)
somlist = list()
for(i in 1:4){
	somlist[[i]] = temp1_som %>%
		filter(week==myweeks[i]) %>%
		mutate(pos10kb = round(pos/25000,0)) %>%
		arrange(chr,pos10kb,drug,rep,desc(abs(frqdiff))) %>%
		group_by(chr,pos10kb,drug,rep) %>%
		slice(1) %>%
		ungroup() %>%
		mutate(weekmost = myweeks[i])
	}

### old
#most_som = do.call(rbind, somlist) %>%
#	arrange(chr,drug,rep,pos,weekmost) %>%
#	select(c(chr,drug,rep,pos,weekmost,hap3)) %>% 
#	left_join(temp1_som, 
#		join_by(chr == chr, drug == drug, rep == rep, pos == pos, weekmost == week, hap3 == hap3)
#		) %>%	
#	mutate(lab=paste(chr,pos,drug,rep,weekmost,sep="-")) %>%
#	select(lab,week,frqdiff) %>%
#	pivot_wider(names_from=week,values_from=frqdiff) %>%
#	remove_rownames() %>%
#	column_to_rownames(var = "lab")

### new
most_som = do.call(rbind, somlist) %>%
	arrange(chr,drug,rep,pos,weekmost) %>%
	select(c(chr,drug,rep,pos,weekmost,hap3)) %>% 
	left_join(bind_rows(temp1_som,temp2_som)) %>%	
	mutate(lab=paste(chr,pos,drug,rep,weekmost,sep="-")) %>%
	select(lab,week,frq) %>%
	pivot_wider(names_from=week,values_from=frq) %>%
	relocate("0", .after = lab) %>%
	remove_rownames() %>%
	column_to_rownames(var = "lab")

		
library(som)
out = som(most_som, xdim=5, ydim=5, topol="rect", neigh="gaussian")

som_plot = most_som %>%
	mutate(x = out$visual$x) %>%
	mutate(y = out$visual$y) %>%
	rownames_to_column(var = "lab") %>%
	separate(lab, c("chr", "pos","drug","rep","weekmost"),sep="-",remove=FALSE) %>%
	select(-c(chr,pos)) %>%
	pivot_longer(cols=5:10, names_to='week', values_to='frq') %>%
	mutate(week = as.numeric(week))

# unlike ylim, under coord_cartesian, smoothing respects ALL datapoints
PP = ggplot(som_plot, aes(x=week, y=frq)) +
	geom_line(aes(group = lab), linewidth=0.2, alpha=0.01) +
	geom_smooth(method = "loess", span = 0.25) +
	theme_light() +
	labs(x = "week", y = "MCH freq") +
	coord_cartesian(ylim = c(0, 1)) +
	theme(axis.text.x=element_blank()) +
	facet_grid(x~y, labeller = label_both)
	
#	guides(colour = guide_legend(override.aes = list(linewidth=3,alpha=1)))
		
temp_bar = som_plot %>%
	group_by(x,y,drug,weekmost) %>%
	summarize(count=n()) %>%
	mutate(count = ifelse(is.na(count), 0, count))
	
PB = ggplot(temp_bar, aes(x=drug, y=count, fill=factor(weekmost,levels=c(c("1", "5", "9", "16"))))) +
	geom_bar(stat='identity') +
	theme_light() +
	facet_grid(x~y, labeller = label_both) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	guides(fill=guide_legend(title="@week"))
#	theme(axis.text.x=element_blank())


library(patchwork)

patched = PP + PB +
	plot_layout(ncol = 1) +
	plot_annotation(tag_levels = 'A')

tiff("Figure_SOM_pool.tiff",width = 7.5, height = 10, units = 'in', res = 600)
patched
dev.off()   

