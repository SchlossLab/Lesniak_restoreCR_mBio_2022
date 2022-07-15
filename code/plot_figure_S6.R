###############################################################################
#
# Figure S6
#	main idea - community at the time of challenge for cefoperazone-treated mice
# 	plot - relative abundace of otus present in cef mice 
#
# input files:
#   code/functions.R
#	data/process/restore_metadata_clean.tsv
#	data/mothur/sample.final.0.03.subsample.shared
#	data/process/restore_taxonomy_clean.tsv
#
# output files:
#	submission/Figure_S6.tiff
#
# Nick Lesniak 2022-02-25
###############################################################################
# setup environment
###############################################################################
source(here::here('code/functions.R'))

###############################################################################
# load data
###############################################################################
taxonomy <- read_taxonomy()

shared <- read_shared() %>% 
	select(-label, -numOtus) %>% 
	pivot_longer(-Group, names_to = 'OTU', values_to = 'counts') %>% 
	mutate(relative_abundance = counts/2480 * 100)%>% 
	left_join(taxonomy, by = 'OTU')

cef_metadata <- read_metadata() %>% 
	filter(antibiotic == 'Cefoperazone',
		treatment %in% c('-1', '-2'))

###############################################################################
# setup data
###############################################################################
cef_shared <- cef_metadata %>% 
	filter(day == 0) %>% 
	select(treatment, sample_id) %>% 
	left_join(shared, by = c('sample_id' = 'Group'))

otus_present <- cef_shared %>% 
	group_by(OTU) %>% 
	summarise(n = sum(counts > 1)) %>% 
	filter(n > 0) %>% 
	pull(OTU)

cef_shared_plot_data <- cef_shared %>% 
	filter(OTU %in% otus_present) %>% 
	mutate(taxa_label = tax_otu_label,
		otu_order = paste0(gsub(' \\(.*', '', taxa_label), OTU),
		relative_abundance = ifelse(relative_abundance == 0, 0.035, 
			relative_abundance),
		taxa_label = gsub('_', ' ', taxa_label),
		taxa_label = paste0('*', taxa_label),
		taxa_label = gsub(' \\(OTU','* (OTU', taxa_label),
		taxa_label = ifelse(grepl('unclassified', taxa_label), 
			gsub(' unclassified\\*', '**', taxa_label),
			taxa_label),
		treatment = ifelse(treatment == '-1', 'FCT', '1:10^1'))

###############################################################################
# plot data
###############################################################################
cef_shared_plot <- cef_shared_plot_data %>% 
	ggplot(aes(x = reorder(taxa_label, desc(otu_order)), y = relative_abundance, 
			shape = treatment)) +
		geom_point(position=position_jitterdodge(
				dodge.width = 0.5,
				jitter.width = 0.2),
			color = '#3A9CBC') + 
		geom_hline(yintercept = 1/24.8, linetype = 'dashed', lwd = 0.3) + 
		theme_bw() + 
		labs(x = NULL, y = 'Relative Abundance (%)', shape = NULL) + 
		theme(
			axis.text.y = ggtext::element_markdown(),
			panel.grid.major.y = element_line(color = 'gray95'),
			panel.grid.major.x = element_line(color = 'gray85'),
			panel.grid.minor.x = element_blank(),
			legend.text = ggtext::element_markdown(),
			legend.position = 'top',
			legend.margin=margin(2,0,0,0),
	    	legend.box.margin=margin(2,-10,-10,-10)) + 
		guides(shape = guide_legend(reverse = TRUE)) + 
		scale_shape_manual(values = c(16,1)) + 
		scale_y_log10(#lim = c(0.04,20),
			breaks=c(0.01, 0.1, 1, 10, 100),
			labels=c("0","0.1","1","10","100")) +
		coord_flip()

###############################################################################
# load data
###############################################################################
ggsave(here('submission/Figure_S6.tiff'),
	cef_shared_plot,
	height = 8, width = 5, unit = 'in',
	compression = 'lzw')

system(paste("convert -compress lzw", 
		here('submission/Figure_S6.tiff'), 
		here('submission/Figure_S6.tiff')))
###############################################################################