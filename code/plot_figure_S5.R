###############################################################################
#
# Figure S4
#	main idea - 
# 	plot - 
#
# Nick Lesniak 2022-01-12
###############################################################################
# setup environment
###############################################################################
source(here::here('code/functions.R'))

###############################################################################
# load data
###############################################################################
metadata <- read_metadata() %>% 
	select(sample_id, mouse_id)

taxonomy <- read_taxonomy()

shared <- read_shared() %>% 
	select(-label, -numOtus) %>% 
	pivot_longer(-Group, names_to = 'OTU', values_to = 'counts') %>% 
	mutate(relative_abundance = counts/2480 * 100)%>% 
	left_join(taxonomy, by = 'OTU') %>% 
	select(Group, OTU, relative_abundance, taxa_label = tax_otu_label)

lefse_design <- data.table::fread(here(paste0('data/process/lefse/strep_antibiotic_OTU.design')))
	colnames(lefse_design) <- c('Group', 'class')

lefse_df <- data.table::fread(here(paste0('data/process/lefse/strep_antibiotic_OTU.0.03.lefse_summary')),
		fill = T) %>% 
	filter(!is.na(LDA)) %>% 
	pull(OTU)

###############################################################################
#  setup data
###############################################################################

plot_data <- shared %>% 
	filter(OTU %in% lefse_df) %>% 
	right_join(lefse_design, by = c('Group')) %>% 
	left_join(metadata, by = c('Group' = 'sample_id')) %>% 
	#group_by(mouse_id, OTU, tax_otu_label, class) %>% # removed mouse_id from interactive plot
	#summarise(relative_abundance = sum(relative_abundance)) %>% 
	select(mouse_id, taxa_label, OTU, class, relative_abundance) %>% 
	mutate(otu_order = paste0(gsub(' \\(.*', '', taxa_label), OTU),
		class = factor(class, levels = c('initial', 'post_antibiotic')),
		taxa_label = gsub('_', ' ', taxa_label),
		taxa_label = paste0('*', taxa_label),
		taxa_label = gsub(' \\(OTU','* (OTU', taxa_label),
		taxa_label = ifelse(grepl('unclassified', taxa_label), 
			gsub(' unclassified\\*', '**', taxa_label),
			taxa_label)) %>% 
	pivot_wider(names_from = 'class', values_from = 'relative_abundance') %>% 
	group_by(taxa_label, otu_order) %>% 
		mutate(median_initial = median(initial, na.rm = T),
			median_post_antibiotic = median(post_antibiotic, na.rm = T)) %>% 
	pivot_longer(cols = c('initial', 'post_antibiotic'), 
		names_to = 'outcome', values_to = 'relative_abundance') %>% 
	mutate(high_low = median_initial > 0.12 | median_post_antibiotic > 0.1)

###############################################################################
#  plot data
###############################################################################
# plot difference between initial and streptomycin treatment

create_plot <- function(plot_dataframe){
	plot_dataframe %>% 
		ggplot(aes(reorder(taxa_label, desc(otu_order)))) + 
			# limit of detection
			geom_hline(yintercept = 1/24.8, linetype = 'dashed', lwd = 0.3) + 
			#  barbell geom
			geom_segment(aes(y = median_initial, yend = median_post_antibiotic, 
				xend = reorder(taxa_label, desc(otu_order))),
				arrow = arrow(type = 'closed', angle = 10), color = 'darkgrey') +
			geom_point(aes(y = median_initial), color = 'black', size = 3) + 
			geom_point(aes(y = median_post_antibiotic), color = '#D37A1F', size = 3) + 
			geom_point(aes(y = relative_abundance, color = outcome), 
				alpha = 0.1, position = position_dodge(width = 0.8)) + 
			scale_color_manual(values = c('black', '#D37A1F'),
				breaks = c('initial', 'post_antibiotic'),
				labels = c('Initial', 'After streptomycin')) + 
			# plot layout
			coord_flip() + theme_bw() + 
			labs(x = NULL, y = 'Relative Abundance (%)', color = NULL) + 
			guides(color = guide_legend(override.aes = list(linetype = 0, 
					size = 2, alpha = 1))) +
			theme(panel.grid.minor.x = element_blank(),
				legend.position = 'top', 
				legend.justification='left',
	        	legend.direction='horizontal',
				axis.text.y = ggtext::element_markdown(),
				axis.text.x = ggtext::element_markdown())
}

strep_antibiotic_high_plot <- plot_data %>% 
	filter(high_low) %>% 
	create_plot() + 
		scale_y_log10(limits = c(0.03,100),
			breaks = c(0.01, 0.1, 1, 10, 100),
			labels = c('10^-2', '10^-1', '10^0', '10^1', '10^2')) +
		theme(legend.position = 'none')
	
strep_antibiotic_low_plot <- plot_data %>% 
	filter(!high_low) %>% 
	create_plot() + 
		scale_y_log10(limits = c(0.03,10),
			breaks = c(0.01, 0.1, 1, 10, 100),
			labels = c('10^-2', '10^-1', '10^0', '10^1', '10^2'))

###############################################################################
#  save plot
###############################################################################
ggsave(here('submission/Figure_S5.tiff'),
	cowplot::plot_grid(strep_antibiotic_high_plot,
		strep_antibiotic_low_plot,
		nrow = 1),
	height = 11, width = 10, unit = 'in',
	compression = 'lzw')

###############################################################################
