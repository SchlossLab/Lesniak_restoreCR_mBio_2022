###############################################################################
#
# Figure 5
#	main idea - 
# 	plot - 
#
# Nick Lesniak 2022-01-12
###############################################################################
# setup environment
###############################################################################
source(here::here('code/functions.R'))

plot_lefse <- function(file_name){
	lefse_design <- data.table::fread(here(paste0('data/process/lefse/', 
		file_name, '.design')))
		colnames(lefse_design) <- c('Group', 'class')

	lefse_df <- data.table::fread(here(paste0('data/process/lefse/', file_name, '.0.03.lefse_summary')),
		fill = T) %>% 
	filter(!is.na(LDA)) %>% 
	mutate(fctr_class = factor(Class),
		order = LDA + 10 * as.numeric(fctr_class))

	plot_data <- shared %>% 
		right_join(lefse_df, by = c('OTU' = 'OTU')) %>% 
		right_join(lefse_design, by = c('Group')) %>% 
		group_by(Group, OTU, tax_otu_label, order, class) %>% # removed mouse_id from interactive plot
		summarise(relative_abundance = sum(relative_abundance)) %>% 
		mutate(taxa_label = tax_otu_label,
			otu_order = paste0(gsub(' \\(.*', '', taxa_label), OTU),
			relative_abundance = ifelse(relative_abundance == 0, 0.035, 
				relative_abundance),
			class = factor(class, levels = c('initial', 'post_antibiotic',
				'uncolonized', 'cleared', 'colonized')),
			taxa_label = gsub('_', ' ', taxa_label),
			taxa_label = paste0('*', taxa_label),
			taxa_label = gsub(' \\(OTU','* (OTU', taxa_label),
			taxa_label = ifelse(grepl('unclassified', taxa_label), 
				gsub(' unclassified\\*', '**', taxa_label),
				taxa_label))

	plot_data %>% 
		ggplot(aes(x = reorder(taxa_label, desc(otu_order)), y = relative_abundance, 
				color = class)) +
			stat_summary(fun.data = 'median_hilow', aes(group = class),
				fun.args = (conf.int=0.5), position = position_dodge(width = .5)) +
			geom_hline(yintercept = 1/24.8, linetype = 'dashed', lwd = 0.3) + 
			theme_bw() + 
			theme(
				axis.text.y = ggtext::element_markdown(),
				panel.grid.major.y = element_line(color = 'gray95'),
				panel.grid.major.x = element_line(color = 'gray85'),
				panel.grid.minor.x = element_blank(),
				legend.position = 'top',
				legend.margin=margin(2,0,0,0),
    	    	legend.box.margin=margin(2,-10,-10,-10)) + 
			guides(color = guide_legend(reverse = TRUE))
}

###############################################################################
# load data
###############################################################################
metadata <- read_metadata()

taxonomy <- read_taxonomy()

shared <- read_shared() %>% 
	select(-label, -numOtus) %>% 
	pivot_longer(-Group, names_to = 'OTU', values_to = 'counts') %>% 
	mutate(relative_abundance = counts/2480 * 100)%>% 
	left_join(taxonomy, by = 'OTU')

###############################################################################
#  setup data
###############################################################################

###############################################################################
#  analyze data
###############################################################################

###############################################################################
#  plot data
###############################################################################
strep_day0_clearance_OTU_plot <- plot_lefse('strep_day0_clearance_OTU') + 
			scale_color_manual(
				breaks = c('uncolonized', 'cleared', 'colonized'),
				labels = c('', 'Cleared colonization', 'Remain colonized'),
				values = c("#7c8a4f", "#c2dcb8", "#096013"))
strep_day10_clearance_OTU_plot <- plot_lefse('strep_day10_clearance_OTU') + 
	theme(legend.position = 'none') + 
			scale_color_manual(
				breaks = c('uncolonized', 'cleared', 'colonized'),
				labels = c('', 'Cleared colonization', 'Remain colonized'),
				values = c("#7c8a4f", "#c2dcb8", "#096013"))
strep_day0_colonization_OTU_plot <- plot_lefse('strep_day0_colonization_OTU') + 
			scale_color_manual(
				breaks = c('uncolonized', 'cleared', 'colonized'),
				labels = c('Uncolonized', '', 'Colonized'),
				values = c("#62f065", "#95b833", "#165f28"))
###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_5.jpg'),
  cowplot::plot_grid(cowplot::plot_grid(NULL, 
  						strep_day0_colonization_OTU_plot + 
  							labs(x = NULL, y = NULL, color = NULL) + 
							scale_y_log10(
								breaks=c(0.01, 0.1, 1, 10, 100),
								labels=c("0","0.1","1","10","100")) +
 							coord_cartesian(ylim = c(0.035,20)) + 
							coord_flip(),
  						rel_widths = c(1, 65)),
  					 cowplot::plot_grid(NULL, 
  					 	strep_day0_clearance_OTU_plot + 
  					 		labs(x = NULL, y = NULL, color = NULL) + 
							scale_y_log10(lim = c(0.035,20),
								breaks=c(0.01, 0.1, 1, 10, 100),
								labels=c("0","0.1","1","10","100")) +
							coord_flip(),
  						rel_widths = c(1, 65)),
					 strep_day10_clearance_OTU_plot + 
					 		labs(x = NULL, y = 'Relative Abundance (%)', color = NULL) + 
								scale_y_log10(lim = c(0.035,20),
								breaks=c(0.01, 0.1, 1, 10, 100),
								labels=c("0","0.1","1","10","100")) +
 							coord_flip(), 
                     ncol = 1,
                     rel_heights = c(8, 9, 26), 
                     labels = c('A', 'B', 'C')),
  height = 9, width = 5, unit = 'in')
###############################################################################
