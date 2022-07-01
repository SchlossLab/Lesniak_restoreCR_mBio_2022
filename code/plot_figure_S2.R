###############################################################################
#
# Figure S2
#	main idea - compare structure of dilutions
# 	plot - alpha and beta for dilutions
#		qpcr of dilutions
#
# Nick Lesniak 2022-01-12
###############################################################################
# setup environment
###############################################################################
source(here::here('code/functions.R'))

###############################################################################
# load data
###############################################################################
metadata <- read_metadata()

taxonomy <- read_taxonomy()

shared <- read_shared() %>% 
  select(-label, -numOtus) %>% 
  pivot_longer(-Group, names_to = 'OTU', values_to = 'counts') %>% 
  mutate(relative_abundance = counts/2480 * 100)

beta_df <- read_beta()

alpha_df <- read_alpha()

qpcr_data <- read_tsv(here::here('data/raw/inocula_16S_qpcr_022521_Abs_Quant.txt'),
		col_types = c('ccdddccccd'))

###############################################################################
#  setup data
###############################################################################
metadata <- metadata %>% 
	mutate(exp = as.numeric(ifelse(treatment == 'inoculum',
						gsub('.*exp_*(\\d*)(_d2)*', '\\1', sample_id),
						exp)),
		   dilution = ifelse(treatment == 'inoculum',
		   				gsub('in.*1to10e(\\d*)_exp.*', '\\1', sample_id),
		   				treatment))

dilution_samples <- metadata %>% 
	filter(treatment == 'inoculum') %>% 
	select(sample_id, day, dilution) %>% 
	separate(sample_id, sep = '_', remove = F,
			 into = c('sample_type', 'dilution_factor', 'experiment', 'dose')) %>% 
	mutate(dilution_factor = gsub('1to10e', '-', dilution_factor),
		   experiment = as.numeric(gsub('exp', '', experiment)),
		   dose = as.numeric(gsub('d', '', dose)))

# select only untreated mice or inocula samples
beta_samples <- metadata %>% 
	filter(treatment == 'inoculum' | timepoint == 'initial') %>% 
	pull(sample_id)

# for each experiement compare to initial 
beta_plot_data <- beta_df %>% 
	filter(rows %in% beta_samples & columns %in% beta_samples) %>% 
	left_join(select(metadata, sample_id, 
					 dilution_r = dilution, exp_r = exp), 
			  by = c('rows' = 'sample_id')) %>% 
	left_join(select(metadata, sample_id, 
			         dilution_c = dilution, exp_c = exp), 
			  by = c('columns' = 'sample_id')) %>%  
	mutate(dilution_r = ifelse(grepl('minus', rows), '0', dilution_r),
		   dilution_c = ifelse(grepl('minus', columns), '0', dilution_c)) %>% 
	filter(dilution_r == '0' | dilution_c == '0') %>% 
	mutate(comparison = case_when(dilution_r == '0' & dilution_c == '0' ~ '0',
		   						  dilution_r == '0' ~ dilution_c, 
		   						  T ~ dilution_r),
		   comparison = -as.numeric(comparison))
# zymo mock used for positive control
# ~ 1.4 x 10^10 cells/ml
# qpcr reaction used 

dilution_shared <- dilution_samples %>% 
	select(sample_id, day, dilution) %>% 
	inner_join(shared, by = c('sample_id' = 'Group')) %>% 
	left_join(taxonomy, by = c('OTU')) %>% 
	rename(taxonomic_group = Genus) %>% 
	group_by(sample_id, taxonomic_group, day, dilution) %>% 
	summarise(relative_abundance = sum(relative_abundance)) %>% 
	mutate(present = relative_abundance > 0.1) %>% 
	group_by(taxonomic_group) %>% 
	mutate(n_present = sum(present),
		   median_abundance = median(relative_abundance),
		   taxonomic_group = gsub('_unclassified', '*', taxonomic_group)) %>% 
	filter(n_present > 4)

taxonomic_order <- dilution_shared %>% 
	select(taxonomic_group, median_abundance) %>% 
	distinct %>% 
	arrange(desc(median_abundance)) %>% 
	pull(taxonomic_group)

###############################################################################
#  analyze data
###############################################################################
dilution_correlation <- dilution_shared %>% 
	filter(dilution < 5) %>% 	
	mutate(dilution = 10^-as.numeric(dilution)) %>% 
	group_by(taxonomic_group) %>% 
	nest() %>% 
	mutate(spearman_p = map(.x = data, .f = ~ cor.test(.x$relative_abundance, .x$dilution, 
			method = 'spearman', exact = F)$p.value),
		spearman_rho = map(.x = data, .f = ~ cor.test(.x$relative_abundance, .x$dilution, 
			method = 'spearman', exact = F)$estimate)) %>% # compare cleared vs colonized
	unnest(c(spearman_p, spearman_rho)) %>% 
	mutate(pvalue = p.adjust(spearman_p, method = 'BH', n = nrow(.))) %>% # correct p values
	filter(pvalue < 0.05) %>% # select only those above 0.05 after pvalue correction
	unnest(data)

###############################################################################
#  plot data
###############################################################################
dilution_sobs_plot <- alpha_df %>% 
	filter(method == 'ave') %>% 
	inner_join(dilution_samples, by = c('group' = 'sample_id')) %>% 
	ggplot(aes(x = as.numeric(dilution_factor), y = sobs)) +
		geom_point(alpha = 0.3) + 
		ylim(c(0,170)) + 
		scale_x_continuous(breaks = c(-6, -5, -4, -3, -2, -1),
  						   labels = c('1:10^5', '1:10^4', '1:10^3', 
  						   			  '1:10^2', '1:10', 'FCT')) + 
		theme_bw() + 
		labs(x = 'Fecal community dilution', y = expression(~S[obs])) + 
		theme(axis.text.x = ggtext::element_markdown(),
          strip.background = element_rect(color = 'white', fill = 'white'))

dilution_invsimp_plot <- alpha_df %>% 
	filter(method == 'ave') %>% 
	inner_join(dilution_samples, by = c('group' = 'sample_id')) %>% 
	ggplot(aes(x = as.numeric(dilution_factor), y = invsimpson)) +
		geom_point(alpha = 0.3) + 
		ylim(c(0,35)) + 
		scale_x_continuous(breaks = c(-6, -5, -4, -3, -2, -1),
  						   labels = c('1:10^5', '1:10^4', '1:10^3', 
  						   			  '1:10^2', '1:10', 'FCT')) + 
		theme_bw() + 
		labs(x = 'Fecal community dilution', y = 'Inverse Simpson') + 
		theme(axis.text.x = ggtext::element_markdown(),
          strip.background = element_rect(color = 'white', fill = 'white'))

dilution_beta_plot <- beta_plot_data %>% 
	ggplot(aes(x = comparison, y = distances, group = comparison)) + 
	    stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
		scale_x_continuous(breaks = c(-6, -5, -4, -3, -2, -1, 0),
  					  labels = c('1:10^5', '1:10^4', '1:10^3', 
  						   		 '1:10^2', '1:10', 'FCT', 'Feces')) +  
		theme_bw() + 
		ylim(c(0,1)) + 
		labs(x = 'Fecal community dilution', y = expression(~Theta[YC])) + 
		theme(axis.text.x = ggtext::element_markdown())

dilution_qpcr_plot <- qpcr_data %>% 
	mutate(dilution = case_when(dilution == 'PBS' ~ '7',
			dilution == 'Positive\nControl' ~ '0',
			dilution == 'Negative\nControl' ~ '8',
			T ~ dilution),
		dilution = -as.numeric(dilution)) %>% 
	filter(dilution %in% -(1:7)) %>% 
	ggplot(aes(x = dilution, y = Cq)) + 
		stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5)) +
		scale_x_continuous(breaks = c(#-8, 
									  -7, -6, -5, -4, -3, -2, -1#, 0
									  ),
  					  labels = c(#'Negative<br />Control', 
  					  			 'PBS', '1:10^5', '1:10^4', '1:10^3', 
  						   		 '1:10^2', '1:10', 'FCT'#, 'Positive<br />Control'
  						   		 )) +  
		scale_y_continuous(trans = 'reverse') +
		theme_bw() + 
		theme(axis.text.x = ggtext::element_markdown()) + 
		labs(x = 'Fecal community dilution')

dilution_correlation_plot <- dilution_correlation %>%
	mutate(dilution = log10(dilution)) %>% 
	ggplot(aes(x = dilution, y = relative_abundance)) + 
		geom_point(alpha = 0.3) + 
		facet_wrap(.~factor(taxonomic_group, levels = taxonomic_order), 
			scale ='free_y', nrow = 1) + 
		theme_bw() + 
		guides(color = 'none') + 
		scale_x_continuous(breaks = c(-4, -3, -2, -1),
  						   labels = c('1:10^3', '1:10^2', '1:10', 'FCT')) + 
		labs(x = 'Fecal community dilution', 
			 y = 'Relative Abundance (%)') + 
		theme(axis.text.x = ggtext::element_markdown(),
		  strip.background = element_rect(color = 'white', fill = 'white'))
###############################################################################
#  save plot
###############################################################################
ggsave(here('submission/Figure_S2.tiff'),
	cowplot::plot_grid(
		cowplot::plot_grid(dilution_sobs_plot, dilution_invsimp_plot,
	                     dilution_beta_plot, dilution_qpcr_plot, 
	                     ncol = 2, 
	                     rel_heights = c(1, 1),
	                     labels = c('A', 'B', 'C', 'D')),
		cowplot::plot_grid(dilution_correlation_plot, labels = c('E')),
		ncol = 1, 
		rel_heights = c(3, 2)),
	height = 8, width = 8, unit = 'in',
	compression = 'lzw')
###############################################################################
