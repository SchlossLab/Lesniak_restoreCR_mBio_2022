###############################################################################
#
# Figure 4 and S3
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
	filter(antibiotic %in% c('Cefoperazone', 'Streptomycin'), 
		day < 11) %>% 
	mutate(treatment = case_when(treatment == 'PBS' ~ 'PBS',
      treatment == '-1' ~ 'FCT',
			T ~ paste0('1:10^', -(as.numeric(treatment) + 1))),
		treatment = factor(treatment, levels = c('FCT', paste0('1:10^', (1:5)), 'PBS')))

taxonomy <- read_taxonomy()

shared <- read_shared() %>% 
  select(-label, -numOtus) %>% 
  pivot_longer(-Group, names_to = 'otu', values_to = 'counts') %>% 
  mutate(relative_abundance = counts/2480 * 100)

beta_df <- read_beta()

alpha_df <- read_alpha() %>% 
	filter(method == 'ave')

###############################################################################
#  setup data
###############################################################################
min_rel_abund <- 100 * 1/2480
lod_df <- data.frame(x = 0.75, y = min_rel_abund)
dilution_colors <- c("#58b5e1", "#19477d", "#8b6fed", "#edb1ff", "#432ab7", "#ec4dd8", "#7c3b6b")

alpha_data <- metadata %>% 
	filter(#antibiotic == 'Streptomycin' |
		#(antibiotic == 'Cefoperazone' & treatment %in% c(-1, -2)),
		day %in% c(-9, 0, 10)) %>% 
	left_join(select(alpha_df, group, sobs, invsimpson), 
		by = c('sample_id' = 'group')) 

beta_data <- beta_df %>% 
	right_join(metadata %>% 
			filter(day %in% c(-9, 0, 10)) %>% 
			select(sample_id, mouse_id_r = mouse_id, treatment, antibiotic,
				day_r = day), 
		by = c('rows' = 'sample_id')) %>% 
	right_join(metadata %>% 
			filter(day %in% c(-9, 0, 10)) %>% 
			select(sample_id, mouse_id_c = mouse_id, day_c = day), 
		by = c('columns' = 'sample_id')) %>% 
	filter(mouse_id_r == mouse_id_c,
		day_c == -9 | day_r == -9) %>% 
	mutate(day = ifelse(day_r == -9, day_c, day_r))

###############################################################################
#  analyze data
###############################################################################

###############################################################################
#  plot data
###############################################################################
plot_diversity <- function(input_antibiotic){
	sobs_plot <- alpha_data %>% 
		filter(antibiotic == input_antibiotic) %>% 
		ggplot(aes(x = as.factor(day), y  = sobs, color = treatment)) + 
			scale_color_manual(values = dilution_colors) + 
			theme_bw() + 
			ylim(c(0, 150)) + 
			theme(axis.title.x = ggtext::element_markdown(),
				legend.position = 'none') + 
			labs(x = NULL,
				y = expression(~S[obs]), color = NULL)

	invsimpson_plot <- alpha_data %>% 
		filter(antibiotic == input_antibiotic) %>% 
		ggplot(aes(x = as.factor(day), y  = invsimpson, color = treatment)) + 
			theme_bw() + 
			scale_color_manual(values = dilution_colors) + 
			#scale_shape_manual(values = c(1:6,8)) + 
			ylim(c(0, 25)) + 
			theme(axis.title.x = ggtext::element_markdown(),
				legend.position = 'none') + 
			labs(x = NULL,
				y = 'Inverse Simpson')

	beta_plot <- beta_data %>% 
		filter(antibiotic == input_antibiotic) %>% 
		ggplot(aes(x = as.factor(day), y = distances, color = treatment)) + 
			ylim(c(0,1)) + 
			scale_color_manual(values = dilution_colors) + 
			theme_bw() +
			theme(axis.title.x = ggtext::element_markdown(),
				legend.text = ggtext::element_markdown(),
				legend.position = 'bottom') + 
			labs(x = 'Day', 
				y = expression(atop(~Theta[YC] ~ distance, "to day -9")), 
				color = NULL)

	if(input_antibiotic == 'Cefoperazone'){
		return(list(sobs_plot = sobs_plot + 
					geom_jitter(position=position_jitterdodge(0.1)), 
				invsimpson_plot = invsimpson_plot + 
					geom_jitter(position=position_jitterdodge(0.1)), 
				beta_plot = beta_plot + 
					geom_jitter(position=position_jitterdodge(0.1))))
	} else if(input_antibiotic == 'Streptomycin'){
		return(list(sobs_plot = sobs_plot + 
					stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
						position=position_dodge(0.5)), 
				invsimpson_plot = invsimpson_plot + 
					stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
						position=position_dodge(0.5)), 
				beta_plot = beta_plot + 
					stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
						position=position_dodge(0.5)) + 
					guides(color = guide_legend(override.aes = list(shape = 16, linetype = 0)))))

	}

}

cef_plots <- plot_diversity('Cefoperazone')
strep_plots <- plot_diversity('Streptomycin')

annotation_df <- data.frame(
	x1=c(0.85, 1.85), 
	x2=c(1.85, 2.15), 
	xnote = c(1.3, 2),
	y1 = c(0.95, 0.9),
	ynote = c(0.96, 0.91),
	annotations='*', outcome = 'NA')

beta_diff_plot <- beta_data %>% 
    left_join(select(metadata, sample_id, outcome), 
        by = c('rows' = 'sample_id')) %>%        
	filter(antibiotic == 'Streptomycin',
        outcome %in% c('colonized', 'cleared')) %>% 
	ggplot(aes(x = factor(day), y = distances, color = outcome)) + 
		ylim(c(0,1)) + 
		scale_color_manual(
			breaks = c('cleared', 'colonized'),
			labels = c('Cleared colonization', 'Remain colonized'),
			values = c("#c2dcb8", "#096013")) + 
		theme_bw() +
		theme(axis.title.x = ggtext::element_markdown(),
			legend.text = ggtext::element_markdown(),
			legend.position = 'bottom',
			strip.background = element_rect(color = 'white', fill = 'white'),
			panel.grid = element_blank()) + 
		labs(x = 'Day', 
			y = expression(atop(~Theta[YC] ~ distance, "to day -9")), 
			color = NULL) + 
        stat_summary(fun.data = 'median_hilow', fun.args = (conf.int=0.5),
         position=position_dodge(0.5)) + 
        guides(color = guide_legend(override.aes = list(shape = 16, linetype = 0))) + 
		geom_text(data = annotation_df, aes(x = xnote, y = ynote, label = annotations), 
			color = 'black') +
		geom_segment(data = annotation_df, aes(x = x1, xend = x2, y = y1, yend = y1), 
			color = 'black', size = 0.25)

###############################################################################
#  test differences
###############################################################################
#temp_df <- beta_data %>% 
#	left_join(select(metadata, sample_id, outcome), 
#		by = c('rows' = 'sample_id')) %>%        
#	filter(antibiotic == 'Streptomycin',
#		outcome %in% c('colonized', 'cleared')) %>% 
#	select(day, outcome, distances) %>% 
#	pivot_wider(names_from = 'outcome', values_from = 'distances')
# multiple comparison correction w/bonferroni, since few comparisons
#p.adjust(c(wilcox.test(temp_df$cleared[[1]], temp_df$cleared[[2]])$p.value,
#	wilcox.test(temp_df$colonized[[1]], temp_df$colonized[[2]])$p.value,
#	wilcox.test(temp_df$cleared[[1]], temp_df$colonized[[1]])$p.value,
#	wilcox.test(temp_df$cleared[[2]], temp_df$colonized[[2]])$p.value))
# colonized 0 vs 10 - P = 0.092
# colonized 10 vs cleared 10 - P =  0.006
# cleared 0 vs 10 - P = 2.3e-7
# colonized 0 vs cleared 0 - P =  0.092

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_4.jpg'),
	cowplot::plot_grid(
		strep_plots$sobs_plot,
		cowplot::plot_grid(NULL, strep_plots$invsimpson_plot,
			nrow = 1, rel_widths = c(.26, 10)),
		cowplot::plot_grid(NULL, strep_plots$beta_plot + 
				theme(legend.position = 'none'),
			nrow = 1, rel_widths = c(2.2, 10)),
		cowplot::plot_grid(get_legend(strep_plots$beta_plot)),
		ncol = 1,
		rel_heights = c(6.4, 6.4, 7, 1.4),
		labels = c('A', 'B', 'C', '')),
	height = 7, width = 3.75, unit = 'in')

ggsave(here('results/figures/Figure_S3.jpg'),
	cowplot::plot_grid(
		cef_plots$sobs_plot,
		cowplot::plot_grid(NULL, cef_plots$invsimpson_plot,
			nrow = 1, rel_widths = c(.26, 10)),
		cowplot::plot_grid(NULL, cef_plots$beta_plot + 
				theme(legend.position = 'none'), 
			nrow = 1, rel_widths = c(2.2, 10)),
		cowplot::plot_grid(get_legend(cef_plots$beta_plot)),
		ncol = 1,
		rel_heights = c(6.4, 6.4, 7, 1.4),
		labels = c('A', 'B', 'C', '')),
	height = 7, width = 3.75, unit = 'in')

ggsave(here('results/figures/Figure_S4.jpg'),
	beta_diff_plot,
	height = 3, width = 4, unit = 'in')
###############################################################################
