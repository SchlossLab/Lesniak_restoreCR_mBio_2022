###############################################################################
#
# Figure 3
#	main idea - test effect of fmt dilutions
# 	plot - cfu on 1 and 10 dpi for cef and strep when using diluted fmt
#		relative abundance distribution for otus correlated with dilution
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

###############################################################################
#  setup data
###############################################################################
lod_df <- data.frame(day = '1 dpi', cfu = 100, treatment = '10^-5')

metadata <- metadata %>% 
	mutate(exp = as.numeric(ifelse(treatment == 'inoculum',
						gsub('.*exp_*(\\d*)(_d2)*', '\\1', sample_id),
						exp)),
		   dilution = ifelse(treatment == 'inoculum',
		   				gsub('in.*1to10e(\\d*)_exp.*', '\\1', sample_id),
		   				treatment))

dilution_cfu_data <- metadata %>% 
  filter(day %in% c(1, 10),
		 antibiotic %in% c('Cefoperazone', 'Streptomycin')) %>% 
	mutate(cfu = case_when(cfu == 0 ~ 60, T ~ cfu),
		   day = paste(day, 'dpi')) # shift 0 counts to just below limit of detection line

###############################################################################
#  plot data
###############################################################################
dilution_cfu_plot <- dilution_cfu_data %>% 
	filter(treatment %in% paste0('-', 2:6)) %>%
	mutate(treatment = paste0('10^', as.numeric(treatment) + 1)) %>%  
	ggplot(aes(x = treatment, y = cfu, color = antibiotic)) + 
		geom_point(position = position_jitterdodge()) + 
			scale_y_log10(breaks = c(10^2, 10^4, 10^6, 10^8),
	  		labels = c('10^2', '10^4', '10^6', '10^8')) + # scale y axis log10 and label 10^x
			scale_x_discrete(breaks = c('10^-1', '10^-2', '10^-3', '10^-4', '10^-5'),
	  		labels = c('1:10', '1:10^2', '1:10^3', '1:10^4', '1:10^5')) + 
		labs(x = NULL, y = expression(italic('C. difficile')~' CFU/g Feces'), shape = 'Treatment') + 
		scale_color_manual(values = abx_color$color[1:2], breaks = abx_color$antibiotic[1:2]) + 
		facet_grid(antibiotic~day) + 
		  geom_hline(yintercept = 101, linetype = 'dashed', size = 0.25, color = 'grey') + 
			geom_label(data = lod_df, label = "LOD", fill = 'white', color = 'white', size = 2) + 
			geom_text(data = lod_df, label = "LOD", color = 'grey', size = 3) + 
		theme_bw() + 
		theme(axis.text.y = ggtext::element_markdown(),
			axis.text.x = ggtext::element_markdown(),
			strip.background = element_rect(color = 'white', fill = 'white'),
			panel.spacing.y = unit(2, "lines")) + 
		guides(color = 'none')

###############################################################################
#  save plot
###############################################################################
ggsave(here('results/figures/Figure_3.jpg'),
	cowplot::plot_grid(
		cowplot::plot_grid(NULL, NULL, 
					 ncol = 1, 
					 rel_heights = c(1, 1),
					 labels = c('A', 'B')),
		dilution_cfu_plot,
		nrow = 1, 
		rel_widths = c(1, 20)),
  height = 4, width = 5, unit = 'in')
###############################################################################
