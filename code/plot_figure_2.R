###############################################################################
#
# Figure 2
#	main idea - test Fecal community in different antibiotic induced susceptibility
# 	plot - cfu w/ and w/o Fecal community
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
lod_df <- data.frame(day = 1, cfu = 100, treatment = 'PBS')

plot_data <- metadata %>% 
  filter(day %in% c(1:10),
         treatment %in% c('PBS', '-1')) %>% 
	mutate(cfu = case_when(cfu == 0 ~ 60, T ~ cfu),  # shift 0 counts to just below limit of detection line
	       antibiotic = factor(antibiotic, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin')),
	       treatment = ifelse(treatment == '-1', 'FCT', treatment))

###############################################################################
#  analyze data
###############################################################################

###############################################################################
#  plot data
###############################################################################
cfu_plot <- plot_data %>% 
  ggplot(aes(x = day, y = cfu, color = antibiotic, shape = treatment)) + 
	geom_hline(yintercept = 101, linetype = 'dashed', 
		size = 1, color = 'white') + 
	geom_label(data = lod_df, label = "LOD", fill = 'white', color = 'white', size = 2) + 
	geom_text(data = lod_df, label = "LOD", color = 'grey', size = 3) + 
    geom_point(position = position_jitterdodge()) + 
	scale_y_log10(breaks = c(10^2, 10^4, 10^6, 10^8),
  				  labels = c('10^2', '10^4', '10^6', '10^8')) + # scale y axis log10 and label 10^x
    labs(x = 'Days post infection', 
    	 y = expression(italic('C. difficile')~' CFU/g Feces'), shape = 'Treatment') + 
    scale_color_manual(values = abx_color$color, breaks = abx_color$antibiotic) + 
    scale_shape_manual(values = c(1, 16), 
    				   breaks = c('PBS', 'FCT'), 
    				   labels = c('PBS', 'FCT')) + 
    scale_x_continuous(breaks = c(1:10)) + 
    facet_grid(antibiotic~factor(treatment, levels = c('PBS', 'FCT'))) + 
    theme_bw() + 
    theme(axis.text.y = ggtext::element_markdown(),
          strip.background = element_rect(fill = 'white', color = 'white'),
          panel.grid.minor = element_blank()) + 
    guides(color = 'none', shape = 'none')

###############################################################################
#  save plot
###############################################################################
ggsave(here('submission/Figure_2.tiff'),
  cfu_plot, 
  height = 3.5, width = 4, unit = 'in',
  compression = 'lzw')
###############################################################################
