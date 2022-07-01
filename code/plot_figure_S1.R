###############################################################################
#
# Figure S1
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
metadata <- read_metadata()

###############################################################################
#  setup data
###############################################################################
metadata <- metadata %>% 
    mutate(cfu = as.numeric(cfu) + 50) %>% 
  	filter(antibiotic == 'Streptomycin',
  	       day %in% -1:10,
  	       !is.na(cfu)) %>% 
	mutate(treatment = case_when(treatment == 'PBS' ~ 'PBS',
      treatment == '-1' ~ 'FCT',
			T ~ paste0('1:10^', -(as.numeric(treatment) + 1))),
		treatment = factor(treatment, levels = c('FCT', paste0('1:10^', (1:5)), 'PBS')))

###############################################################################
#  plot data
###############################################################################
individual_strep_plot <- metadata %>% 
    ggplot(aes(x = day, y = cfu)) + 
		scale_y_log10(breaks = c(10^2, 10^4, 10^6, 10^8),
	  		labels = c('10^2', '10^4', '10^6', '10^8')) + # scale y axis log10 and label 10^x
  		geom_hline(data = hline_df, aes(yintercept = yint), 
  			size = 1.1, linetype = 'dashed', color = 'white') + 
  		geom_label(data = hline_df, aes(x = (-1 + .5), y = yint, label = "LOD"), 
  			fill = "white", color = 'gray', label.size = NA, inherit.aes = FALSE) + 
  		geom_line(aes(group = mouse_id), color = '#D37A1F', alpha = 0.3) + 
  		labs(x = 'Day', y = expression(italic('C. difficile')~' CFU/g Feces')) + 
  		scale_x_continuous(breaks=-1:10) + 
  		theme_bw() + 
		theme(axis.text.y = ggtext::element_markdown(),
      		legend.position = 'none',
            panel.grid.minor.x = element_blank(),
            strip.background =element_rect(fill='#D37A1F'),
            strip.text.x = element_text(colour = 'white'),
            strip.text.y = ggtext::element_markdown(colour = 'white')) + 
      facet_grid(treatment~.) 

###############################################################################
#  save plot
###############################################################################
ggsave(here('submission/Figure_S1.tiff'),
	individual_strep_plot,
	height = 6, width = 4.5, unit = 'in',
	compression = 'lzw')
###############################################################################
