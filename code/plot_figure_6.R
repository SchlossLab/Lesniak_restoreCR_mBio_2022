###############################################################################
#
# Figure 6
#	main idea - 
# 	plot - 
#
# Nick Lesniak 2022-01-31
###############################################################################
# setup environment
###############################################################################

library(GGally)
# inorder to have mixed text formatting in network labels
# changed geom_text() to geom_richtext() - lines 904 and 1068 
# used trace(ggnet2, edit='nano') to insert changes into ggnet2 function
# also edited R function saved in code/R/functions/ggnet2.R
library(ggplot2)
library(ggtext)
library(network)
library(igraph)
library(sna)
library(intergraph)
source('code/ggnet2.R')
#library(scales)

###############################################################################
# load data
###############################################################################
network_data <- readRDS('data/process/cdiff_network.rds')

###############################################################################
#  setup data
###############################################################################

###############################################################################
#  analyze data
###############################################################################

###############################################################################
#  plot data
###############################################################################
graph_network <- function(n, mode){
	ggnet2(n, mode = mode,
			size = 0.2,
			label = T, 
			#vjust = 1.3, 
			label.size = 7/.pt,
			#edge.label = 'width', 
			#edge.size = 'width', 
			edge.color = 'color', 
			edge.label.size = 7/.pt,
			layout.exp = 0.2) +
		ylim(-0.1, 1.1) + xlim(-0.1, 1.1) + 
		theme(axis.title = element_blank(), 
			axis.text =element_blank(),
				axis.ticks = element_blank(),
				legend.position ='bottom') + 
		guides(size = guide_legend(override.aes = list(color = '#A40019'))) + 
		#	 scale_size_manual(name = '',
		#                  breaks = c(5, 7.02, 2.67, 3.43, 3.38, 3.28), 
		#                  labels = c('5', '8.58', '0.07', '0.17', '0.13', '0.125'), 
		#                  values = c(5, 7.02, 2.67, 3.43, 3.38, 3.28)) +
		guides(size = F)
}
# legible layouts for cdiff
# fruchtermanreingold
# seeds 38 49 53 59 77 80 93
# target
# seeds 6 34 76 85 91 99
set.seed(6)
cdiff_graph <- graph_network(network_data, 'target')

###############################################################################
#  save plot
###############################################################################
ggsave(here::here('results/figures/Figure_6.jpg'),
	cdiff_graph,
	height = 4, width = 5, unit = 'in')
###############################################################################