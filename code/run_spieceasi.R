##############
#
# create network of cdiff interactions
#	What interactions associate with clearance of C. difficile colonization?
# 
# Nick Lesniak 01-28-2022
#
#  need files:
#	data/process/
#	data/mothur/
#	data/process/
#
###############################################################################
# setup environment
###############################################################################
source(here::here('code/functions.R'))
library(SpiecEasi)
library(igraph)
library(network)
library(sna)
library(intergraph)
library(scales)

seed <- 18

# arguments for spiec easi, 
se_pargs <- list(rep.num=99, seed=seed, ncores=4)
# function to set offset angle to arrange node labels outside circle
#  kjhealy/polar-labels.r https://gist.github.com/kjhealy/834774
#radian.rescale <- function(x, start=0, direction=1) {
#  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
#  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
#}

###############################################################################
# load data
###############################################################################
metadata <- read_metadata() %>% 
	filter(antibiotic == 'Streptomycin',
		day %in% 1:5)

taxonomy <- read_taxonomy()

shared <- read_shared() %>% 
	select(-label, -numOtus) %>% 
	filter(Group %in% metadata$sample_id)

###############################################################################
#  setup data
###############################################################################
network_labels <- taxonomy %>% 
	mutate(otu_number = gsub('OTU ', '', otu_label),
		tax_otu_label = gsub('_unclassified', '', tax_otu_label), # remove unclassified note
		tax_otu_label = paste0('*', tax_otu_label), # italicize
		tax_otu_label = gsub(' \\(', '*<br/>(', tax_otu_label), # add line breaks between tax and otu
		tax_otu_label = gsub('_', ' ', tax_otu_label), #convert underscores to spaces
		tax_otu_label = gsub(' ([XI1V]+)\\*<', '\\* \\1<', tax_otu_label)#, # unitalicize roman numerals
		#tax_otu_label = case_when(grepl('OTU 10)|OTU 45)', tax_otu_label) ~ # bold otus in multiple networks
		#	paste0('**', tax_otu_label, '**'),
		#	T ~ tax_otu_label)
		) %>% 
	select(OTU, tax_otu_label) %>% 
	bind_rows(data.frame(OTU = c('Cdiff', 'Resistant', 'Cleared', 'Persistent'),
		tax_otu_label = c('*C. difficile*', 'Resistant', 'Cleared', 'Persistent')))

se_df <- metadata %>% 
	select(sample_id, cfu, outcome) %>% 
	mutate(cfu = ifelse(cfu == 0, log10(0.1), log10(cfu)),
		cfu = round(2480 * ((cfu - min(cfu, na.rm = T)) / 
				(max(cfu, na.rm = T) - min(cfu, na.rm = T))), 0)) %>% 
	filter(!is.na(cfu)) %>% 
	rename(Cdiff = cfu) %>% 
	inner_join(shared, by = c('sample_id' = 'Group')) %>% 
	select(-sample_id)

# select only OTUs present in 10% of samples
otus_present <- se_df %>% 
	summarise_all(function(x){
			sum(x >= 1) >= .1 * nrow(shared)
		}) %>% 
	gather(OTU, present) %>% 
	filter(present == T) %>% 
	pull(OTU)

se_df <- se_df %>% 
	select(all_of(otus_present))

SpiecEasi_data <- se_df %>% 
	select(-outcome) %>% 
	as.matrix

feature_names <- network_labels %>% 
	right_join(data.frame(order = 1:ncol(SpiecEasi_data),
		OTU = colnames(SpiecEasi_data)), by = c('OTU')) %>% 
	arrange(order) %>% 
	pull(tax_otu_label)

###############################################################################
#  analyze data
###############################################################################

# SPIEC-EASI: data transformation, sparse inverse covariance estimation and model selection
se_model <- spiec.easi(SpiecEasi_data, method = 'mb', lambda.min.ratio = 1e-3, nlambda = 500,
	sel.criterion = 'bstars', pulsar.select = TRUE, pulsar.params = se_pargs)

# set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(SpiecEasi_data, 1))+3
names(vsize) <- feature_names
# determine edge weights
se_beta <- symBeta(getOptBeta(se_model), mode='maxabs')
se_edges <- Matrix::summary(se_beta)
# network degree distributions
se_network <- adj2igraph(getRefit(se_model))
se_dd <- degree.distribution(se_network)
# determine network stability (closest to 0.05 is best, increase nlambda if not close)
se_stability <- getStability(se_model)
if(se_stability < 0.045){ stop(paste0('Stability low (', se_stability, 
	'), increase nlambda or decrease lambda.min.ratio arguments'))}

## setup model output for graphing network
se_interaction_matrix <- getRefit(se_model)
# name matrix positions
colnames(se_interaction_matrix) <- rownames(se_interaction_matrix) <- feature_names

get_network <- function(focal_feature){
	# subset network to only OTUs directly interecting with C. difficile
	first_order_otus <- c(names(which(se_interaction_matrix[,focal_feature] > 0)), focal_feature)
	focal_interactions <- se_interaction_matrix[first_order_otus, first_order_otus]
	#labels <- c('Cdiff', 
	#	network_labels[as.numeric(tail(head(colnames(focal_interactions), -3), -1))], 
	#	"Cleared", "Persistent", "Resistant")
	#colnames(focal_interactions) <- rownames(focal_interactions) <- labels

	# add edge weights
	# names matrix positions
	colnames(se_beta) <- rownames(se_beta) <- feature_names
	# subset network to OTUs interacting with C difficile
	wt_focal_interactions <- se_beta[first_order_otus, first_order_otus]
	# filter interactions with wieght near 0
	temp_matrix <- as.matrix(wt_focal_interactions)
	temp_matrix[-0.1 < temp_matrix & temp_matrix < 0.1] <- 0
	wt_focal_interactions <- as(temp_matrix, "dgCMatrix")
	#colnames(wt_focal_interactions) <- rownames(wt_focal_interactions) <- labels
	# create igraph object with edge weights
	wt_first_order_network <- adj2igraph(wt_focal_interactions, 
		vertex.attr = list(name = colnames(wt_focal_interactions)))

	# filter otus in matrix based on edge weights
	subset_nodes <- unique(as.vector(as_edgelist(wt_first_order_network)))
	focal_interactions <- se_interaction_matrix[subset_nodes, subset_nodes]
	#labels <- c('*C. difficile*', 
	#	network_labels[as.numeric(tail(head(colnames(focal_interactions), -3), -1))], 
	#	"Cleared", "Persistent", "Resistant")
	#colnames(focal_interactions) <- rownames(focal_interactions) <- labels
	wt_focal_interactions <- se_beta[subset_nodes, subset_nodes]
	wt_first_order_network <- adj2igraph(wt_focal_interactions, 
		vertex.attr = list(name = colnames(wt_focal_interactions)))

	# setup network attributes to create igraph network graph for output
	vsize_order <- c()
	for( i in 1:length(subset_nodes)){
	  vsize_order <- c(vsize_order, vsize[names(vsize) == subset_nodes[i]])
	}
	focal_vsize <- vsize_order
	#vsize['Cdiff'] <- 7.0711
	edge_wt <- abs(E(wt_first_order_network)$weight)
	edge_direction <- ifelse(E(wt_first_order_network)$weight < 0, 'red', 'blue')
	#lab.locs <- radian.rescale(x=1:length(first_order_otus), direction=-1, start=0)
	focal_network <- adj2igraph(focal_interactions, 
			vertex.attr = list(name = colnames(focal_interactions), 
				size = round(focal_vsize^2/10, 2), 
				#color = abx_col, 
				label.color='black', label.cex = 0.7, label.dist = 2#, 
				#label.degree = lab.locs
				),
			edge.attr = list(width = round(edge_wt*10, 2), color = edge_direction))
	return(focal_network)
}

network_data <- get_network('*C. difficile*')

###############################################################################
#  save data
###############################################################################
saveRDS(network_data, here::here('data/process/cdiff_network.rds'))
###############################################################################