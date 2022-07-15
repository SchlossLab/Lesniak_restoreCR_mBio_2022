###############################################################################
#
#	Common functions used across scripts
#
# need files:
# 	produced by code/clean_metadata.R
#		data/process/restore_metadata_clean.tsv
#	produced through running mothur scripts
#		data/mothur/sample.final.0.03.subsample.shared
#		data/mothur/sample.final.thetayc.0.03.lt.ave.nmds.axes
#		data/mothur/sample.final.groups.ave-std.summary
#		data/mothur/sample.final.thetayc.0.03.lt.ave.dist
#	produced by code/tidy_taxonomy.R
#		data/process/restore_taxonomy_clean.tsv
#
# Nick Lesniak 2022-01-12
###############################################################################
# setup environment
###############################################################################
library(tidyverse)
library(cowplot)
library(here)

###############################################################################
# functions to load data
###############################################################################
read_metadata <- function(){
  read_tsv(here('data/process/restore_metadata_clean.tsv'),
           col_types = c('ccdddcDcccdDccc'))  %>% 
    select(-file, -group) %>% 
    filter(!mouse_id %in% c('1A', '1B', '2A', '2B')) # remove mice from testing 2 day delay for clindamycin challenge
}

read_shared <- function(){
  read_tsv(here('data/mothur/sample.final.0.03.subsample.shared'),
           col_types = cols(
             .default = col_double(),
             Group = col_character()
           ))
}

read_nmds <- function(){
  read_tsv(here('data/mothur/sample.final.thetayc.0.03.lt.ave.nmds.axes'),
           col_types = c('cdd'))
}

read_taxonomy <- function(){
  read_tsv(here('data/process/restore_taxonomy_clean.tsv'),
           col_types = cols(.default = col_character()))
}

read_alpha <-  function(){
  read_tsv(here('data/mothur/sample.final.groups.ave-std.summary'),
           col_types = c('dccdddddd'))
}

read_beta <- function(){
	source(here('code/read_dist.R'))
	read_dist(here('data/mothur/sample.final.thetayc.0.03.lt.ave.dist'))
}

hline_df <- data.frame(measure = c('CFU (Log10)'),
                       yint = 100,
                       x_pos = 10,#min(data$day)
                       day = 10
                       )

#### color by antibiotic to match previous publication
abx_color <- tibble(antibiotic = c('Streptomycin', 'Cefoperazone', 'Clindamycin'),
                    color = c('#D37A1F', '#3A9CBC', '#A40019'))