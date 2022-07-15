
###############################################################################
#
# Run LEfSe
#
# need files:
# 	data/process/metadata_tidy.tsv 
#	data/mothur/sample.final.0.03.subsample.shared 
#	data/process/final.taxonomy.tidy.tsv
#
# output files:
#	data/process/lefse/strep_day0_clearance_OTU.0.03.lefse_summary
#	data/process/lefse/strep_day0_clearance_OTU.design
#	data/process/lefse/strep_day10_clearance_OTU.0.03.lefse_summary
#	data/process/lefse/strep_day10_clearance_OTU.design
#	data/process/lefse/strep_day0_colonization_OTU.0.03.lefse_summary
#	data/process/lefse/strep_day0_colonization_OTU.design
#	data/process/lefse/strep_antibiotic.0.03.lefse_summary
#	data/process/lefse/strep_antibiotic.design
#
###############################################################################
# setup environment
###############################################################################
source(here::here('code/functions.R'))

###############################################################################
# load data
###############################################################################
metadata <- read_metadata() %>% 
	filter(antibiotic %in% c('Streptomycin', 'Cefoperazone'))

taxonomy <- read_taxonomy()

shared <- read_shared()

###############################################################################
# create function to run lefse
###############################################################################
run_lefse <- function(sample_df_name, tax_level){
	i <- sample_df_name
	current_df <- get(i)

	cfu_column <- metadata %>% 
		select(Group = sample_id, cfu) %>% 
		mutate(cfu = ifelse(cfu == 0, log10(60), log10(cfu)),
			cfu = round(2480 * ((cfu - min(cfu, na.rm = T)) / 
											(max(cfu, na.rm = T) - min(cfu, na.rm = T))), 0)) %>% 
		rename(C_difficile_CFU = cfu)

	current_shared <- shared %>% 
		filter(Group %in% current_df$Group) %>% 
		pivot_longer(cols = -c('label', 'Group', 'numOtus'),
			names_to = "OTU", values_to = 'value') %>% 
		left_join(taxonomy, by = c('OTU')) %>% 
		group_by(label, Group, numOtus, .data[[tax_level]]) %>% 
		summarise(value = sum(value)) %>% 
		pivot_wider(names_from = .data[[tax_level]], values_from = value) %>% 
		ungroup %>% 
	# add cdifficile cfu
		left_join(cfu_column, by = c('Group'))

	# remove otus that either have 0 or only present in 5 or fewer samples
	present_otus <- current_shared %>% 
		select(-label, -Group, -numOtus) %>% 
		map_dbl(~ sum(. > 0, na.rm = T)) %>% 
		which(x = (. > round(nrow(current_df)*.1, 1))) %>% 
		names
	current_shared <- current_shared %>% 
		select(label, Group, numOtus, one_of(present_otus)) %>% 
		mutate(numOtus = length(present_otus))

	# write files to be used in mothur for lefse analysis
	write_tsv(file = paste0('data/process/lefse/', i, '_', tax_level, '.shared'), 
		x = filter(current_shared, Group %in% current_df$Group))
	write_tsv(file = paste0('data/process/lefse/', i, '_', tax_level, '.design'), 
		x = filter(current_df, Group %in% current_shared$Group))
	print(tax_level)
	# run lefse
	system(paste0('code/mothur/mothur "#set.dir(input=data/process/lefse, output=data/process/lefse);
		lefse(shared=', i, '_', tax_level, '.shared, design=', i, '_', tax_level, '.design)"'))
}


###############################################################################
# setup data for comparisons
###############################################################################
#Streptomycin  
#	difference between colonized and uncolonized for 10^-1, 10^-2, 10^-3  
#		difference in day 0 relative abundance  
strep_day0_colonization <- metadata %>% 
	filter(antibiotic == 'Streptomycin',
		treatment %in% -(1:3),
		day == 0) %>% 
	mutate(outcome = ifelse(outcome == 'uninfected', 'uncolonized', 'colonized')) %>% 
	select(Group = sample_id, outcome)
#	Difference between cleared and persistent for 10^-3, 10^-4, 10^-5  
#		difference in day 0 relative abundance  
strep_day0_clearance <- metadata %>% 
	filter(antibiotic == 'Streptomycin',
		treatment %in% -(3:5),
		day == 0,
		outcome %in% c('colonized', 'cleared')) %>% 
	select(Group = sample_id, outcome)
#		difference in endpoint relative abundance  
strep_day10_clearance <- metadata %>% 
	filter(antibiotic == 'Streptomycin',
		treatment %in% -(3:5),
		day == 10,
		outcome %in% c('colonized', 'cleared')) %>% 
	select(Group = sample_id, outcome)
#		difference with streptomycin treatment
strep_antibiotic <- metadata %>% 
	filter(antibiotic == 'Streptomycin',
		day %in% c(-9, -2)) %>% 
	select(Group = sample_id, timepoint)

# Clindamycin treatment
# no significant differences at either OTU/Genus on day 0
#clinda_endpoint <- metadata %>% 
#	filter(antibiotic == 'Clindamycin',
#		treatment %in% c('PBS', '-1'),
#		day > 4,
#		cfu < 1000,
#		!grepl('ntod', mouse_id)) %>% 
#	group_by(mouse_id, treatment) %>% 
#	mutate(min_day = min(day)) %>% 
#	filter(day == min_day) %>% 
#	ungroup %>% 
#	select(Group = sample_id, treatment)
# endpoint differences not comparable to antibiotic treatments

# Cefoperazone
#cef_day0_clearance <- metadata %>% 
#	filter(day == 0,
#		antibiotic == 'Cefoperazone',
#		treatment %in% c('PBS', '-1')) %>% 
#	select(Group = sample_id, treatment)
#cef_day10_clearance <- metadata %>% 
#	filter(day == 10,
#		antibiotic == 'Cefoperazone',
#		treatment %in% c('PBS', '-1')) %>% 
#	select(Group = sample_id, treatment)
# only two samples per group, only sig diff is CFU at day 10

###############################################################################
# load data
###############################################################################

run_lefse('strep_day0_clearance', 'OTU')
run_lefse('strep_day10_clearance', 'OTU')
run_lefse('strep_day0_colonization', 'OTU')
run_lefse('strep_antibiotic', 'OTU')

###############################################################################