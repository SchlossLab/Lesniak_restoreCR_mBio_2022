###
#
# Clean up issues in metadata file
#
# Dependencies: 
#	* data/raw/:
#		'restore_exp_04__05_06_19.xlsx', # clinda - 1Dfecal treatment 10,100,1000
#		'restore_exp_05_strep_24_06_19.xlsx', # strep - 2Drec2Dfecal treatment 10,100,1000
#		'restore_exp_06_strep_2019_07_08.xlsx') # strep - 2Drec2Dfecal treatment 10,100,1000
#		'restore_exp_07_strep_07222019.xlsx', # strep - 2Drec2Dfecal treatment 10,100,1000
#		'restore_exp_08_strep_08052019.xlsx', # strep - 2Drec2Dfecal treatment 10-10^6
#		'restore_exp_09_clinda_11122019.xlsx', # clinda - 1Dfecal treatment 10,100,1000
#		'restore_exp_09_clinda_pilot_1112019.xlsx', # clinda - 2Drec only one cage colonized @ 10^3
#		'restore_exp_10_cef_12032019.xlsx', # cef - 2Drec2Dfecal treatment 10-10^6
#		'restore_exp_11_strep_01132020.xlsx') # strep - 2Drec2Dfecal treatment 10-10^6
#
# Output: 
#	* data/process/restore_metadata_clean.tsv
#
###

#load dependencies 
library(readxl)
library(tidyverse)

#load metadata file
file_list_w_earmark <- c('restore_exp_04_clinda_05_06_19.xlsx', # clinda - 1D fecal treatment 10,100,1000
	'restore_exp_05_strep_24_06_19.xlsx', # strep - 2Drec2Dfecal treatment 10,100,1000
	'restore_exp_06_strep_2019_07_08.xlsx') # strep - 2Drec2Dfecal treatment 10,100,1000
file_list_wo_earmark <- c('restore_exp_07_strep_07222019.xlsx', # strep - 2Drec2Dfecal treatment 10,100,1000
	'restore_exp_08_strep_08052019.xlsx', # strep - 2Drec2Dfecal treatment 10-10^6
	'restore_exp_09_clinda_11122019.xlsx', # clinda - 1Dfecal treatment 10,100,1000
	'restore_exp_09_clinda_pilot_1112019.xlsx', # clinda - 2Drec only one cage colonized @ 10^3
	'restore_exp_10_cef_12032019.xlsx', # cef - 2Drec2Dfecal treatment 10-10^6
	'restore_exp_11_strep_01132020.xlsx') # strep - 2Drec2Dfecal treatment 10-10^6

meta_data <- bind_rows(
	map_dfr(file_list_w_earmark, function(x){
			read_xlsx(paste0('data/raw/', x), sheet = 'clean_cfu_df', na = 'NA',
				col_types = c('text', 'text', 'text', 'numeric', 'numeric', 'numeric', 
					'text', 'date', 'text')) %>% 
		mutate(file = x)
		}) %>% 
		select(-ear_mark),
	map_dfr(file_list_wo_earmark, function(x){
			read_xlsx(paste0('data/raw/', x), sheet = 'clean_cfu_df', na = c('NA', '#DIV/0!'),
				col_types = c('text', 'text', 'numeric', 'numeric', 'numeric', 
					'text', 'date', 'text')) %>% 
		mutate(file = x)})) %>% 
	mutate(sample_id = paste0(gsub('-', '', mouse_id), '_d', gsub('-', 'minus', day)),
			# sample_id = <treatment>_<experiment number>_<ear tag>_d<day>
		date = as.Date(date)) 

read_cage_dob_exp <- function(file_name, last_cell, mouse_id_col){
	read_xlsx(paste0('data/raw/', file_name), sheet = 'Inventory', 
			na = 'NA', range = paste0('A1:', last_cell), 
			col_types = c('text', 'numeric', 'numeric', 'skip', 'text', 'text', 
				'date', 'text', 'text', 'date', 'numeric', 'date', 'numeric', 'numeric', 'text')) %>% 
		mutate(dob = as.Date(DOB)) %>% 
		select(mouse_id = all_of(mouse_id_col), exp, dob, cage = `Cage Card #`)	
	}


cage_dob_exp_df <- bind_rows(
	read_cage_dob_exp(file_list_w_earmark[1], 'O11', 'unique_id'),
	read_cage_dob_exp(file_list_w_earmark[2], 'O9', 'unique_id'),
	read_cage_dob_exp(file_list_w_earmark[3], 'O9', 'unique_id'),
	read_cage_dob_exp(file_list_wo_earmark[1], 'O18', 'unique_id'),
	read_cage_dob_exp(file_list_wo_earmark[2], 'O16', 'unique_id'),
	read_cage_dob_exp(file_list_wo_earmark[3], 'O17', 'Mouse ID'),
	read_cage_dob_exp(file_list_wo_earmark[5], 'O16', 'Mouse ID'),
	read_cage_dob_exp(file_list_wo_earmark[6], 'O17', 'Mouse ID'),
	tibble(mouse_id = c('1A', '1B', '2A', '2B'),
		exp = 9, dob = as.Date('2019-08-12'), cage = c('1_9', '1_9', '2_9', '2_9')))


# create sample names for the fecal dilutions used as inocula
# inocula for two days, except exp 4 and 9
# exp 7 had two different full (10^-1) inocula, F1 is full from Exp5 and F2 is full from Exp6
inoculum_samples <- meta_data %>% 
	select(group, treatment, file) %>% 
	filter(!grepl('NT', group)) %>% 
	mutate(group = gsub('\\.0','', group)) %>% 
	distinct %>% 
	mutate(group = gsub('-', '', group),
		group = gsub('10\\^', '', group),
		group = gsub('F', '1', group),
		group = gsub('M', '2', group),
		group = gsub('L', '3', group),
		exp = str_extract(file, 'exp_\\d{2}'),
		# replace _EX5/6 with inoculum1/2 to match tube labels
		sample_id = case_when(grepl('EX5', group) ~ paste0('inoculum1_1to10e', gsub('_EX5', '', group), '_', exp),
			grepl('EX6', group) ~ paste0('inoculum2_1to10e', gsub('_EX6', '', group), '_', exp),
			T ~ paste0('inoculum_1to10e', group, '_', exp)),
		# add day inocula was used
		day = ifelse(exp %in% c('exp_09', 'exp_04'), -1, -2),
		group = 'inoculum')

inoculum_samples <- inoculum_samples %>% 
	filter(day == -2) %>% 
	mutate(day = -1,
		sample_id = paste0(gsub('exp_', 'exp', sample_id), '_d2')) %>% 
	bind_rows(inoculum_samples) %>% 
	select(-exp)

meta_data <- bind_rows(meta_data, inoculum_samples) %>% 
  mutate(treatment = case_when(group %in% c('F', 'F_EX5', 'F_EX6', '-1.0') ~ '-1',
                               group %in% c('M', '-2.0') ~ '-2',
                               group %in% c('L', '-3.0') ~ '-3',
                               group %in% c('10^-4', '-4.0') ~ '-4',
                               group %in% c('10^-5', '-5.0') ~ '-5',
                               group %in% c('10^-6', '-6.0') ~ '-6',
                               group %in% c('NTod', 'NT') ~ 'PBS',
                               T ~ group)) %>% 
  left_join(cage_dob_exp_df, by = 'mouse_id')

timepoint_df <- tibble(antibiotic = c(rep('Clindamycin', 3),
                                      rep('Cefoperazone', 5),
                                      rep('Streptomycin', 5)),
                       day = c(-2,-1, 0, rep(c(-9,-4,-2,-1, 0), 2)),
                       timepoint = c('initial', 'post_antibiotic', 'infection',
                                     rep(c('initial', 'stop_antibiotic', 
                                           'post_antibiotic', 'post_first_fmt', 'infection'), 
                                     	 2)))

# remove clinda mice ntod (testing no day added for FCT) 
# and 1A/B 2A/B (testing 2 days between clinda and Cdiff challenge)
meta_data <- meta_data %>% 
	filter(!grepl('ntod', mouse_id),
		!mouse_id %in% c('1A', '1B', '2A', '2B'))

final_cfu <- meta_data %>% 
	filter(!is.na(mouse_id),
		day <= 10,
		!is.na(cfu)) %>% 
	group_by(mouse_id) %>% 
	mutate(endpoint = max(day)) %>% 
	filter(endpoint == day) %>% 
	select(mouse_id, final_cfu = cfu) %>% 
	ungroup

max_cfu <- meta_data %>% 
	filter(!is.na(mouse_id),
		day <= 10) %>% 
	group_by(mouse_id) %>% 
	summarise(max_cfu = max(cfu, na.rm = T))

meta_data <- meta_data %>% 
	left_join(final_cfu, by = c('mouse_id')) %>% 
	left_join(max_cfu, by = c('mouse_id')) %>% 
	group_by(mouse_id) %>% 
	mutate(outcome = case_when(is.na(mouse_id) ~ 'NA',
		   				    max_cfu < 1000 ~ 'uninfected',
		   					final_cfu < 2000 ~ 'cleared',
		   					T ~ 'colonized')) %>% 
	ungroup %>% 	
	left_join(timepoint_df, by = c('antibiotic', 'day')) %>% 
	mutate(timepoint = ifelse(day > 0, 'dpi', timepoint)) %>% 
	select(-max_cfu, -final_cfu)

write_tsv(meta_data, 'data/process/restore_metadata_clean.tsv')