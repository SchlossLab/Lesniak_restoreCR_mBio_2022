################################################################################
#
# make_files_file.R
#
# This script will build the files file that will be used in the make.contigs
# command from within mothur.
#
# Dependencies...
# * fastq files stored in data/mothur
#
# Output...
# * data/mothur/restore_cr.files
#
################################################################################

library(tidyverse)

# get the list of fastq file names
fastqs <- list.files("data/mothur/", pattern="*.fastq.gz$")

# create a sample column using the samples cage, mouse ear tag and day to match meta
files_df <- tibble(seq_id = fastqs) %>% 
	mutate(seq_R = case_when(grepl('R1', seq_id) ~ 'R1', 
			grepl('R2', seq_id) ~ 'R2',
			T ~ 'NA')) %>% # create column for fastq pairs
	mutate(group = gsub('\\_S\\d.*', '', seq_id)) %>% # remove sequencing filename suffix
	pivot_wider(names_from = 'seq_R', values_from = 'seq_id') %>% # spread fastqs into R1/R2 columns
	select(group, R1, R2)

write_tsv(files_df, "data/mothur/restore_cr.files", col_names = F)
