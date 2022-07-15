################################################################################
#
#	Part 1: Get the reference files
#
#	Here we give instructions on how to get the necessary reference files that
#	are used throughout the rest of the analysis. These are used to calculate
#	error rates, generate alignments, and classify sequences. We will pull down
#	the mock community reference (zymo_mock.fasta), the silva reference alignment
#	(silva.bacteria.align), and the RDP training set data (trainset16_022016).
#	Finally, we use the zymo_mock.align to get the alignment coordinates for the
#	V4 data. These data will be stored in the data/references/ folder.
#
#	The targets in this part are all used as dependencies in other rules
#
################################################################################

#Location of files
REFS = data/references/
RAW = data/raw/
MOTHUR = data/mothur/
PROC = data/process/

# create all data and reference files
#get v4 region of the silva reference alignment, rdp training set data and hmp mock
$(REFS)silva.v4.align $(REFS)trainset16_022016.v4.fasta $(REFS)trainset16_022016.v4.tax $(REFS)zymo_mock.align : code/get_references.batch
	bash code/get_references.batch

references : $(REFS)silva.v4.align $(REFS)trainset16_022016.v4.fasta $(REFS)trainset16_022016.v4.tax $(REFS)zymo_mock.align

################################################################################
#
#	Part 2: Get fastq files
#
################################################################################

# get the fastq files
$(RAW)get_data : code/get_fastqs.sh $(MOTHUR)SRR_Acc_List.txt
	bash code/get_fastqs.sh\
	touch $(RAW)get_data

# build the files file. 
$(MOTHUR)restore_cr.files : code/make_files_file.R $(RAW)get_data
	Rscript code/make_files_file.R


################################################################################
#
#	Part 3: Run data through mothur 
#
################################################################################

SUB = 2480

# install mothur in the project code directory
code/mothur/mothur : code/install_mothur.sh
	sh code/code/install_mothur.sh

run_mothur=code/mothur/mothur

# here we go from the raw fastq files and the files file to generate a fasta,
# taxonomy, and count_table file that has had the chimeras removed as well as
# any non bacterial sequences
# then we get the sequencing error as seen in the mock community samples
# then we go from the good sequences and generate a shared file and a cons.taxonomy 
# file based on OTU data. Lastly, we rarefy the number of reads to $(SUB) sequences 
# per library for the alpha and beta diversity analyses and modeling and 
# for the shared file
$(MOTHUR)complete.sample.final.shared $(MOTHUR)sample.final.shared $(PROC)mock.sample.final.shared $(MOTHUR)final.taxonomy $(MOTHUR)sample.error.count : code/get_good_seqs.batch\
									code/get_error.batch\
									code/get_shared_otus.batch\
									$(REFS)silva.v4.align\
									$(REFS)trainset16_022016.v4.fasta\
									$(REFS)trainset16_022016.v4.tax\
									$(REFS)HMP_MOCK.v4.fasta
	$(run_mothur) code/get_good_seqs.batch
	$(run_mothur) code/get_error.batch
	$(run_mothur) code/get_shared_otus.batch
	#rename and split shared files
	mv $(MOTHUR)restore_cr.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared $(MOTHUR)complete.sample.final.shared
	$(run_mothur) "#remove.groups(shared=$(MOTHUR)complete.sample.final.shared, groups=mock-mock3-mock4-mock5-mock6-mock8-mock9-mock10-mock11-mock12-mock13-mock15-mock16-mock32-mock53)"
	mv $(MOTHUR)complete.sample.final.0.03.pick.shared $(MOTHUR)sample.final.shared
	$(run_mothur) "#set.dir(input=data/mothur, output=data/mothur);
		summary.single(shared=sample.final.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=$(SUB));
		dist.shared(shared=sample.final.shared, calc=thetayc, subsample=$(SUB));
		sub.sample(shared=sample.final.shared, size=$(SUB));
		get.groups(shared=complete.sample.final.shared, groups=mock-mock3-mock4-mock5-mock6-mock8-mock9-mock10-mock11-mock12-mock13-mock14-mock15-mock16-mock32-mock53)"
	mv $(MOTHUR)complete.sample.final.0.03.pick.shared $(MOTHUR)mock.sample.final.shared
	mv $(MOTHUR)restore_cr.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy $(MOTHUR)final.taxonomy
	mv $(MOTHUR)restore_cr.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.error.count $(MOTHUR)sample.error.count
################################################################################
#
#	Part 4: Write the paper
#
################################################################################


################################################################################
# Process raw data for analysis
################################################################################

# Clean metadata for this project
$(PROC)restore_metadata_clean.tsv : $(RAW)restore_exp_04__05_06_19.xlsx\
									$(RAW)restore_exp_05_strep_24_06_19.xlsx\
									$(RAW)restore_exp_06_strep_2019_07_08.xlsx\
									$(RAW)restore_exp_07_strep_07222019.xlsx\
									$(RAW)restore_exp_08_strep_08052019.xlsx\
									$(RAW)restore_exp_09_clinda_11122019.xlsx\
									$(RAW)restore_exp_09_clinda_pilot_1112019.xlsx\
									$(RAW)restore_exp_10_cef_12032019.xlsx\
									$(RAW)restore_exp_11_strep_01132020.xlsx\
									code/clean_metadata.R
	Rscript code/clean_metadata.R

# Convert taxonomy file into a dataframe with OTU labels
$(PROC)final.taxonomy.tidy.tsv : $(MOTHUR)final.taxonomy\
									code/tidy_taxonomy.R
	Rscript code/tidy_taxonomy.R


################################################################################
#
# Run LEfSE and RF
#
################################################################################

# Run LEfSe analysis
$(PROC)lefse/strep_day0_clearance.0.03.lefse_summary $(PROC)lefse/strep_day0_clearance.design $(PROC)lefse/strep_day10_clearance.0.03.lefse_summary $(PROC)lefse/strep_day10_clearance.design $(PROC)lefse/strep_day0_colonization.0.03.lefse_summary $(PROC)lefse/strep_day0_colonization.design $(PROC)lefse/strep_antibiotic.0.03.lefse_summary $(PROC)lefse/strep_antibiotic.design :code/functions.R\
									$(PROC)restore_metadata_clean.tsv\
									$(PROC)restore_taxonomy_clean.tsv\
									$(MOTHUR)sample.final.0.03.subsample.shared\
									code/run_lefse.R
	Rscript code/run_lefse.R

# Run SpiecEasi
$(PROC)cdiff_network.rds : code/functions.R\
									$(PROC)restore_metadata_clean.tsv\
									$(PROC)restore_taxonomy_clean.tsv\
									$(MOTHUR)sample.final.0.03.subsample.shared\
									code/run_spieceasi.R
	Rscript code/run_spieceasi.R

################################################################################
# Create figures
################################################################################

# Figure 1 is a schematic of the experimental design created with Affinity Designer

# Figure 2
submission/Figure_2.tiff : code/plot_figure_2.R\
							code/functions.R\
							$(PROC)restore_metadata_clean.tsv
	Rscript code/plot_figure_2.R

# Figure 3 
submission/Figure_3.tiff : code/plot_figure_3.R\
							code/functions.R\
							$(PROC)restore_metadata_clean.tsv
	Rscript code/plot_figure_3.R

# Figure 4, S3, and S4
submission/Figure_4.tiff submission/Figure_S3.tiff submission/Figure_S4.tiff : code/plot_figure_4.R\
							code/functions.R\
							$(PROC)restore_metadata_clean.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(MOTHUR)sample.final.groups.ave-std.summary\
							$(MOTHUR)sample.final.thetayc.0.03.lt.ave.dist\
							$(PROC)restore_taxonomy_clean.tsv
	Rscript code/plot_figure_4.R

# Figure 5
submission/Figure_5.tiff : code/plot_figure_5.R\
							code/functions.R\
							$(PROC)restore_metadata_clean.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(PROC)restore_taxonomy_clean.tsv\
							$(PROC)lefse/strep_day0_clearance_OTU.0.03.lefse_summary\
							$(PROC)lefse/strep_day0_clearance_OTU.design\
							$(PROC)lefse/strep_day10_clearance_OTU.0.03.lefse_summary\
							$(PROC)lefse/strep_day10_clearance_OTU.design\
							$(PROC)lefse/strep_day0_colonization_OTU.0.03.lefse_summary\
							$(PROC)lefse/strep_day0_colonization_OTU.design
	Rscript code/plot_figure_5.R

# Figure 6
submission/Figure_6.tiff : code/plot_figure_6.R\
							code/functions.R\
							code/ggnet2.R\
							$(PROC)cdiff_network.rds
	Rscript code/plot_figure_6.R

# Figure S1
submission/Figure_S1.tiff : code/plot_figure_S1.R\
							code/functions.R\
							$(PROC)restore_metadata_clean.tsv\
	Rscript code/plot_figure_S1.R

# Figure S2
submission/Figure_S2.tiff : code/plot_figure_S2.R\
							code/functions.R\
							$(PROC)restore_metadata_clean.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(MOTHUR)sample.final.groups.ave-std.summary\
							$(MOTHUR)sample.final.thetayc.0.03.lt.ave.dist\
							$(PROC)restore_taxonomy_clean.tsv\
							$(RAW)inocula_16S_qpcr_022521_Abs_Quant.txt
	Rscript code/plot_figure_S2.R

# Figure S5
submission/Figure_S5.tiff : code/plot_figure_S5.R\
							code/functions.R\
							$(PROC)restore_metadata_clean.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(PROC)restore_taxonomy_clean.tsv\
							$(PROC)lefse/strep_antibiotic.0.03.lefse_summary\
							$(PROC)lefse/strep_antibiotic.design
	Rscript code/plot_figure_S5.R

# Figure S6
submission/Figure_S6.tiff : code/plot_figure_S6.R\
							code/functions.R\
							$(PROC)restore_metadata_clean.tsv\
							$(MOTHUR)sample.final.0.03.subsample.shared\
							$(PROC)restore_taxonomy_clean.tsv
	Rscript code/plot_figure_S6.R

################################################################################
# Create manuscript
################################################################################


submission/manuscript.pdf submission/manuscript.docx : submission/manuscript.Rmd\
		submission/figure_1.tiff\
		submission/figure_2.tiff\
		submission/figure_3.tiff\
		submission/figure_4.tiff\
		submission/figure_5.tiff\
		submission/figure_6.tiff\
		submission/figure_S1.tiff\
		submission/figure_S2.tiff\
		submission/figure_S3.tiff\
		submission/figure_S4.tiff\
		submission/figure_S5.tiff\
		submission/figure_S6.tiff\
		submission/mbio.csl\
		submission/references.bib
	R -e 'library(rmarkdown);render("submission/manuscript.Rmd", output_format="all")'

