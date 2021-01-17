################################################################################
# (C) 2017 Tyler W. H. Backman
# Makefile to load PKS data
################################################################################

SHELL := /bin/bash
dataUrl = https://clustercad.jbei.org/static/data
username = pkscluster
password = pkscluster

all: data/generate_aa_plots_placeholder 

clean:
	-rm -rf data/mibigdata_clustercad_production.tgz
	-rm -rf data/antismash
	-rm -rf data/blast
	-rm -rf data/mibig
	-rm data/database_loaded_placeholder
	-rm data/clustercad_corrections_placeholder 
	-rm data/clustercad_clean_placeholder
	-rm data/process_aa_analysis_placeholder
	-rm data/generate_aa_plots_placeholder
	-rm data/aa_sequence_analysis/plotCache.tgz

# download PKS data from online
data/mibigdata_clustercad_production.tgz:
	mkdir -p data
	wget --user=$(username) --password=$(password) $(dataUrl)/mibigdata_clustercad_production.tgz -O $@ --no-check-certificate
	tar xvfz data/mibigdata_clustercad_production.tgz

# download plot cache
data/aa_sequence_analysis/plotCache.tgz:
	mkdir -p data/aa_sequence_analysis/plotCache 
	wget --user=$(username) --password=$(password) $(dataUrl)/plotCache.tgz -O data/aa_sequence_analysis/plotCache.tgz --no-check-certificate
	cd data/aa_sequence_analysis/ && tar xvfz plotCache.tgz

# aggregate chemical structures
data/compounds/all_known_products.p: aggregate_final_structures.py data/mibigdata_clustercad_production.tgz
	mkdir -p data/compounds
	./$< 

# load database
data/database_loaded_placeholder: antismash_to_database.py data/compounds/all_known_products.p
	./$<
	touch $@

# load corrections
data/clustercad_corrections_placeholder: clustercad_corrections.py data/database_loaded_placeholder
	./$<
	touch $@

# clean database
data/clustercad_clean_placeholder: clustercad_clean.py data/clustercad_corrections_placeholder
	./$<
	touch $@

# generate blast database (Reviewed and all)
data/blast/clustercad_subunits_reviewed: blast_database.py data/clustercad_clean_placeholder 
	mkdir -p data/blast
	./$< --name='clustercad_subunits_reviewed' --include_all=0
	cd data/blast && makeblastdb -in clustercad_subunits_reviewed -parse_seqids -dbtype prot 

# generate blast database (This is just a header.. for makefile purposes)
data/blast/clustercad_subunits_all: blast_database.py data/clustercad_clean_placeholder 
	mkdir -p data/blast
	./$< --name='clustercad_subunits_all' --include_all=1
	cd data/blast && makeblastdb -in clustercad_subunits_all -parse_seqids -dbtype prot

# load AA annotations into database
data/process_aa_analysis_placeholder: process_aa_analysis.py data/blast/clustercad_subunits_reviewed data/blast/clustercad_subunits_all
	./$<
	touch $@

# loop generating plots until an error occurs, then stop
# script will exit with error code 1 if there are no
# plots left to generate
# this is a temporary hack to overcome a memory leak in the plot
# code
data/generate_aa_plots_placeholder: generate_aa_plots.py data/process_aa_analysis_placeholder data/aa_sequence_analysis/plotCache.tgz
	mkdir -p data/aa_sequence_analysis/plotCache
	while [ $$? -eq 0 ]; do \
		./$<; \
	done; \
	touch $@ 
