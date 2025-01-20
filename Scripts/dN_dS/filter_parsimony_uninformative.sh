# Run R script to check that all genes have at least 10 parsimony informative sites

# For each gene
cat list_enough_taxa.txt | while read gene
do
	
	# copy prepared R script
	cp filter_parsimony_uninformative.R filter_parsimony_uninformative_gene.R
	# replace initial gene ID with current gene ID in duplicated script
	sed -i '' "s/99895at32523/${gene}/g" filter_parsimony_uninformative_gene.R
	# run R script to filter based on number of parsimony informative sites
	Rscript filter_parsimony_uninformative_gene.R

done