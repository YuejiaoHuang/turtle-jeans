# Run R script to check that all genes have at least 10 segregating sites

# For each gene
cat list_enough_taxa.txt | while read gene
do
	
	# copy prepared R script
	cp filter_uninformative.R filter_uninformative_gene.R
	# replace initial gene ID with current gene ID in duplicated script
	sed -i '' "s/99895at32523/${gene}/g" filter_uninformative_gene.R
	# run R script to filter based on number of segregating sites
	Rscript filter_uninformative_gene.R

done