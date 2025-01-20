# Script to copy all relevant codeml output into one folder, rename by gene

# For each gene
cat list_enough_taxa.txt | while read gene
do

	# copy 2ML.dN and 2ML.dS into new folder, rename according to gene ID
	cp codeml_out/${gene}/2ML.dN dN_dS_data/${gene}_dN.dist
	cp codeml_out/${gene}/2ML.dS dN_dS_data/${gene}_dS.dist
done