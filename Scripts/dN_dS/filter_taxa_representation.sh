# Check that all genes have at least 9 taxa

# For each gene
cat list_all.txt | while read gene
do
	# Get number of taxa in alignment as number of IDs (number of ">")
	n_taxa=($(grep ">" BUSCO_alignments_no_outgroups/${gene}.fa | wc -l))
	# if number of taxa > 8:
	if [ ${n_taxa} -gt 8 ]; then
		# keep gene
		cp BUSCO_alignments_no_outgroups/${gene}.fa BUSCO_alignments_enough_taxa/${gene}.fa
	fi
done
