# Script to remove a list of taxa (here outgroups) from multiple sequence alignment

# For each gene
cat list_all.txt | while read gene
do
	# Remove all outgroups listed in outgroups.txt
	# grep -v inverts: keep all IDs that do not match to outgroups.txt
	seqkit grep -v -f outgroups.txt Trimmed_BUSCO_alignments/${gene}.trim > BUSCO_alignments_no_outgroups/${gene}.fa
done
