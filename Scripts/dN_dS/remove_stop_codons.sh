# Script to get rid of internal stop codons with macse (replace with NNN)

# For each gene
cat list_enough_taxa.txt | while read gene
do
	
	# run macse to replace internal stop codons in given alignment with NNN
	java -jar macse_v2.07.jar -prog exportAlignment -align BUSCO_alignments_correct_n_codons/${gene}.fa -codonForInternalStop NNN -codonForFinalStop --- -codonForInternalFS --- -charForRemainingFS - -out_NT BUSCO_alignments_no_stop/${gene}.fa -out_AA no_stop_AA/${gene}.fa
done