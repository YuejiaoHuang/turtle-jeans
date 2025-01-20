# Script to make number of bases multiple of 3
# forces alignment into reading-frame, but removing first 0-2 bases

# For each gene
cat list_enough_taxa.txt | while read gene
do
	# Write whole sequence in one line and then go line by line
	seqkit seq -w 0 BUSCO_alignments_enough_taxa/${gene}.fa | while read line
	do
		# If line is a header starting with >
		if [[ $line == ">"* ]]; then
			# Just write header
            echo "$line"
        # Else the line is a full sequence
        else
        	# Get the number of bases in the sequence
            n_sites=${#line}
            # Check how many bases need to be trimmed, which is the remainder of dividing by 3
            to_trim=$((n_sites % 3))
            # If 1 base needs to be trimmed
            if (( to_trim == 1 )); then
            	# Write full sequence starting from second base
                echo "${line:1}"
            # If 2 bases need to be trimmed
            elif (( to_trim == 2 )); then
            	# Write full sequence starting from third base
                echo "${line:2}"
            # Else the remainder is 0 and no bases need to be removed
            else
            	# Write full sequence without any changes
                echo "$line"
            fi
        fi
    # Save in new fasta file
    done > BUSCO_alignments_correct_n_codons/${gene}.fa
done
