# start with Trimmed_BUSCO_alignments - not codon aware trimming :(

# remove outgroups with remove_taxa.sh
./remove_taxa.sh

# remove genes with too little taxa (>=9)
./filter_taxa_representation.sh

# create new gene list to work with from now on
ls BUSCO_alignments_enough_taxa > list_enough_taxa.txt
sed -i '' "s/.fa//g" list_enough_taxa.txt

# remove 0-2 codons from the start of the alignment to make sure that number of codons is multiple of 3
./trim_n_codons.sh

# even though the alignments were not created with a codon-aware aligner, we assume for the next steps that all alignments are mostly in reading frame

# run macse to remove internal stop codons
./remove_stop_codons.sh

# run R trimming for codons with too many gaps (from macse) or where too many different amino acids (based on species)
./trim_bad_codons.sh

# Different filtering options, but parsimony informative sites should be best:

# filter out genes that are too short (< 500 sites)
# ./filter_short_genes.sh

# filter out genes that are not informative (< 10 segregating sites)
# ./filter_uninformative.sh > n_seg_sites.txt

# filter out genes that are not informative (< 10 parsimony informative sites)
./filter_parsimony_uninformative.sh > n_parsimony_sites.txt

# get final list of genes that should be used as input for codeml
ls BUSCO_alignments_trimmed_final_pis > list_final.txt
sed -i '' "s/.fa//g" list_final.txt
# diff list_enough_taxa.txt list_final.txt

# run codeML to get pairwise dN and dS
./codeml.sh

# copy all dN and dS into on folder
./sort_dN_dS_data.sh

# calculate pure distance matrices on final trimmed alignments (independent of dN/dS)
./calculate_pdm.sh