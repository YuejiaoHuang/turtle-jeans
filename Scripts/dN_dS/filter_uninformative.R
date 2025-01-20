# Script that takes multiple sequence alignment and copies it to new location
# if it contains at least 10 segregating sites.

library(ape)

# read in current alignment file that has been trimmed based on codons
al <- read.dna('BUSCO_alignments_trimmed/99895at32523.fa', format="fasta")

# get indices of segregating sites of alignment
n_seg_sites <- seg.sites(al)

# print number of segregating sites 
# (used later to get file with n for each gene and look at distribution)
print(length(n_seg_sites))

# filter out all genes that have less 10 segregating sites
if (length(n_seg_sites) > 9) {
  
  # filtering is done by only coping those alignments to the new location that fulfill criterion
  write.dna(al, 'BUSCO_alignments_trimmed_final/99895at32523.fa', 
            format="fasta", colsep = "")
  
}
