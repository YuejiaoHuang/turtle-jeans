# Script that takes multiple sequence alignment and calculates pairwise distance matrix

library(ape)

# read in current alignment file
al <- read.dna('BUSCO_alignments_trimmed_final_pis/1000at32523.fa', format = "fasta")

# compute pairwise distance matrix
dist <- dist.dna(al, model = "raw")

# convert to matrix
dist <- as.matrix(dist)

# write pariwise distance matrix into new directory
write.table(dist, 'dist_matrices/1000at32523.dist', row.names=TRUE, col.names=TRUE)

# how to read it into R again
# dist <- read.table('dist_matrices/1000at32523.dist')
