#################
### LOAD DATA ###
#################

setwd("/Users/jule/Desktop/turtle-jeans")

# load gene trees
locus.trees <- read.tree('Data/genetrees.nwk')

# load gene names
names <- readLines("Data/list.txt")


### Drop outgroups
locus.trees <- drop.tip.multiPhylo(locus.trees, c("homSap", "galGal", "allMis"))

# filter based on at least 9 tips (50 % of taxa present)
names <- names[Ntip(locus.trees) >= 9]
locus.trees <- locus.trees[Ntip(locus.trees) >= 9]



# RUNNN phylteR to get outliers
results <- phylter(locus.trees_combined, gene.names = names_combined)


# get matrix for each gene that contains pairwise distances between species
matrices <- results$Initial$mat.data

# Optional scaling/normalization

## OUTLIER SPECIES
# for taxon 1: compare row for taxon1 between all genes
# ...

## OUTLIER GENES
# for build distribution for each gene (based on values for all pairwise distances)
#     plot all pariwise distances from upper diagonal of matrix
#     fit parameters for that distribution
# which distribution is different (based on comparing parameters)?


