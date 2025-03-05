library(tidyverse)
library(ape)
library(phylter)
library(TreeTools)



#################
### LOAD DATA ###
#################

#setwd("/Users/jule/Desktop/turtle-jeans")


# load gene trees
locus.trees <- read.tree('Data/genetrees.nwk')

# load gene names
names <- readLines("Data/list.txt")


### Drop outgroups
locus.trees <- drop.tip.multiPhylo(locus.trees, c("homSap", "galGal", "allMis"))
tips <- AllTipLabels(locus.trees)

# filter based on at least 9 tips (50 % of taxa present)
names <- names[Ntip(locus.trees) >= 9]
locus.trees <- locus.trees[Ntip(locus.trees) >= 9]

# RUNNN phylteR to get outliers
results <- phylter(locus.trees, gene.names = names)

# OUTLIERS
outliers_og_subset <- results$Final$Outliers[,1]

names(outliers_og_subset) <- results$Final$Outliers[,2]

table(results$Final$Outliers[,2])

outliers <- unique(outliers_og_subset)

write.csv(outliers,'Results/phylter_outliers.csv')

# get matrix for each gene that contains pairwise distances between species
matrices <- results$Initial$mat.data

summary(results)





### MAKE VENNS 
