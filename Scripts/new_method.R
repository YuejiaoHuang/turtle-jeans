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


# get matrix for each gene that contains pairwise distances between species
matrices <- results$Initial$mat.data

# Optional scaling/normalization

## OUTLIER SPECIES
# for taxon 1: compare row for taxon1 between all genes
# ...

# for taxon 1: take mean distance to other taxa for each gene
# plot all means
# take outliers from that distribution
# -> Outlier genes for taxa 1

all <- data.frame()

for(gene_name in names){
  print(gene_name)
  matrix <- matrices[[gene_name]]
  sums <- rowSums(matrix)
  means <- sums/(ncol(matrix)-1)
  frame <- data.frame(name = names(means), mean = means, gene = gene_name,
                      row.names = NULL)
  all <- rbind(all,frame)
}


library(reshape2)
reshape <- dcast(all, name ~ gene, value.var = 'mean')
all2 <- reshape[,-1]
rownames(all2) <- reshape[,1]

for(i in tips){
  print(i)
  plot(density(all$mean[all$name == i]),main = i)
}

## OUTLIER GENES
# for build distribution for each gene (based on values for all pairwise distances)
#     plot all pariwise distances from upper triangle of matrix
#     fit parameters for that distribution ISSUE: WHICH DISTRIBUTIONS, HOW TO COMPARE DISTRIBUTIONS
# which distribution is different (based on comparing parameters)?

# get upper triangle
matrices_example <- matrices$`100088at32523`
matrices_example_vector<-as.vector(matrices_example[upper.tri(matrices_example)])

data_frame(val = matrices_example_vector) %>%
  ggplot(., aes(val)) + 
  geom_density()
