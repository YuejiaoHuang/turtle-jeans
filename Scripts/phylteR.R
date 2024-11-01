#install.packages("phylter")
library(phylter)
library(tidyverse)
library(fishualize)
library(phytools)

# load gene trees
locus.trees <- read.tree('gene_species_trees.nwk')
locus.trees_filtered <- locus.trees[Ntip(locus.trees) > 10]

# load gene names
names <- readLines("gene_trees/list.txt")
names_filtered <- names[Ntip(locus.trees) > 10]

# name trees
names(locus.trees_filtered) <- names_filtered


# RUNNN
results <- phylter(locus.trees_filtered, gene.names = names_filtered, k=6)

# premade plots
# Get a summary: nb of outliers, gain in concordance, etc.
summary(results)

# Show the number of species in each gene, and how many per gene are outliers
plot(results, "genes") 

# Show the number of genes where each species is found, and how many are outliers
plot(results, "species") 

# save
write.phylter(results, file = "phylter.out")

# get dists
distances <- results$Initial$RV



# Make outlier categories df
results$Final$Outliers

outlier_df <- as_tibble_col(names_filtered, column_name = "Name")
outlier_df$Outlier <- "No outlier"
outlier_df <- outlier_df %>% 
  mutate(Outlier = case_when(Name %in% results$Final$Outliers[,1] ~ "Outlier",
                             Name == "species_tree" ~ "Species tree",
                             TRUE ~ Outlier))


# PCA

groups <- as.factor(outlier_df$Outlier)

library(factoextra)
res.pca <- prcomp(distances, scale = TRUE)

# premade plot
# fviz_pca_ind(res.pca, geom = c("point"), col.ind = groups, palette = c("#4694B8FF", "#93AC38FF", "#FCEF0FFF"), legend.title = "Outlier")

colours <- c("#A2BAC5FF", "#727149FF", "#4E4139FF")
# ggplot PCA
pca <- res.pca$x[,1:2]
pca <- as_tibble(pca)
pca$Category <- groups
ggplot(data=pca, aes(x=PC1, y=PC2, colour=Category)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual("Category", values=colours) 

### 2x densitytree plots for 2 sets of outliers :'(

outlier_trees <- locus.trees_filtered[names_filtered %in% results$Final$Outliers[,1]]
densityTree(outlier_trees, colors="#669ABFFF", alpha=0.9, method="plotTree", 
            fix.depth=T, use.edge.length=F, compute.consensus=F, 
            use.gradient=FALSE, show.axis=TRUE,type="cladogram",nodes="intermediate")

write.tree(outlier_trees, file = "outlier_trees.nwk")


library(phangorn)
library(ggplot2)
library(gridExtra)
library(grid)

# First, save the ggplot as a separate plot object
pca_plot <- ggplot(data=pca, aes(x=PC1, y=PC2, colour=Category)) +
  geom_point() +
  theme_minimal() +
  xlab("PC1 (87.8%)") +
  ylab("PC2 (2.5%)") +
  scale_color_manual("Category", values=colours) 

# Start the PNG output
png("pca_plot.png", width = 700, height = 700, res = 150)
# Create a grid layout with grid.arrange
grid.arrange(
  pca_plot
)
# Close the PNG device
dev.off()


png("densiTree_plot.png", width = 700, height = 700, res = 150)
densiTree(outlier_trees, col = c(rep('steelblue', 48), "black"), alpha = 0.15, type = "cladogram", use.edge.length = F, width = 2,scaleX = 1)
#densiTree(merged_trees, col = c("black", rep('steelblue', 7)), alpha = 0.7, use.edge.length = F, width = 2,scaleX = 1)
dev.off()
