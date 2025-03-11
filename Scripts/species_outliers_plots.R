library(tidyverse)
library(fishualize)
library(ggvenn)
library(pheatmap)
library(reshape2)
library(ape)
library(ggtree)
library(deeptime)
library(phangorn)
library(ggrepel)

setwd("/Users/jule/Desktop/turtle-jeans")

colours_classes <- fish(n=1,option="Balistoides_conspicillum", end=0.9, 
                        begin=0.9)
colours_classes[2] <- fish(n=1,option="Balistoides_conspicillum", end=0.8, 
                           begin=0.8)
colours_classes[3] <- fish(n=1,option="Balistoides_conspicillum", end=0.4, 
                           begin=0.4)
colours_classes[4] <- fish(n=1,option="Balistoides_conspicillum", end=0.2, 
                           begin=0.2)

gene_names <- readLines("Scripts/dN_dS/list_final.txt")

# outliers_dn <- read.csv("Results/outliers_genes_dn_q95.csv")
# colnames(outliers_dn) <- c("n", "gene")
# outliers_ds <- read.csv("Results/outliers_genes_ds_q95.csv")
# colnames(outliers_ds) <- c("n", "gene")
# outliers_dnds <- read.csv("Results/outliers_genes_dnds_q95.csv")
# colnames(outliers_dnds) <- c("n", "gene")
# outliers_pure <- read.csv("Results/outliers_genes_pure_q95.csv")
# colnames(outliers_pure) <- c("n", "gene")

outliers_dn_species <- read.csv("Results/outliers_species_dn_q95.csv")
colnames(outliers_dn_species) <- c("n", "species", "gene")
outliers_ds_species <- read.csv("Results/outliers_species_ds_q95.csv")
colnames(outliers_ds_species) <- c("n", "species", "gene")
outliers_dnds_species <- read.csv("Results/outliers_species_dnds_q95.csv")
colnames(outliers_dnds_species) <- c("n", "species", "gene")
outliers_pure_species <- read.csv("Results/outliers_species_pure_q95.csv")
colnames(outliers_pure_species) <- c("n", "species", "gene")


#########################
### SPECIES TREE PLOT ###
#########################

# colours_classes5 <- fish(n=5,option="Balistoides_conspicillum", end=0.95, 
#                          begin=0.3,direction=-1)
# colours_classes4 <- colours_classes5[1:4]
# 
# c("#0F3D5CFF", "#4C98B8FF", "#7EA77DFF", "#9DB327FF", "#DEE100FF")
# 
# colour_marine <- colours_classes4[1]
# colour_aquatic <- colours_classes4[2]
# colour_aquatic_terrestrial <- colours_classes4[3]
# colour_terrestrial <- colours_classes4[4]
# 

meta_turtles <- read_tsv("Data/metadata_habitat_reptraits.tsv")
meta_turtles$Habitat_factor <- factor(meta_turtles$Microhabitat, 
                                      levels=c("Marine", "Aquatic", 
                                               "Aquatic_Terrestrial", "Terrestrial", 
                                               "Outgroup"))
meta_turtles <- meta_turtles %>% filter(Microhabitat != "Outgroup")

# load species tree with branch lengths from Thomson et al. (2021)
species_tree_plot <- read.nexus("Data/bd.mcc.median_heights.tre")
species_tree_plot <- drop.tip(species_tree_plot, 
                              setdiff(species_tree_plot$tip.label, meta_turtles$ID))
species_tree_plot$tip.label <- meta_turtles$Species[
  match(species_tree_plot$tip.label, meta_turtles$ID)]

# plot for our tree
plot_tree <- ggtree::ggtree(species_tree_plot)
plot_tree <- plot_tree +
  ggtree::geom_tiplab(size=3, offset=0.5, fontface = "italic") + 
  theme_tree2() +
  coord_geo(xlim = c(-240, 100), ylim = c(-0.5, Ntip(species_tree_plot)+2),
            neg = TRUE, abbrv = list(TRUE, FALSE), dat=list("epochs", "periods"),
            pos = list("bottom", "bottom"), size = "auto",
            height = list(unit(1, "lines"), unit(1, "lines"))) +
  scale_x_continuous(breaks = seq(-240, 0, 20), labels = abs(seq(-240, 0, 20)))
revts(plot_tree)

ggsave("Results/species_tree_branch_lengths.pdf", width = 8, height = 5)



###############################
### OVERLAP BETWEEN SPECIES ###
###############################

species_list <- meta_turtles %>% filter(Microhabitat != "Outgroup")
species_list <- species_list$ID

df_species_id <- meta_turtles %>% dplyr::select(c("Species", "ID"))

species_ordered_list <- c("Malaclemys terrapin",
                          "Chrysemys picta",
                          "Terrapene mexicana",
                          "Platysternon megacephalum",
                          "Cuora mccordi",
                          "Cuora amboinensis",
                          "Gopherus agassizii",
                          "Chelonoidis abingdonii",
                          "Dermochelys coriacea",
                          "Chelonia mydas",
                          "Dermatemys mawii",
                          "Chelydra serpentina",
                          "Carettochelys insculpta",
                          "Pelodiscus sinensis",
                          "Pelusios castaneus",
                          "Podocnemis expansa",
                          "Emydura subglobosa",
                          "Mesoclemmys tuberculata")

# DN
data_heatmap_dn <- outliers_dn_species %>% dplyr::select(c("species", "gene"))

matrix_heatmap_dn <- matrix(0, nrow = length(species_list), ncol = length(species_list))
rownames(matrix_heatmap_dn) <- species_list
colnames(matrix_heatmap_dn) <- species_list

for (i in 1:length(species_list)) {
  for (j in 1:length(species_list)) {
    genes_i <- data_heatmap_dn$gene[data_heatmap_dn$species == species_list[i]]
    genes_j <- data_heatmap_dn$gene[data_heatmap_dn$species == species_list[j]]
    matrix_heatmap_dn[i, j] <- length(intersect(genes_i, genes_j))
    
    # if same species
    if (i == j) {
      other_genes <- data_heatmap_dn$gene[data_heatmap_dn$species != species_list[i]]
      unique_genes <- setdiff(genes_i, other_genes)
      matrix_heatmap_dn[i, j] <- length(unique_genes)
    }
  }
}

df_heatmap_dn <- melt(matrix_heatmap_dn)
df_heatmap_dn <- merge(df_heatmap_dn, df_species_id, by.x="Var1", by.y="ID", all.x=T)
colnames(df_heatmap_dn) <- c("Var1", "Var2", "value", "Species1")
df_heatmap_dn <- merge(df_heatmap_dn, df_species_id, by.x="Var2", by.y="ID", all.x=T)
colnames(df_heatmap_dn) <- c("Var2", "Var1", "value", "Species1", "Species2")
df_heatmap_dn$Species1 <- factor(df_heatmap_dn$Species1, 
                                 levels=species_ordered_list)
df_heatmap_dn$Species2 <- factor(df_heatmap_dn$Species2, 
                                 levels=species_ordered_list[rev(1:length(species_ordered_list))])

ggplot(df_heatmap_dn, aes(x = Species1, y = Species2, fill = value)) +
  geom_tile() +
  geom_text(aes(label=value), color = "white") +
  labs(title = "dN",
       x = "Species",
       y = "Species") +
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#72315C", "#A6A867")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("Results/heatmap_dn.pdf", width = 8, height = 5)


# DS
data_heatmap_ds <- outliers_ds_species %>% dplyr::select(c("species", "gene"))

matrix_heatmap_ds <- matrix(0, nrow = length(species_list), ncol = length(species_list))
rownames(matrix_heatmap_ds) <- species_list
colnames(matrix_heatmap_ds) <- species_list

for (i in 1:length(species_list)) {
  for (j in 1:length(species_list)) {
    genes_i <- data_heatmap_ds$gene[data_heatmap_ds$species == species_list[i]]
    genes_j <- data_heatmap_ds$gene[data_heatmap_ds$species == species_list[j]]
    matrix_heatmap_ds[i, j] <- length(intersect(genes_i, genes_j))
    
    # if same species
    if (i == j) {
      other_genes <- data_heatmap_ds$gene[data_heatmap_ds$species != species_list[i]]
      unique_genes <- setdiff(genes_i, other_genes)
      matrix_heatmap_ds[i, j] <- length(unique_genes)
    }
  }
}

df_heatmap_ds <- melt(matrix_heatmap_ds)
df_heatmap_ds <- merge(df_heatmap_ds, df_species_id, by.x="Var1", by.y="ID", all.x=T)
colnames(df_heatmap_ds) <- c("Var1", "Var2", "value", "Species1")
df_heatmap_ds <- merge(df_heatmap_ds, df_species_id, by.x="Var2", by.y="ID", all.x=T)
colnames(df_heatmap_ds) <- c("Var2", "Var1", "value", "Species1", "Species2")
df_heatmap_ds$Species1 <- factor(df_heatmap_ds$Species1, 
                                 levels=species_ordered_list)
df_heatmap_ds$Species2 <- factor(df_heatmap_ds$Species2, 
                                 levels=species_ordered_list[rev(1:length(species_ordered_list))])

ggplot(df_heatmap_ds, aes(x = Species1, y = Species2, fill = value)) +
  geom_tile() +
  geom_text(aes(label=value), color = "white") +
  labs(title = "dS",
       x = "Species",
       y = "Species") +
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#72315C", "#A6A867")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("Results/heatmap_ds.pdf", width = 8, height = 5)

# DNDS
data_heatmap_dnds <- outliers_dnds_species %>% dplyr::select(c("species", "gene"))

matrix_heatmap_dnds <- matrix(0, nrow = length(species_list), ncol = length(species_list))
rownames(matrix_heatmap_dnds) <- species_list
colnames(matrix_heatmap_dnds) <- species_list

for (i in 1:length(species_list)) {
  for (j in 1:length(species_list)) {
    genes_i <- data_heatmap_dnds$gene[data_heatmap_dnds$species == species_list[i]]
    genes_j <- data_heatmap_dnds$gene[data_heatmap_dnds$species == species_list[j]]
    matrix_heatmap_dnds[i, j] <- length(intersect(genes_i, genes_j))
    
    # if same species
    if (i == j) {
      other_genes <- data_heatmap_dnds$gene[data_heatmap_dnds$species != species_list[i]]
      unique_genes <- setdiff(genes_i, other_genes)
      matrix_heatmap_dnds[i, j] <- length(unique_genes)
    }
  }
}

df_heatmap_dnds <- melt(matrix_heatmap_dnds)
df_heatmap_dnds <- merge(df_heatmap_dnds, df_species_id, by.x="Var1", by.y="ID", all.x=T)
colnames(df_heatmap_dnds) <- c("Var1", "Var2", "value", "Species1")
df_heatmap_dnds <- merge(df_heatmap_dnds, df_species_id, by.x="Var2", by.y="ID", all.x=T)
colnames(df_heatmap_dnds) <- c("Var2", "Var1", "value", "Species1", "Species2")

df_heatmap_dnds$Species1 <- factor(df_heatmap_dnds$Species1, 
                                   levels=species_ordered_list)
df_heatmap_dnds$Species2 <- factor(df_heatmap_dnds$Species2, 
                                   levels=species_ordered_list[rev(1:length(species_ordered_list))])

ggplot(df_heatmap_dnds, aes(x = Species1, y = Species2, fill = value)) +
  geom_tile() +
  geom_text(aes(label=value), color = "white") +
  labs(title = "dN/dS",
       x = "Species",
       y = "Species") +
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#72315C", "#A6A867")) +
  # scale_fill_gradientn(name = "Overlapping genes", colours=c("#632A50", "#BDBF09")) +
  # scale_fill_gradientn(name = "Overlapping genes", colours=c("#632A50", "#60935D")) +
  # scale_fill_gradientn(name = "Overlapping genes", colours=c("#632A50", "#188FA7")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("Results/heatmap_dnds.pdf", width = 8, height = 5)

# PURE
data_heatmap_pure <- outliers_pure_species %>% dplyr::select(c("species", "gene"))

matrix_heatmap_pure <- matrix(0, nrow = length(species_list), ncol = length(species_list))
rownames(matrix_heatmap_pure) <- species_list
colnames(matrix_heatmap_pure) <- species_list

for (i in 1:length(species_list)) {
  for (j in 1:length(species_list)) {
    genes_i <- data_heatmap_pure$gene[data_heatmap_pure$species == species_list[i]]
    genes_j <- data_heatmap_pure$gene[data_heatmap_pure$species == species_list[j]]
    matrix_heatmap_pure[i, j] <- length(intersect(genes_i, genes_j))
    
    # if same species
    if (i == j) {
      other_genes <- data_heatmap_pure$gene[data_heatmap_pure$species != species_list[i]]
      unique_genes <- setdiff(genes_i, other_genes)
      matrix_heatmap_pure[i, j] <- length(unique_genes)
    }
  }
}

df_heatmap_pure <- melt(matrix_heatmap_pure)
df_heatmap_pure <- merge(df_heatmap_pure, df_species_id, by.x="Var1", by.y="ID", all.x=T)
colnames(df_heatmap_pure) <- c("Var1", "Var2", "value", "Species1")
df_heatmap_pure <- merge(df_heatmap_pure, df_species_id, by.x="Var2", by.y="ID", all.x=T)
colnames(df_heatmap_pure) <- c("Var2", "Var1", "value", "Species1", "Species2")

df_heatmap_pure$Species1 <- factor(df_heatmap_pure$Species1, 
                                   levels=species_ordered_list)
df_heatmap_pure$Species2 <- factor(df_heatmap_pure$Species2, 
                                   levels=species_ordered_list[rev(1:length(species_ordered_list))])

ggplot(df_heatmap_pure, aes(x = Species1, y = Species2, fill = value)) +
  geom_tile() +
  geom_text(aes(label=value), color = "white") +
  labs(title = "PURE",
       x = "Species",
       y = "Species") +
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#72315C", "#A6A867")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("Results/heatmap_pure.pdf", width = 8, height = 5)



### PAIRWISE DIVERGENCE TIMES

# just tips
div_times <- cophenetic.phylo(species_tree_plot)
div_times[lower.tri(div_times)] <- NA

df_div_times <- melt(div_times)
colnames(df_div_times) <- c("Species1", "Species2", "divergence_time")
# change from distance to divergence time
df_div_times <- df_div_times %>% filter(!is.na(divergence_time))
df_div_times$divergence_time <- df_div_times$divergence_time/2
df_div_times$Species_pair <- with(df_div_times, paste0(Species1,", ", Species2))

df_heatmap_dnds$Species_pair <- with(df_heatmap_dnds, paste0(Species1,", ", Species2))
df_heatmap_dnds_reduced <- df_heatmap_dnds%>% 
  dplyr::select(Species_pair, value)
df_div_times <- merge(df_div_times, df_heatmap_dnds_reduced, by="Species_pair")
df_div_times <- df_div_times %>% filter(Species1 != Species2)

ggplot(df_div_times, aes(x=divergence_time, y=value, label=Species_pair)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Species pairs",
       x = "Divergence time",
       y = "Number of overlapping genes") +
  


### INTERNAL NODES

outliers_dnds_species <- merge(outliers_dnds_species, df_species_id, by.x="species", by.y="ID", all.x=T)
# divergence times all nodes
div_times_nodes <- dist.nodes(species_tree_plot)[19,]
div_times_nodes <- floor(max(div_times_nodes) - div_times_nodes)

# get all internal nodes
internal_nodes <- (length(species_tree_plot$tip.label) + 1):max(species_tree_plot$edge)

get_tips <- function(tree, node) {
  tips <- Descendants(tree, node, type = "tips")[[1]]
  # return species names for all tips
  return(tree$tip.label[tips])
}

get_overlapping_genes <- function(tips, genes_species) {
  tip_genes <- genes_species %>% filter(Species %in% tips)
  split_gene_lists <- split(tip_genes$gene, tip_genes$Species)
  overlapping_genes <- Reduce(intersect, split_gene_lists)
  return(length(overlapping_genes))
}

results <- data.frame(node = integer(), divergence_time = numeric(), num_overlapping_genes = integer())

for (node in internal_nodes) {
  tips <- get_tips(species_tree_plot, node)
  overlapping_genes <- get_overlapping_genes(tips, outliers_dnds_species)
  results <- rbind(results, data.frame(
    node = node,
    divergence_time = div_times_nodes[node],
    num_overlapping_genes = overlapping_genes
  ))
}

ggplot(results, aes(x=divergence_time, y=num_overlapping_genes)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Internal nodes",
     x = "Divergence time",
     y = "Number of outlier genes present in all associated tips")
