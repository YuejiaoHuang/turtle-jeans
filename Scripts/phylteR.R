#install.packages("phylter")
library(phylter)
library(tidyverse)
library(fishualize)
library(phytools)
library(AnnotationForge)
library(clusterProfiler)
library(tidyverse)
library(ape)
library(phangorn)
library(factoextra)
library(ggplot2)
library(TreeDist)
library(cowplot)
library(ggpubr)
library(ComplexHeatmap)
library(phangorn)
library(grid)
library(factoextra)
library(enrichplot)
library(ggtree)

#################
### LOAD DATA ###
#################

setwd("/Users/jule/Desktop/turtle-jeans")

# load gene trees
locus.trees <- read.tree('Data/genetrees.nwk')

# load gene names
names <- readLines("Data/list.txt")

# filter based on at least 10 tips (50 % of taxa present)
names <- names[Ntip(locus.trees) > 10]
locus.trees <- locus.trees[Ntip(locus.trees) > 10]

# load species tree
species_tree <- read.tree("Data/concord.cf.tree_concat.nwk")

# metadata for turtles
meta_turtles <- read_tsv("Data/metadata_habitat.tsv")
meta_turtles$Habitat_factor <- factor(meta_turtles$Habitat, 
                                      levels=c("Marine", "Rivers", "Lakes", 
                                               "Swamps", "Terrestrial", "Outgroup"))

meta_turtles2 <- read_tsv("Data/metadata_habitat_reptraits.tsv")
meta_turtles2$Habitat_factor <- factor(meta_turtles2$Microhabitat, 
                                      levels=c("Marine", "Aquatic", "Aquatic_Terrestrial", "Terrestrial", "Outgroup"))

# load species tree for plotting
species_tree_plot <- read.tree("Data/species_tree_plot.nwk")
# rename tips from ID to Name (nicer to read)
species_tree_plot$tip.label <- meta_turtles2$Species[match(species_tree_plot$tip.label, meta_turtles2$ID)]


###############
### COLOURS ###
###############

outlier_colour <- "#727149FF"
no_outlier_colour <- "#A2BAC5FF"
species_tree_colour <- "#685369"


##################################
### GENE TREES VS SPECIES TREE ###
##################################

### calculate the pariwise normalized RF distance 
normal_rf <- c()
for (i in seq(length(locus.trees))) {
  rf <- TreeDist::RobinsonFoulds(species_tree, locus.trees[i], normalize = T)
  normal_rf <- c(normal_rf,rf)
}

# Find the indices where RF equal to 1
high_rf_indices <- which(normal_rf == 1)

# Extract the RF distances that are greater than 0.9
high_rf_distances <- normal_rf[high_rf_indices]

# Get the names of the detected outlier genes with a large RF distance
og_candidate_gene_trees_vs_species_trees <- names[high_rf_indices]

# Density plot
normal_rf <- tibble(names=names, normal_rf=normal_rf)

# outliers as defined by 99% confidence interval
outliers_cutoff <- quantile(normal_rf$normal_rf, 0.99)
density_data_normal_rf <- density(normal_rf$normal_rf)
density_df_normal_rf <- tibble(normal_rf=density_data_normal_rf$x,
                               density=density_data_normal_rf$y)
density_df_normal_rf <- density_df_normal_rf %>% 
  mutate(outlier = case_when(normal_rf > outliers_cutoff ~ "Outlier",
                             TRUE ~ "No outlier"))

ggplot(density_df_normal_rf, aes(x=normal_rf, y=density, fill=outlier)) +
  geom_line(color = no_outlier_colour) +
  geom_ribbon(aes(ymin = 0, ymax = density), alpha = 0.8) +
  scale_fill_manual("Gene category",
                    values = c("No outlier" = no_outlier_colour, "Outlier" = outlier_colour)) + 
  geom_vline(xintercept = outliers_cutoff, linetype = "dashed", 
             color = "grey40", size = 0.4) +
  annotate("text", x = outliers_cutoff + 0.01, y = 0.3,  # Position text slightly below the peak
           label = "99% Cutoff", color = "grey40", hjust = -0.1, size=4) +
  xlab("Normalized RF distance") +
  ylab("Density") +
  theme_minimal()

ggsave("Results/density.pdf", width=7, height=7)


################################
### GENE TREES VS GENE TREES ###
################################

# also include species tree in analysis
locus.trees <- c(locus.trees, species_tree)
names <- c(names, "species_tree")

# RUNNN phylteR to get outliers
results <- phylter(locus.trees, gene.names = names, k=6)

# premade plots
# Get a summary: nb of outliers, gain in concordance, etc.
summary(results)

# Show the number of species in each gene, and how many per gene are outliers
plot(results, "genes") 

# Show the number of genes where each species is found, and how many are outliers
plot(results, "species") 

# get dists
distances <- results$Initial$RV



# Make outlier categories df
results$Final$Outliers

outlier_df <- as_tibble_col(names, column_name = "Name")
outlier_df$Outlier <- "No outlier"
outlier_df <- outlier_df %>% 
  mutate(Outlier = case_when(Name %in% results$Final$Outliers[,1] ~ "Outlier",
                             Name == "species_tree" ~ "Species tree",
                             TRUE ~ Outlier))

groups <- as.factor(outlier_df$Outlier)

# PCA
res.pca <- prcomp(distances, scale = TRUE)

outlier_colour <- "#8E8D71"
# ggplot PCA
pca <- res.pca$x[,1:2]
pca <- as_tibble(pca)
pca$Category <- groups
ggplot(data=pca, aes(x=PC1, y=PC2, colour=Category)) +
  geom_point() +
  theme_minimal() +
  xlab("PC1 (87.8%)") +
  ylab("PC2 (2.5%)") +
  scale_color_manual("Gene category", 
                     values=c(no_outlier_colour, outlier_colour, species_tree_colour)) 

ggsave("Results/PCA.pdf", width=7, height=7)

### densitytree plot

outlier_trees <- locus.trees[names %in% results$Final$Outliers[,1]]
densiTree(outlier_trees, col = c(rep('steelblue', 48), "black"), alpha = 0.15, type = "cladogram", 
          use.edge.length = F, width = 2,scaleX = 1)


##############################
### GO ENRICHMENT ANALYSIS ###
##############################

### turtles annotation db
go_all <- read.table("/Users/jule/Downloads/odb10v1_OG_xrefs.tab", 
                     header = F, col.names = c("GID","ontology","GO","Support_number"), sep = "\t")

### subset DB associated with biological process
go_bp <- go_all %>% 
  filter(ontology == "biological_process",
         str_detect(GO,"GO")) %>% 
  select(GO,GID) %>% 
  mutate(goterm = Term(GO))
go_bp_name <- go_bp %>% select(GO,GID)
go_bp_term <- go_bp %>% select(GO,goterm)

### subset DB associated with cellular component
go_cc <- go_all %>% 
  filter(ontology == "cellular_component",
         str_detect(GO,"GO")) %>% 
  select(GO,GID) %>% 
  mutate(goterm = Term(GO))
go_cc_name <- go_cc %>% select(GO,GID)
go_cc_term <- go_cc %>% select(GO,goterm)


### subset DB associated with molecular function
go_mf <- go_all %>% 
  filter(ontology == "molecular_function",
         str_detect(GO,"GO")) %>% 
  select(GO,GID) %>% 
  mutate(goterm = Term(GO))
go_mf_name <- go_mf %>% select(GO,GID)
go_mf_term <- go_mf %>% select(GO,goterm)


### DB with all ontologies

go_all_name <- go_all %>% filter(str_detect(GO,"GO")) %>% select(GO,GID)
go_all_term <- go_all %>% filter(str_detect(GO,"GO")) %>% select(GO,GID) %>% mutate(goterm = Term(GO)) %>% select(GO,goterm)

bp_enrich <- enricher(gene = og_candidate_gene_trees_vs_species_trees,
                      # pvalueCutoff = 100,
                      pAdjustMethod = "fdr",
                      # qvalueCutoff = 100,
                      TERM2GENE = go_bp_name,
                      TERM2NAME = go_bp_term)
cc_enrich <- enricher(gene = og_candidate_gene_trees_vs_species_trees,
                      # pvalueCutoff = 100,
                      pAdjustMethod = "fdr",
                      # qvalueCutoff = 100,
                      TERM2GENE = go_cc_name,
                      TERM2NAME = go_cc_term)
mf_enrich <- enricher(gene = og_candidate_gene_trees_vs_species_trees,
                      # pvalueCutoff = 100,
                      pAdjustMethod = "fdr",
                      # qvalueCutoff = 100,
                      TERM2GENE = go_mf_name,
                      TERM2NAME = go_mf_term)

dotplot(bp_enrich, showCategory=30)

bp_enrich2 <- pairwise_termsim(bp_enrich)
treeplot(bp_enrich2)

### visualization

# Add an 'ontology' column to each data frame to indicate the GO category
bp_df <- as.data.frame(bp_enrich) %>% mutate(ontology = "BP")
cc_df <- as.data.frame(cc_enrich) %>% mutate(ontology = "CC")
mf_df <- as.data.frame(mf_enrich) %>% mutate(ontology = "MF")

# Combine the data frames into one
combined_result <- bind_rows(bp_df, cc_df, mf_df)

# Select only necessary columns for plotting (for example, "Description", "Count", and "ontology")
# Make sure to adjust column names if different in your enricher output
plot_data <- combined_result %>%
  select(Description, Count,p.adjust, ontology) %>% 
  remove_missing(.)

ggplot(bp_df, aes(x=geneID, y=Description, fill=pvalue)) +
  geom_tile()

ggsave("./Results/GO_heatmap.pdf", height = 7, width = 9)


# OUTLIERS
outliers_og_subset <- results$Final$Outliers[,1]

names(outliers_og_subset) <- results$Final$Outliers[,2]

table(results$Final$Outliers[,2])

unique(outliers_og_subset)


Outliers_all <- results$Final$Outliers %>% 
  as_data_frame(.) %>% 
  left_join(.,meta_turtles, by = c("V2" = "Species")) %>% 
  remove_missing(.)

# Outliers_Freshwater <- Outliers_all %>% filter(Habitat == "Freshwater")
Outliers_Rivers <- Outliers_all %>% filter(Habitat == "Rivers")
Outliers_Swamps <- Outliers_all %>% filter(Habitat == "Swamps")
Outliers_Lakes <- Outliers_all %>% filter(Habitat == "Lakes")
Outliers_Terrestrial<- Outliers_all %>% filter(Habitat == "Terrestrial")
Outliers_Marine <- Outliers_all %>% filter(Habitat == "Marine")

enrich_plot <- function(df,plot_name,h,w){
  bp <- enricher(gene = unique(df$V1),
                 pAdjustMethod = "fdr",
                 TERM2GENE = go_bp_name,
                 TERM2NAME = go_bp_term)
  cc <- enricher(gene = unique(df$V1),
                 pAdjustMethod = "fdr",
                 TERM2GENE = go_cc_name,
                 TERM2NAME = go_cc_term)
  mf <- enricher(gene = unique(df$V1),
                 pAdjustMethod = "fdr",
                 TERM2GENE = go_mf_name,
                 TERM2NAME = go_mf_term)
  # Add an 'ontology' column to each data frame to indicate the GO category
  bp_df1 <- as.data.frame(bp) %>% mutate(ontology = "BP")
  cc_df1 <- as.data.frame(cc) %>% mutate(ontology = "CC")
  mf_df1 <- as.data.frame(mf) %>% mutate(ontology = "MF")
  
  # Combine the data frames into one
  combined_result1 <- bind_rows(bp_df1, cc_df1, mf_df1)
  
  # Select only necessary columns for plotting (for example, "Description", "p.adjust", "Count", and "ontology")
  plot_data1 <- combined_result1 %>%
    select(Description, Count, p.adjust, ontology) %>% 
    remove_missing(.)
  
  # Create the bar plot with facet_wrap
  ggplot(plot_data1, aes(x = Description, y = Count, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ ontology, scales = "free_y",nrow = 3) +
    theme_bw() +
    xlab(NULL) +
    ylab("Count") +
    scale_fill_gradient2(low = "grey", mid = "yellow", high = "red") +
    theme(strip.text = element_text(size = 12),  # Adjust facet label size
          axis.text = element_text(size = 10),   # Adjust axis text size
          axis.title = element_text(size = 14))  # Adjust axis title size
  
  ggsave(plot_name, height = h, width = w)
  
  
}
# enrich_plot(Outliers_Rivers,"Outliers_Rivers.png",10,9)
enrich_plot(Outliers_Swamps,"Outliers_Swamps.png",10,9)
enrich_plot(Outliers_Lakes,"Outliers_Lakes.png",10,9)
enrich_plot(Outliers_Marine,"Outliers_Marine.png",10,9)
# enrich_plot(Outliers_Terrestrial,"Outliers_Terrestrial.png",10,9)

################
### TREEPLOT ###
################


colours_classes5 <- fish(n=5,option="Balistoides_conspicillum", end=0.95, 
                         begin=0.3,direction=-1)
colours_classes6 <- c(colours_classes5, "grey")

group_habitat <- list(Marine=meta_turtles %>% filter(Habitat_factor=="Marine") %>% select(Species) %>% .$Species,
                      Rivers=meta_turtles %>% filter(Habitat_factor=="Rivers") %>% select(Species) %>% .$Species,
                      Lakes=meta_turtles %>% filter(Habitat_factor=="Lakes") %>% select(Species) %>% .$Species,
                      Swamps=meta_turtles %>% filter(Habitat_factor=="Swamps") %>% select(Species) %>% .$Species,
                      Terrestrial=meta_turtles %>% filter(Habitat_factor=="Terrestrial") %>% select(Species) %>% .$Species,
                      Outgroup=meta_turtles %>% filter(Habitat_factor=="Outgroup") %>% select(Species) %>% .$Species)

species_tree_plot <- groupOTU(species_tree_plot, group_habitat)
# plot for our tree
plot_tree <- ggtree::ggtree(species_tree_plot) + ggtree::xlim_tree(13)
plot_tree <- plot_tree %<+% meta_turtles +
  ggtree::geom_tiplab(size=3, offset=0.5, fontface = "italic") + 
  ggtree::geom_tippoint(aes(color=Habitat_factor)) + 
  scale_color_manual("Habitat", values=colours_classes6) +
  theme_tree2() +
  vexpand(0.01, direction = -1)
plot_tree

# should we remove the outgroups???
ggsave("Results/species_tree.pdf", width = 8, height = 5)

colours_classes6 <- c( "#7EA77DFF", "#0F3D5CFF", "grey", "#4C98B8FF", "#9DB327FF", "#DEE100FF")

plot_tree <- ggtree::ggtree(species_tree_plot, aes(color=group)) + ggtree::xlim_tree(13)
plot_tree <- plot_tree %<+% meta_turtles +
  ggtree::geom_tiplab(size=3, offset=0.5, fontface = "italic") + 
  ggtree::geom_tippoint(aes(color=Habitat_factor)) + 
  scale_color_manual("Habitat", values=colours_classes6) +
  theme_tree2() +
  vexpand(0.01, direction = -1)
plot_tree

# should we remove the outgroups???
ggsave("Results/species_tree2.pdf", width = 8, height = 5)


colours_classes5 <- c( "#4C98B8FF", "#9DB327FF", "#0F3D5CFF", "grey", "#DEE100FF")

group_habitat <- list(Marine=meta_turtles2 %>% filter(Habitat_factor=="Marine") %>% select(Species) %>% .$Species,
                      Aquatic=meta_turtles2 %>% filter(Habitat_factor=="Aquatic") %>% select(Species) %>% .$Species,
                      Aquatic_Terrestrial=meta_turtles2 %>% filter(Habitat_factor=="Aquatic_Terrestrial") %>% select(Species) %>% .$Species,
                      Terrestrial=meta_turtles2 %>% filter(Habitat_factor=="Terrestrial") %>% select(Species) %>% .$Species,
                      Outgroup=meta_turtles2 %>% filter(Habitat_factor=="Outgroup") %>% select(Species) %>% .$Species)

species_tree_plot <- groupOTU(species_tree_plot, group_habitat)

plot_tree <- ggtree::ggtree(species_tree_plot, aes(color=group)) + ggtree::xlim_tree(13)
plot_tree <- plot_tree %<+% meta_turtles2 +
  ggtree::geom_tiplab(size=3, offset=0.5, fontface = "italic") + 
  ggtree::geom_tippoint(aes(color=Habitat_factor)) + 
  scale_color_manual("Habitat", values=colours_classes5) +
  theme_tree2() +
  vexpand(0.01, direction = -1)
plot_tree
