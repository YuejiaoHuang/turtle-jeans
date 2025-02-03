library(tidyverse)
library(clusterProfiler)
library(AnnotationForge)
library(enrichplot)
library(fishualize)
library(ggvenn)
library(pheatmap)
library(reshape2)
library(ape)
library(ggtree)


setwd("/Users/jule/Desktop/turtle-jeans")

colours_classes <- fish(n=1,option="Balistoides_conspicillum", end=0.9, 
                        begin=0.9)
colours_classes[2] <- fish(n=1,option="Balistoides_conspicillum", end=0.8, 
                           begin=0.8)
colours_classes[3] <- fish(n=1,option="Balistoides_conspicillum", end=0.4, 
                           begin=0.4)
colours_classes[4] <- fish(n=1,option="Balistoides_conspicillum", end=0.2, 
                           begin=0.2)

gene_names <- readLines("Data/list.txt")

outliers_dn <- read.csv("Results/outliers_genes_dn_q95.csv")
colnames(outliers_dn) <- c("n", "gene")
outliers_ds <- read.csv("Results/outliers_genes_ds_q95.csv")
colnames(outliers_ds) <- c("n", "gene")
outliers_dnds <- read.csv("Results/outliers_genes_dnds_q95.csv")
colnames(outliers_dnds) <- c("n", "gene")
outliers_pure <- read.csv("Results/outliers_genes_pure_q95.csv")
colnames(outliers_pure) <- c("n", "gene")

outliers_dn_species <- read.csv("Results/outliers_species_dn_q95.csv")
colnames(outliers_dn_species) <- c("n", "species", "gene")
outliers_ds_species <- read.csv("Results/outliers_species_ds_q95.csv")
colnames(outliers_ds_species) <- c("n", "species", "gene")
outliers_dnds_species <- read.csv("Results/outliers_species_dnds_q95.csv")
colnames(outliers_dnds_species) <- c("n", "species", "gene")
outliers_pure_species <- read.csv("Results/outliers_species_pure_q95.csv")
colnames(outliers_pure_species) <- c("n", "species", "gene")

outliers_dn_species_allgene <- read.csv("Results/allgenes_species_dn_q95.csv")
colnames(outliers_dn_species_allgene) <- c("n", "species", "gene")
outliers_ds_species_allgene <- read.csv("Results/allgenes_species_dnds_q95.csv")
colnames(outliers_ds_species_allgene) <- c("n", "species", "gene")
outliers_dnds_species_allgene <- read.csv("Results/allgenes_species_ds_q95.csv")
colnames(outliers_dnds_species_allgene) <- c("n", "species", "gene")
outliers_pure_species_allgene <- read.csv("Results/allgenes_species_pure_q95.csv")
colnames(outliers_pure_species_allgene) <- c("n", "species", "gene")

meta_turtles <- read_tsv("Data/metadata_habitat_reptraits.tsv")
meta_turtles$Habitat_factor <- factor(meta_turtles$Microhabitat, 
                                      levels=c("Marine", "Aquatic", 
                                               "Aquatic_Terrestrial", "Terrestrial", 
                                               "Outgroup"))

# merge outliers, species and habitat
Outliers_all_dn <- outliers_dn_species %>%
  left_join(.,meta_turtles, by = c("species" = "ID")) %>% 
  remove_missing(.)

Outliers_Marine_dn <- Outliers_all_dn %>% 
  filter(Microhabitat == "Marine")
Outliers_Aquatic_dn <- Outliers_all_dn %>% 
  filter(Microhabitat == "Aquatic")
Outliers_Aquatic_Terrestrial_dn <- Outliers_all_dn %>% 
  filter(Microhabitat == "Aquatic_Terrestrial")
Outliers_Terrestrial_dn <- Outliers_all_dn %>% 
  filter(Microhabitat == "Terrestrial")

Outliers_all_ds <- outliers_ds_species %>%
  left_join(.,meta_turtles, by = c("species" = "ID")) %>% 
  remove_missing(.)

Outliers_Marine_ds <- Outliers_all_ds %>% 
  filter(Microhabitat == "Marine")
Outliers_Aquatic_ds <- Outliers_all_ds %>% 
  filter(Microhabitat == "Aquatic")
Outliers_Aquatic_Terrestrial_ds <- Outliers_all_ds %>% 
  filter(Microhabitat == "Aquatic_Terrestrial")
Outliers_Terrestrial_ds <- Outliers_all_ds %>% 
  filter(Microhabitat == "Terrestrial")

Outliers_all_dnds <- outliers_dnds_species %>%
  left_join(.,meta_turtles, by = c("species" = "ID")) %>% 
  remove_missing(.)

Outliers_Marine_dnds <- Outliers_all_dnds %>% 
  filter(Microhabitat == "Marine")
Outliers_Aquatic_dnds <- Outliers_all_dnds %>% 
  filter(Microhabitat == "Aquatic")
Outliers_Aquatic_Terrestrial_dnds <- Outliers_all_dnds %>% 
  filter(Microhabitat == "Aquatic_Terrestrial")
Outliers_Terrestrial_dnds <- Outliers_all_dnds %>% 
  filter(Microhabitat == "Terrestrial")

Outliers_all_pure <- outliers_pure_species %>%
  left_join(.,meta_turtles, by = c("species" = "ID")) %>% 
  remove_missing(.)

Outliers_Marine_pure <- Outliers_all_pure %>% 
  filter(Microhabitat == "Marine")
Outliers_Aquatic_pure <- Outliers_all_pure %>% 
  filter(Microhabitat == "Aquatic")
Outliers_Aquatic_Terrestrial_pure <- Outliers_all_pure %>% 
  filter(Microhabitat == "Aquatic_Terrestrial")
Outliers_Terrestrial_pure <- Outliers_all_pure %>% 
  filter(Microhabitat == "Terrestrial")


# merge ALLGENES outliers, species and habitat
Outliers_all_dn_allgene <- outliers_dn_species_allgene %>%
  left_join(.,meta_turtles, by = c("species" = "ID")) %>% 
  remove_missing(.)

Outliers_Marine_dn_allgene <- Outliers_all_dn_allgene %>% 
  filter(Microhabitat == "Marine")
Outliers_Aquatic_dn_allgene <- Outliers_all_dn_allgene %>% 
  filter(Microhabitat == "Aquatic")
Outliers_Aquatic_Terrestrial_dn_allgene <- Outliers_all_dn_allgene %>% 
  filter(Microhabitat == "Aquatic_Terrestrial")
Outliers_Terrestrial_dn_allgene <- Outliers_all_dn_allgene %>% 
  filter(Microhabitat == "Terrestrial")

Outliers_all_ds_allgene <- outliers_ds_species_allgene %>%
  left_join(.,meta_turtles, by = c("species" = "ID")) %>% 
  remove_missing(.)

Outliers_Marine_ds_allgene <- Outliers_all_ds_allgene %>% 
  filter(Microhabitat == "Marine")
Outliers_Aquatic_ds_allgene <- Outliers_all_ds_allgene %>% 
  filter(Microhabitat == "Aquatic")
Outliers_Aquatic_Terrestrial_ds_allgene <- Outliers_all_ds_allgene %>% 
  filter(Microhabitat == "Aquatic_Terrestrial")
Outliers_Terrestrial_ds_allgene <- Outliers_all_ds_allgene %>% 
  filter(Microhabitat == "Terrestrial")

Outliers_all_dnds_allgene <- outliers_dnds_species_allgene %>%
  left_join(.,meta_turtles, by = c("species" = "ID")) %>% 
  remove_missing(.)

Outliers_Marine_dnds_allgene <- Outliers_all_dnds_allgene %>% 
  filter(Microhabitat == "Marine")
Outliers_Aquatic_dnds_allgene <- Outliers_all_dnds_allgene %>% 
  filter(Microhabitat == "Aquatic")
Outliers_Aquatic_Terrestrial_dnds_allgene <- Outliers_all_dnds_allgene %>% 
  filter(Microhabitat == "Aquatic_Terrestrial")
Outliers_Terrestrial_dnds_allgene <- Outliers_all_dnds_allgene %>% 
  filter(Microhabitat == "Terrestrial")

Outliers_all_pure_allgene <- outliers_pure_species_allgene %>%
  left_join(.,meta_turtles, by = c("species" = "ID")) %>% 
  remove_missing(.)

Outliers_Marine_pure_allgene <- Outliers_all_pure_allgene %>% 
  filter(Microhabitat == "Marine")
Outliers_Aquatic_pure_allgene <- Outliers_all_pure_allgene %>% 
  filter(Microhabitat == "Aquatic")
Outliers_Aquatic_Terrestrial_pure_allgene <- Outliers_all_pure_allgene %>% 
  filter(Microhabitat == "Aquatic_Terrestrial")
Outliers_Terrestrial_pure_allgene <- Outliers_all_pure_allgene %>% 
  filter(Microhabitat == "Terrestrial")

##############################
### GO ENRICHMENT ANALYSIS ###
##############################

### turtles annotation db
go_all <- read.table("Data/odb10v1_OG_xrefs.tab", 
                     header = F, col.names = c("GID","ontology","GO","Support_number"), sep = "\t")

#### BUSCO DB
go_busco <- go_all %>% filter(GID %in% gene_names)

### subset DB associated with biological process
go_bp_busco <- go_busco %>% 
  filter(ontology == "biological_process",
         str_detect(GO,"GO")) %>% 
  dplyr::select(GO,GID) %>% 
  mutate(goterm = Term(GO))
go_bp_name_busco <- go_bp_busco %>% dplyr::select(GO,GID)
go_bp_term_busco <- go_bp_busco %>% dplyr::select(GO,goterm)

### subset DB associated with cellular component
go_cc_busco <- go_busco %>% 
  filter(ontology == "cellular_component",
         str_detect(GO,"GO")) %>% 
  dplyr::select(GO,GID) %>% 
  mutate(goterm = Term(GO))
go_cc_name_busco <- go_cc_busco %>% dplyr::select(GO,GID)
go_cc_term_busco <- go_cc_busco %>% dplyr::select(GO,goterm)


### subset DB associated with molecular function
go_mf_busco <- go_busco %>% 
  filter(ontology == "molecular_function",
         str_detect(GO,"GO")) %>% 
  dplyr::select(GO,GID) %>% 
  mutate(goterm = Term(GO))
go_mf_name_busco <- go_mf_busco %>% dplyr::select(GO,GID)
go_mf_term_busco <- go_mf_busco %>% dplyr::select(GO,goterm)

### DB with all ontologies

go_all_name_busco <- go_busco %>% filter(str_detect(GO,"GO")) %>% dplyr::select(GO,GID)
go_all_term_busco <- go_busco %>% filter(str_detect(GO,"GO")) %>% dplyr::select(GO,GID) %>% 
  mutate(goterm = Term(GO)) %>% dplyr::select(GO,goterm)


##############################
### GO GENE TREES ANALYSIS ###
##############################

enrich_results_busco <- function(df){
  bp <- enricher(gene = unique(df$gene),
                 pAdjustMethod = "fdr",
                 TERM2GENE = go_bp_name_busco,
                 TERM2NAME = go_bp_term_busco)
  cc <- enricher(gene = unique(df$gene),
                 pAdjustMethod = "fdr",
                 TERM2GENE = go_cc_name_busco,
                 TERM2NAME = go_cc_term_busco)
  mf <- enricher(gene = unique(df$gene),
                 pAdjustMethod = "fdr",
                 TERM2GENE = go_mf_name_busco,
                 TERM2NAME = go_mf_term_busco)
  # Add an 'ontology' column to each data frame to indicate the GO category
  bp_df <- as.data.frame(bp) %>% mutate(ontology = "BP")
  cc_df <- as.data.frame(cc) %>% mutate(ontology = "CC")
  mf_df <- as.data.frame(mf) %>% mutate(ontology = "MF")
  
  # Combine the data frames into one
  combined_result <- bind_rows(bp_df, cc_df, mf_df)
  
  return(list(bp=bp, cc=cc, mf=mf, combined_result=combined_result))
  
}

go_genes_dn <- enrich_results_busco(outliers_dn)
go_all_genes_dn <- enrich_results_busco(Outliers_all_dn)
go_Marine_dn <- enrich_results_busco(Outliers_Marine_dn)
go_Aquatic_dn <- enrich_results_busco(Outliers_Aquatic_dn)
go_Aquatic_Terrestrial_dn <- enrich_results_busco(Outliers_Aquatic_Terrestrial_dn)
go_Terrestrial_dn <- enrich_results_busco(Outliers_Terrestrial_dn)

go_genes_ds <- enrich_results_busco(outliers_ds)
go_all_genes_ds <- enrich_results_busco(Outliers_all_ds)
go_Marine_ds <- enrich_results_busco(Outliers_Marine_ds)
go_Aquatic_ds <- enrich_results_busco(Outliers_Aquatic_ds)
go_Aquatic_Terrestrial_ds <- enrich_results_busco(Outliers_Aquatic_Terrestrial_ds)
go_Terrestrial_ds <- enrich_results_busco(Outliers_Terrestrial_ds)

go_genes_dnds <- enrich_results_busco(outliers_dnds)
go_all_genes_dnds <- enrich_results_busco(Outliers_all_dnds)
go_Marine_dnds <- enrich_results_busco(Outliers_Marine_dnds)
go_Aquatic_dnds <- enrich_results_busco(Outliers_Aquatic_dnds)
go_Aquatic_Terrestrial_dnds <- enrich_results_busco(Outliers_Aquatic_Terrestrial_dnds)
go_Terrestrial_dnds <- enrich_results_busco(Outliers_Terrestrial_dnds)

go_genes_pure <- enrich_results_busco(outliers_pure)
go_all_genes_pure <- enrich_results_busco(Outliers_all_pure)
go_Marine_pure <- enrich_results_busco(Outliers_Marine_pure)
go_Aquatic_pure <- enrich_results_busco(Outliers_Aquatic_pure)
go_Aquatic_Terrestrial_pure <- enrich_results_busco(Outliers_Aquatic_Terrestrial_pure)
go_Terrestrial_pure <- enrich_results_busco(Outliers_Terrestrial_pure)

# ALLGENES ANALYSES
go_all_genes_dn_allgene <- enrich_results_busco(Outliers_all_dn_allgene)
go_Marine_dn_allgene <- enrich_results_busco(Outliers_Marine_dn_allgene)
go_Aquatic_dn_allgene <- enrich_results_busco(Outliers_Aquatic_dn_allgene)
go_Aquatic_Terrestrial_dn_allgene <- enrich_results_busco(Outliers_Aquatic_Terrestrial_dn_allgene)
go_Terrestrial_dn_allgene <- enrich_results_busco(Outliers_Terrestrial_dn_allgene)

go_all_genes_ds_allgene <- enrich_results_busco(Outliers_all_ds_allgene)
go_Marine_ds_allgene <- enrich_results_busco(Outliers_Marine_ds_allgene)
go_Aquatic_ds_allgene <- enrich_results_busco(Outliers_Aquatic_ds_allgene)
go_Aquatic_Terrestrial_ds_allgene <- enrich_results_busco(Outliers_Aquatic_Terrestrial_ds_allgene)
go_Terrestrial_ds_allgene <- enrich_results_busco(Outliers_Terrestrial_ds_allgene)

go_all_genes_dnds_allgene <- enrich_results_busco(Outliers_all_dnds_allgene)
go_Marine_dnds_allgene <- enrich_results_busco(Outliers_Marine_dnds_allgene)
go_Aquatic_dnds_allgene <- enrich_results_busco(Outliers_Aquatic_dnds_allgene)
go_Aquatic_Terrestrial_dnds_allgene <- enrich_results_busco(Outliers_Aquatic_Terrestrial_dnds_allgene)
go_Terrestrial_dnds_allgene <- enrich_results_busco(Outliers_Terrestrial_dnds_allgene)

go_all_genes_pure_allgene <- enrich_results_busco(Outliers_all_pure_allgene)
go_Marine_pure_allgene <- enrich_results_busco(Outliers_Marine_pure_allgene)
go_Aquatic_pure_allgene <- enrich_results_busco(Outliers_Aquatic_pure_allgene)
go_Aquatic_Terrestrial_pure_allgene <- enrich_results_busco(Outliers_Aquatic_Terrestrial_pure_allgene)
go_Terrestrial_pure_allgene <- enrich_results_busco(Outliers_Terrestrial_pure_allgene)

ggplot(go_genes_dn[[4]], aes(x=geneID, y=Description, fill=pvalue)) +
  scale_fill_gradientn(colours=colours_classes[2:3]) +
  theme_minimal() + 
  geom_tile()
ggsave("Results/dn_tile.pdf", width = 8, height = 6)

ggplot(go_genes_dn[[4]], aes(x=GeneRatio, y=Description, color=pvalue, size=Count)) + 
  geom_point() +
  scale_color_gradientn(colours=colours_classes[2:3]) +
  theme_minimal() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis dN")
ggsave("Results/dn_point.pdf", width = 8, height = 6)


ggplot(go_genes_dnds[[4]], aes(x=geneID, y=Description, fill=pvalue)) +
  scale_fill_gradientn(colours=colours_classes[2:3]) +
  theme_minimal() + 
  geom_tile()
ggsave("Results/dnds_tile.pdf", width = 8, height = 6)

ggplot(go_genes_dnds[[4]], aes(x=GeneRatio, y=Description, color=pvalue, size=Count)) + 
  geom_point() +
  scale_color_gradientn(colours=colours_classes[2:3]) +
  theme_minimal() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis dN/dS")
ggsave("Results/dnds_point.pdf", width = 8, height = 6)

ggplot(go_all_genes_dn[[4]], aes(x=geneID, y=Description, fill=pvalue)) +
  scale_fill_gradientn(colours=colours_classes[2:3]) +
  theme_minimal() + 
  geom_tile()
ggsave("Results/dn_species_tile.pdf", width = 8, height = 6)

ggplot(go_all_genes_dn[[4]], aes(x=GeneRatio, y=Description, color=pvalue, size=Count)) + 
  geom_point() +
  scale_color_gradientn(colours=colours_classes[2:3]) +
  theme_minimal() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis species dN")
ggsave("Results/dn_species_point.pdf", width = 8, height = 6)

ggplot(go_Marine_dnds[[4]], aes(x=geneID, y=Description, fill=pvalue)) +
  scale_fill_gradientn(colours=colours_classes[2:3]) +
  theme_minimal() + 
  geom_tile()
ggsave("Results/dnds_marine_tile.pdf", width = 8, height = 6)

ggplot(go_Marine_dnds[[4]], aes(x=GeneRatio, y=Description, color=pvalue, size=Count)) + 
  geom_point() +
  scale_color_gradientn(colours=colours_classes[2:3]) +
  theme_minimal() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis Marine dN/dS")
ggsave("Results/dnds_marine_point.pdf", width = 8, height = 6)

ggplot(go_Aquatic_dn[[4]], aes(x=geneID, y=Description, fill=pvalue)) +
  scale_fill_gradientn(colours=colours_classes[2:3]) +
  theme_minimal() + 
  geom_tile()
ggsave("Results/dn_aquatic_tile.pdf", width = 8, height = 6)

ggplot(go_Aquatic_dn[[4]], aes(x=GeneRatio, y=Description, color=pvalue, size=Count)) + 
  geom_point() +
  scale_color_gradientn(colours=colours_classes[2:3]) +
  theme_minimal() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis Aquatic dN")
ggsave("Results/dn_aquatic_point.pdf", width = 8, height = 6)





bp_enrich2 <- pairwise_termsim(go_genes_dn[[2]])
treeplot(bp_enrich2)

bp_enrich2 <- pairwise_termsim(go_genes_dnds[[2]])
treeplot(bp_enrich2)


#############
### PLOTS ###
#############

#### SPECIES TREE PLOT
colours_classes5 <- fish(n=5,option="Balistoides_conspicillum", end=0.95, 
                         begin=0.3,direction=-1)
colours_classes4 <- colours_classes5[1:4]
c("#0F3D5CFF", "#4C98B8FF", "#7EA77DFF", "#9DB327FF", "#DEE100FF")

colour_marine <- colours_classes4[1]
colour_aquatic <- colours_classes4[2]
colour_aquatic_terrestrial <- colours_classes4[3]
colour_terrestrial <- colours_classes4[4]

meta_turtles <- meta_turtles %>% filter(Microhabitat != "Outgroup")

# load species tree for plotting
species_tree_plot <- read.tree("Data/species_tree_plot.nwk")
# rename tips from ID to Name (nicer to read)
species_tree_plot$tip.label <- meta_turtles$Species[
  match(species_tree_plot$tip.label, meta_turtles$ID)]
species_tree_plot <- drop.tip(species_tree_plot, 
                              setdiff(species_tree_plot$tip.label, meta_turtles$Species))

group_habitat <- list(Marine=meta_turtles %>% 
                        filter(Habitat_factor=="Marine") %>% 
                        dplyr::select(Species) %>% .$Species,
                      Aquatic=meta_turtles %>% 
                        filter(Habitat_factor=="Aquatic") %>% 
                        dplyr::select(Species) %>% .$Species,
                      Aquatic_Terrestrial=meta_turtles %>% 
                        filter(Habitat_factor=="Aquatic_Terrestrial") %>% 
                        dplyr::select(Species) %>% .$Species,
                      Terrestrial=meta_turtles %>% 
                        filter(Habitat_factor=="Terrestrial") %>% 
                        dplyr::select(Species) %>% .$Species)

species_tree_plot <- groupOTU(species_tree_plot, group_habitat)

# plot for our tree
plot_tree <- ggtree::ggtree(species_tree_plot, aes(color=group)) + ggtree::xlim_tree(13)
plot_tree <- plot_tree %<+% meta_turtles +
  ggtree::geom_tippoint(aes(color=Habitat_factor)) + 
  ggtree::geom_tiplab(size=3, offset=0.5, fontface = "italic") + 
  scale_color_manual("Primary lifestyle", values=colours_classes4,
                     breaks = c("Marine", "Aquatic", "Aquatic_Terrestrial", "Terrestrial")) +
  theme_tree2() +
  vexpand(0.01, direction = -1)
plot_tree

ggsave("Results/species_tree.pdf", width = 8, height = 5)



#### OVERLAP BETWEEN HABITATS

# DN
ggvenn(
  list(Marine = unique(Outliers_Marine_dn$gene), 
       Aquatic = unique(Outliers_Aquatic_dn$gene), 
       Aquatic_Terrestrial = unique(Outliers_Aquatic_Terrestrial_dn$gene),
       Terrestrial = unique(Outliers_Terrestrial_dn$gene)), 
  fill_color = colours_classes4,
  stroke_size = 0, set_name_size = 4, show_percentage = F
)
ggsave('Results/venn_habitat_overlap_dn.pdf', width = 8, height = 8)

# DS
ggvenn(
  list(Marine = unique(Outliers_Marine_ds$gene), 
       Aquatic = unique(Outliers_Aquatic_ds$gene), 
       Aquatic_Terrestrial = unique(Outliers_Aquatic_Terrestrial_ds$gene),
       Terrestrial = unique(Outliers_Terrestrial_ds$gene)), 
  fill_color = colours_classes4,
  stroke_size = 0, set_name_size = 4, show_percentage = F
)
ggsave('Results/venn_habitat_overlap_ds.pdf', width = 8, height = 8)

# DNDS
ggvenn(
  list(Marine = unique(Outliers_Marine_dnds$gene), 
       Aquatic = unique(Outliers_Aquatic_dnds$gene), 
       Aquatic_Terrestrial = unique(Outliers_Aquatic_Terrestrial_dnds$gene),
       Terrestrial = unique(Outliers_Terrestrial_dnds$gene)), 
  fill_color = colours_classes4,
  stroke_size = 0, set_name_size = 4, show_percentage = F
)
ggsave('Results/venn_habitat_overlap_dnds.pdf', width = 8, height = 8)

# PURE
ggvenn(
  list(Marine = unique(Outliers_Marine_pure$gene), 
       Aquatic = unique(Outliers_Aquatic_pure$gene), 
       Aquatic_Terrestrial = unique(Outliers_Aquatic_Terrestrial_pure$gene),
       Terrestrial = unique(Outliers_Terrestrial_pure$gene)), 
  fill_color = colours_classes4,
  stroke_size = 0, set_name_size = 4, show_percentage = F
)
ggsave('Results/venn_habitat_overlap_pure.pdf', width = 8, height = 8)


#### OVERLAP BETWEEN SPECIES
species_list <- meta_turtles %>% filter(Microhabitat != "Outgroup")
species_list <- species_list$ID

df_species_id <- meta_turtles %>% dplyr::select(c("Species", "ID"))

species_ordered_list <- c("Malaclemys terrapin",
                          "Chrysemys picta",
                          "Terrapene mexicana",
                          "Platysternon megacephalum",
                          "Cuora amboinensis",
                          "Cuora mccordi",
                          "Gopherus agassizii",
                          "Chelonoidis abingdonii",
                          "Dermatemys mawii",
                          "Chelydra serpentina",
                          "Dermochelys coriacea",
                          "Chelonia mydas",
                          "Pelodiscus sinensis",
                          "Carettochelys insculpta",
                          "Emydura subglobosa",
                          "Mesoclemmys tuberculata",
                          "Podocnemis expansa",
                          "Pelusios castaneus")

habitat_colours_list <- c(colour_aquatic,
                          colour_aquatic,
                          colour_aquatic_terrestrial,
                          colour_aquatic,
                          colour_terrestrial,
                          colour_aquatic_terrestrial,
                          colour_terrestrial,
                          colour_terrestrial,
                          colour_aquatic,
                          colour_aquatic,
                          colour_marine,
                          colour_marine,
                          colour_aquatic,
                          colour_aquatic,
                          colour_aquatic,
                          colour_aquatic,
                          colour_aquatic,
                          colour_aquatic)

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
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#C7AF5A", "#BB4430", "#6E291D")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, 
                                   colour = habitat_colours_list),
        axis.text.y = element_text(colour = habitat_colours_list[rev(1:length(habitat_colours_list))]))
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
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#C7AF5A", "#BB4430", "#6E291D")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, 
                                   colour = habitat_colours_list),
        axis.text.y = element_text(colour = habitat_colours_list[rev(1:length(habitat_colours_list))]))
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
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#C7AF5A", "#BB4430", "#6E291D")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, 
                                   colour = habitat_colours_list),
        axis.text.y = element_text(colour = habitat_colours_list[rev(1:length(habitat_colours_list))]))
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
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#C7AF5A", "#BB4430", "#6E291D")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, 
                                   colour = habitat_colours_list),
        axis.text.y = element_text(colour = habitat_colours_list[rev(1:length(habitat_colours_list))]))
ggsave("Results/heatmap_pure.pdf", width = 8, height = 5)
