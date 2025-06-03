library(ggvenn)
library(fishualize)
library(svglite)

# library(showtext)
# font_add(family = "Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
# showtext_auto()

library(systemfonts)
# register_font(
#   name = "Arial",
#   plain = "/System/Library/Fonts/Supplemental/Arial.ttf"
# )

library(extrafont)
# font_import(paths = "/System/Library/Fonts/Supplemental", prompt = FALSE)
# loadfonts(device = "pdf")



setwd("/Users/jule/Desktop/turtle-jeans")

colours_classes5 <- fish(n=5,option="Balistoides_conspicillum", end=0.95, 
                         begin=0.3,direction=-1)

colour_ds <- colours_classes5[1]
colour_dn <- colours_classes5[2]
colour_dnds <- colours_classes5[4]
colour_pairwise_dists <- colours_classes5[5]

c("#0F3D5CFF", "#4C98B8FF", "#9DB327FF", "#DEE100FF")

# or would this be nicer
c("#17638DFF", "#4694B8FF", "#9DB327FF", "#DEE100FF")

# genes
gene_dnds <- read.csv('Results/outliers_genes_dnds_q95.csv',
                      col.names = c('n', 'gene'))
gene_dn <- read.csv('Results/outliers_genes_dn_q95.csv',
                    col.names = c('n', 'gene'))
gene_ds <- read.csv('Results/outliers_genes_ds_q95.csv',
                    col.names = c('n', 'gene'))
gene_pairwise_dists <- read.csv('Results/outliers_genes_pure_q95.csv',
                                col.names = c('n', 'gene'))

# create labels
labels_gene <- c(
  ds = paste0("dS (", length(unique(gene_ds$gene)), ")"),
  dn = paste0("dN (", length(unique(gene_dn$gene)), ")"),
  dnds = paste0("dN/dS (", length(unique(gene_dnds$gene)), ")"),
  pairwise_dists = paste0("Raw distances (", length(unique(gene_pairwise_dists$gene)), ")")
)

# match labels with lists
venn_data_gene <- setNames(
  list(
    unique(gene_ds$gene),
    unique(gene_dn$gene),
    unique(gene_dnds$gene),
    unique(gene_pairwise_dists$gene)
  ),
  labels_gene
)

venn_gene <- ggvenn(venn_data_gene, 
                    fill_color = c(colour_ds, colour_dn, 
                                   colour_dnds, colour_pairwise_dists),
                    stroke_size = 0, set_name_size = 5, show_percentage = F, stroke_color = "grey",
                    set_name_color = c(colour_ds, colour_dn, 
                                       colour_dnds, colour_pairwise_dists), text_color = "black"
)

venn_gene <- venn_gene + theme(text = element_text(family = "Arial"),
                               plot.margin = margin(t = 0, r = 60, b = 0, l = 0)) +
  coord_cartesian(clip = "off")
venn_gene

svglite('Results/venn_outliergenes.svg', width = 11, height = 8)
print(venn_gene)
dev.off()

ggsave('Results/venn_outliergenes.pdf',
       width = 11,
       height = 8)

ggsave('Results/venn_outliergenes.png',
       width = 11,
       height = 8)


### REPEAT WITH SPECIES
# genes
species_dnds <- read.csv('Results/outliers_species_dnds_q95.csv',
                         col.names = c('n', 'species', 'gene'))
species_dn <- read.csv('Results/outliers_species_dn_q95.csv',
                       col.names = c('n', 'species', 'gene'))
species_ds <- read.csv('Results/outliers_species_ds_q95.csv',
                       col.names = c('n', 'species', 'gene'))
species_pairwise_dists <- read.csv('Results/outliers_species_pure_q95.csv',
                                   col.names = c('n', 'species', 'gene'))

# create labels
labels_species <- c(
  ds = paste0("dS (", length(unique(species_ds$gene)), ")"),
  dn = paste0("dN (", length(unique(species_dn$gene)), ")"),
  dnds = paste0("dN/dS (", length(unique(species_dnds$gene)), ")"),
  pairwise_dists = paste0("Raw distances (", length(unique(species_pairwise_dists$gene)), ")")
)

# match labels with lists
venn_data_species <- setNames(
  list(
    unique(species_ds$gene),
    unique(species_dn$gene),
    unique(species_dnds$gene),
    unique(species_pairwise_dists$gene)
  ),
  labels_species
)

venn_species <- ggvenn(venn_data_species, 
       fill_color = c(colour_ds, colour_dn, 
                      colour_dnds, colour_pairwise_dists),
       stroke_size = 0, set_name_size = 5, show_percentage = F, stroke_color = "grey",
       set_name_color = c(colour_ds, colour_dn, 
                          colour_dnds, colour_pairwise_dists), text_color = "black"
)
venn_species <- venn_species + 
  theme(text = element_text(family = "Arial"),
                                     plot.margin = margin(t = 0, r = 60, b = 0, l = 0)) +
  coord_cartesian(clip = "off")
venn_species

svglite('Results/venn_outlierspecies.svg', width = 11, height = 8)
print(venn_species)
dev.off()

ggsave('Results/venn_outlierspecies.pdf',last_plot(),
       width = 11,
       height = 8)

ggsave('Results/venn_outlierspecies.png',last_plot(),
       width = 11,
       height = 8)



