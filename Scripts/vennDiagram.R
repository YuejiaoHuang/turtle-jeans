library('ggvenn')

setwd("/Users/jule/Desktop/turtle-jeans")

colours_classes5 <- fish(n=5,option="Balistoides_conspicillum", end=0.95, 
                         begin=0.3,direction=-1)

colour_ds <- colours_classes5[1]
colour_dn <- colours_classes5[2]
colour_dnds <- colours_classes5[4]
colour_pure <- colours_classes5[5]

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
gene_pure <- read.csv('Results/outliers_genes_pure_q95.csv',
                      col.names = c('n', 'gene'))

# create labels
labels_gene <- c(
  ds = paste0("dS (", length(unique(gene_ds$gene)), ")"),
  dn = paste0("dN (", length(unique(gene_dn$gene)), ")"),
  dnds = paste0("dN/dS (", length(unique(gene_dnds$gene)), ")"),
  Pure = paste0("Pure (", length(unique(gene_pure$gene)), ")")
)

# match labels with lists
venn_data_gene <- setNames(
  list(
    unique(gene_ds$gene),
    unique(gene_dn$gene),
    unique(gene_dnds$gene),
    unique(gene_pure$gene)
  ),
  labels_gene
)

ggvenn(venn_data_gene, 
       fill_color = c(colour_ds, colour_dn, 
                      colour_dnds, colour_pure),
       stroke_size = 0, set_name_size = 5, show_percentage = F, stroke_color = "grey",
       set_name_color = c(colour_ds, colour_dn, 
                          colour_dnds, colour_pure), text_color = "grey35"
)

ggsave('Results/venn_outliergenes.pdf',last_plot(),
       width = 8,
       height = 8)


### REPEAT WITH SPECIES
# genes
species_dnds <- read.csv('Results/outliers_species_dnds_q95.csv',
                      col.names = c('n', 'species', 'gene'))
species_dn <- read.csv('Results/outliers_species_dn_q95.csv',
                       col.names = c('n', 'species', 'gene'))
species_ds <- read.csv('Results/outliers_species_ds_q95.csv',
                       col.names = c('n', 'species', 'gene'))
species_pure <- read.csv('Results/outliers_species_pure_q95.csv',
                         col.names = c('n', 'species', 'gene'))

# create labels
labels_species <- c(
  ds = paste0("dS (", length(unique(species_ds$gene)), ")"),
  dn = paste0("dN (", length(unique(species_dn$gene)), ")"),
  dnds = paste0("dN/dS (", length(unique(species_dnds$gene)), ")"),
  Pure = paste0("Pure (", length(unique(species_pure$gene)), ")")
)

# match labels with lists
venn_data_species <- setNames(
  list(
    unique(species_ds$gene),
    unique(species_dn$gene),
    unique(species_dnds$gene),
    unique(species_pure$gene)
  ),
  labels_species
)

ggvenn(venn_data_species, 
       fill_color = c(colour_ds, colour_dn, 
                      colour_dnds, colour_pure),
       stroke_size = 0, set_name_size = 5, show_percentage = F, stroke_color = "grey",
       set_name_color = c(colour_ds, colour_dn, 
                          colour_dnds, colour_pure), text_color = "grey35"
)
# Chart
ggvenn(data = list('dNdS' = species_dnds$gene,
                   'dN' = species_dn$gene,
                   'dS' = species_ds$gene,
                   'Pure' = species_pure$gene),
       fill_color = cols,
       show_percentage = T)+
  scale_y_continuous(expand = expansion(mult = .3))


ggsave('Results/venn_outlierspecies.pdf',last_plot(),
       width = 8,
       height = 8)


