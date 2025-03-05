library('ggvenn')

# genes
gene_dnds <- read.csv('Results/outliers_genes_dnds_q95.csv',
                     col.names = 'gene')
gene_dn <- read.csv('Results/outliers_genes_dn_q95.csv',
                      col.names = 'gene')
gene_ds <- read.csv('Results/outliers_genes_ds_q95.csv',
                      col.names = 'gene')
gene_pure <- read.csv('Results/outliers_genes_pure_q95.csv',
                      col.names = 'gene')


cols <- c('purple','red','blue','green')

# Chart
ggvenn(data = list('dN/dS' = gene_dnds$gene,
                   'dN' = gene_dn$gene,
                   'dS' = gene_ds$gene,
                   'PD' = gene_pure$gene),
       fill_color = cols,
       labs(title = 'Test'),
       show_percentage = F)


ggsave('Results/venn_outliergenes.pdf',last_plot(),
       width = 8,
       height = 8)

### REPEAT WITH SPECIES
# genes
species_dnds <- read.csv('Results/outliers_species_dnds_q95.csv',
                      col.names = 'gene')
species_dn <- read.csv('Results/outliers_species_dn_q95.csv',
                    col.names = 'gene')
species_ds <- read.csv('Results/outliers_species_ds_q95.csv',
                    col.names = 'gene')
species_pure <- read.csv('Results/outliers_species_pure_q95.csv',
                      col.names = 'gene')

# Prepare a palette of 3 colors with R colorbrewer:
cols <- c('purple','red','blue','green')

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




