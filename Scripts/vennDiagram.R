library(VennDiagram)

# gene vs list 
meangene <- read.csv('Results/outliergene_q95_meangene.csv',
                     col.names = 'gene')
meanlist <- read.csv('Results/outliergene_q95_meanlist.csv',
                     col.names='gene')

# Prepare a palette of 3 colors with R colorbrewer:
cols <- c('#2274a5','#E6AF2E')

# Chart
ggvenn(data = list('Gene' = meangene$gene,
                   'List' = meanlist$gene),
       fill_color = cols,
       show_percentage = F)+
  scale_y_continuous(expand = expansion(mult = .3))


ggsave('Results/venn_genelist_quant.pdf',last_plot(),
       width = 4,
       height = 4)


#### GENE
# chi vs quan
chi <- read.csv('Results/outliergene_chi95_meangene.csv',
                     col.names = c('gene'))
quan <- read.csv('Results/outliergene_q95_meangene.csv',
                     col.names = c('gene'))

# Prepare a palette of 3 colors with R colorbrewer:
cols <- c('#A0C4D8','#083047')

# Chart

ggvenn(data = list('ChiSq' = chi$gene,
                   'Quantile' = quan$gene),
       fill_color = cols,
       show_percentage = F)+
  scale_y_continuous(expand = expansion(mult = .3))

ggsave('Results/venn_gene_chiquant.pdf',last_plot())

#### List
# chi vs quan
chi <- read.csv('Results/outliergene_chi95_meanlist.csv',
                col.names = c('gene'))
quan <- read.csv('Results/outliergene_q95_meanlist.csv',
                 col.names = c('gene'))

# Prepare a palette of 3 colors with R colorbrewer:
cols <- c('#DCCAA1','#6E4F08')

# Chart

ggvenn(data = list('ChiSq' = chi$gene,
                   'Quantile' = quan$gene),
       fill_color = cols,
       show_percentage = F)+
  scale_y_continuous(expand = expansion(mult = .3))

ggsave('Results/venn_list_chiquant.pdf',last_plot())





