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

# species driving ---------------------------------------------------------

# Function to identify outlier species in a matrix
find_outlier_rows <- function(mat,conf = 0.95) {
  mat[mat == 0] <- NA
  tab <- t(apply(mat,1,function(row) {
    c(mean = mean(row, na.rm = TRUE),
      median = median(row, na.rm = TRUE),
      sd = sd(row, na.rm = TRUE),
      max = max(row, na.rm = TRUE),
      min = min(row, na.rm = TRUE))
  }))
  
  # Calculate Mahalanobis distances for rows
  center_rows <- mean(tab, na.rm = TRUE)
  cov_rows <- cov(tab, use = "pairwise.complete.obs")
  
  mahalanobis_dist_rows <- mahalanobis(tab, center_rows, cov_rows)
  
  df <- ncol(tab)
  threshold <- chisq(conf,df)
  
  # Identify outlier indices
  outliers <- which(mahalanobis_dist_rows > threshold)
  return(outliers)
}

# Apply to flagged matrices
outlier_matrices <- outliers90$indices  # Replace with indices of outlier matrices
outlier_species90 <- lapply(outlier_matrices, function(i) {
  find_outlier_rows(matrices[[i]])
})

# Apply to flagged matrices
outlier_matrices <- outliers95$indices  # Replace with indices of outlier matrices
outlier_species95 <- lapply(outlier_matrices, function(i) {
  find_outlier_rows(matrices[[i]])
})

outlier_matrices <- outliers99$indices  # Replace with indices of outlier matrices
outlier_species99 <- lapply(outlier_matrices, function(i) {
  find_outlier_rows(matrices[[i]])
})

table(outlier_species95)

print("Outlier species for each matrix:")
print(outlier_species99)

# Combine species from all outlier matrices
unique_species <- names(unlist(outlier_species90))
unique90 <- strsplit(unique_species,"\\.")
unique90 <- sapply(unique90,`[`,2)
length(table(unique90))

unique_species <- names(unlist(outlier_species95))
unique95 <- strsplit(unique_species,"\\.")
unique95 <- sapply(unique95,`[`,2)
length(table(unique95))


unique_species <- names(unlist(outlier_species99))
unique99 <- strsplit(unique_species,"\\.")
unique99 <- sapply(unique99,`[`,2)
length(table(unique99))






# Attempt 2 - outlier genes -----------------------------------------------

# Required Libraries
library(MASS)  # For Mahalanobis distance
library(stats) # For Chi-squared distribution

# Generate all combinations of two turtles
combinations <- combn(tips, 2, simplify = FALSE)

# Merge the names alphabetically with an underscore
combinations_vect <- sapply(combinations, function(x) paste(sort(x), collapse = "_"))

# empty dataframe
df <- as.data.frame(matrix(NA, nrow = length(matrices), ncol = length(combinations_vect)))
colnames(df) <- combinations_vect
rownames(df) <- names(matrices)

#impute from average across matrix
source("Scripts/impMean_matrixaverage.R")
matrices_imputed_meangene <- impMean_matrices(matrices)
table(lapply(matrices_imputed_phylter,dim)[[1]])

# impute from phylter
matrices_imputed_meanlist <- impMean(matrices)
table(lapply(matrices_imputed_phylter,dim)[[1]])

# make combined dataframe
source("Scripts/stitchMatricesToDataFrame.R")
comb_imputed_meangene <- stitchMatricesToDataFrame(matrices_imputed_meangene)
comb_imputed_meanlist <- stitchMatricesToDataFrame(matrices_imputed_meanlist)

#detect outliers
source('Scripts/detectOutliers.R')

# Detect outliers and extract quantiles and chi-square for meangene
outliergene_q90_meangene <- detect_outliers_and_extract_quantiles(matrices_imputed_meangene, comb_imputed_meangene, quantile_threshold = 0.90)
outliergene_q95_meangene <- detect_outliers_and_extract_quantiles(matrices_imputed_meangene, comb_imputed_meangene, quantile_threshold = 0.95)
outliergene_q99_meangene <- detect_outliers_and_extract_quantiles(matrices_imputed_meangene, comb_imputed_meangene, quantile_threshold = 0.99)

outliergene_chi90_meangene <- detect_outliers_and_extract_chisq(matrices_imputed_meangene, comb_imputed_meangene, conf = 0.90)
outliergene_chi95_meangene <- detect_outliers_and_extract_chisq(matrices_imputed_meangene, comb_imputed_meangene, conf = 0.95)
outliergene_chi99_meangene <- detect_outliers_and_extract_chisq(matrices_imputed_meangene, comb_imputed_meangene, conf = 0.99)

# Detect outliers and extract quantiles and chi-square for meanlist
outliergene_q90_meanlist <- detect_outliers_and_extract_quantiles(matrices_imputed_meanlist, comb_imputed_meanlist, quantile_threshold = 0.90)
outliergene_q95_meanlist <- detect_outliers_and_extract_quantiles(matrices_imputed_meanlist, comb_imputed_meanlist, quantile_threshold = 0.95)
outliergene_q99_meanlist <- detect_outliers_and_extract_quantiles(matrices_imputed_meanlist, comb_imputed_meanlist, quantile_threshold = 0.99)

outliergene_chi90_meanlist <- detect_outliers_and_extract_chisq(matrices_imputed_meanlist, comb_imputed_meanlist, conf = 0.90)
outliergene_chi95_meanlist <- detect_outliers_and_extract_chisq(matrices_imputed_meanlist, comb_imputed_meanlist, conf = 0.95)
outliergene_chi99_meanlist <- detect_outliers_and_extract_chisq(matrices_imputed_meanlist, comb_imputed_meanlist, conf = 0.99)

# write CSVs 
write.csv(data.frame(Outlier = outliergene_q90_meangene$indices), 'Results/outliergene_q90_meangene.csv', row.names = FALSE)
write.csv(data.frame(Outlier = outliergene_q95_meangene$indices), 'Results/outliergene_q95_meangene.csv', row.names = FALSE)
write.csv(data.frame(Outlier = outliergene_q99_meangene$indices), 'Results/outliergene_q99_meangene.csv', row.names = FALSE)

write.csv(data.frame(Outlier = outliergene_q90_meanlist$indices), 'Results/outliergene_q90_meanlist.csv', row.names = FALSE)
write.csv(data.frame(Outlier = outliergene_q95_meanlist$indices), 'Results/outliergene_q95_meanlist.csv', row.names = FALSE)
write.csv(data.frame(Outlier = outliergene_q99_meanlist$indices), 'Results/outliergene_q99_meanlist.csv', row.names = FALSE)

write.csv(data.frame(Outlier = outliergene_chi90_meangene$indices), 'Results/outliergene_chi90_meangene.csv', row.names = FALSE)
write.csv(data.frame(Outlier = outliergene_chi95_meangene$indices), 'Results/outliergene_chi95_meangene.csv', row.names = FALSE)
write.csv(data.frame(Outlier = outliergene_chi99_meangene$indices), 'Results/outliergene_chi99_meangene.csv', row.names = FALSE)

write.csv(data.frame(Outlier = outliergene_chi90_meanlist$indices), 'Results/outliergene_chi90_meanlist.csv', row.names = FALSE)
write.csv(data.frame(Outlier = outliergene_chi95_meanlist$indices), 'Results/outliergene_chi95_meanlist.csv', row.names = FALSE)
write.csv(data.frame(Outlier = outliergene_chi99_meanlist$indices), 'Results/outliergene_chi99_meanlist.csv', row.names = FALSE)




# get outlier species - v2 ------------------------------------------------

mat_name <- names(outliergene_q95_meanlist$matrices)[1]
mat <- outliergene_q95_meanlist$matrices[[mat_name]]

# just upper triangle 
# all diagonals to zero
# 
get_upper_triangle_values <- function(mat) {
  # Get the indices of the upper triangle values
  upper_indices <- which(upper.tri(mat), arr.ind = TRUE)
  
  # Extract the upper triangle values
  upper_values <- mat[upper_indices]
  
  # Extract the rownames
  rownames <- rownames(mat)[upper_indices[,1]]
  
  # Extract the columnnames
  columnnames <- colnames(mat)[upper_indices[,2]]
  
  list(values = upper_values,rows = rownames, columns = columnnames)
}

vals <- get_upper_triangle_values(mat)
vals$gene_rows <- paste0(vals$rows,'_',mat_name)

long <- data.frame(values = vals$values,
                   rows = vals$gene_rows,
                   columns = vals$columns)

wide <- reshape(
  long,
  idvar='rows',
  timevar='columns',
  direction='wide'
)

# keep entire thing
names(mat)
rownames(mat) <- paste0(rownames(mat),'_',mat_name)

modified_list <- lapply(names(outliergene_chi95_meanlist$matrices), function(matname){
  mat <- outliergene_chi95_meanlist$matrices[[matname]]
  rownames(mat) <- paste(rownames(mat), matname, sep = "_")
  return(mat)
})

test <- do.call(rbind,modified_list)
df <- as.data.frame(test)
head(df,1)

### zeros?

#detect outliers
source('Scripts/detectOutliers.R')
distances_meanlist <- distances_extraction(df)
detect_outliers_and_extract_chisq <- function(distances_meanlist, conf = 0.95) {
  # Calculate the quantile threshold
  degrees_freedom <- ncol(df) 
  threshold <- qchisq(conf, degrees_freedom)
  
  # Identify outlier indices
  outlier_indices <- names(which(distances_meanlist > threshold))
  
  outlierspecies <- strsplit(outlier_indices,'_')
  outliers_df <- as.data.frame(do.call(rbind, outlierspecies))
  colnames(outliers_df) <- c('species','gene')
  return(outliers_df)
}

outlierspecies_q95_phylter_chi95 <- detect_outliers_and_extract_chisq(distances_meanlist,0.95)
write.csv(outlierspecies_q95_phylter_chi95,'Results/outlierspecies_q95_phylter_chi95.csv')
table(outlierspecies_q95_phylter_chi95$species)

length(unique(outlierspecies_q95_phylter_chi95$gene))




