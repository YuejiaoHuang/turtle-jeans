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
distances_meangene <- distances_extraction(comb_imputed_meangene)
outliergene_chi95_meangene <- detect_outliers_and_extract_quantiles(matrices_imputed_meangene,
                                                                distances = distances_meangene,
                                                                quantile_threshold = 0.95)

distances_meanlist <- distances_extraction(comb_imputed_meanlist)
outliergene_chi95_meanlist<-detect_outliers_and_extract_quantiles(matrices_imputed_meanlist,
                                                                  distance = distances_meanlist,
                                                                  quantile_threshold = 0.95)




# get outlier species - v2 ------------------------------------------------

mat_name <- names(outliergene_chi95_meangene$matrices)[1]
mat_list <- outliergene_chi95_meanlist$matrices[1]
mat <- outliergene_chi95_meangene$matrices[[mat_name]]

# use the entire thing, not just upper triangle 
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









