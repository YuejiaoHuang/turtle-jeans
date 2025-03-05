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

# Attempt 2 - outlier genes -----------------------------------------------

# Required Libraries
library(MASS)  # For Mahalanobis distance
library(stats) # For Chi-squared distribution

### LOAD MATRIX OF INTEREST 
source('Scripts/helperFunctions.R')

#dN_dS_matrices[[1]]
dNdS_checked <- matrixCheck(dN_dS_matrices_nonzero)
dNdS_checked[['121468at32523']]

dN_checked <- matrixCheck(dN_nonzero)
dN_checked[['121468at32523']]

dS_checked <- matrixCheck(dS_nonzero)
dS_checked[['121468at32523']]

pure_checked <- matrixCheck(matrices)
pure_checked[['121468at32523']]
#impute from average across matrix
source("Scripts/impMean_matrixaverage.R")
dNdS_impute <- impMean_matrices(dNdS_checked)
dN_impute <- impMean_matrices(dN_checked)
dS_impute <- impMean_matrices(dS_checked)
pure_impute <- impMean_matrices(pure_checked)

dNdS_impute[['121468at32523']]
dN_impute[['121468at32523']]
dS_impute[['121468at32523']]
pure_impute[['121468at32523']]

# make combined dataframe
source("Scripts/stitchMatricesToDataFrame.R")
dNdS_2d <- stitchMatricesToDataFrame(dNdS_impute)
dN_2d <- stitchMatricesToDataFrame(dN_impute)
dS_2d <- stitchMatricesToDataFrame(dS_impute)
pure_2d <- stitchMatricesToDataFrame(pure_impute)

dNdS_log <- log(dNdS_2d)
dN_log <- log(dN_2d)
dS_log <- log(dS_2d)
pure_log <- log(pure_2d)

#detect outliers
source('Scripts/detectOutliers.R')

# Detect outliers and extract quantiles and chi-square for meangene
dNdS_outliers <- detect_outliers_and_extract_quantiles(dNdS_impute, dNdS_log, quantile_threshold = 0.95)
dN_outliers <- detect_outliers_and_extract_quantiles(dN_impute, dN_log, quantile_threshold = 0.95)
dS_outliers <- detect_outliers_and_extract_quantiles(dS_impute, dS_log, quantile_threshold = 0.95)
pure_outliers <- detect_outliers_and_extract_quantiles(pure_impute, pure_log, quantile_threshold = 0.95)

text <- dNdS_outliers$indices
write.csv(text,'Results/outliers_genes_dnds_q95.csv')

text <- dN_outliers$indices
write.csv(text,'Results/outliers_genes_dn_q95.csv')

text <- dS_outliers$indices
write.csv(text,'Results/outliers_genes_ds_q95.csv')

text <- pure_outliers$indices
write.csv(text,'Results/outliers_genes_pure_q95.csv')


# outlier species - just outlier genes ------------------------------------------------
# keep entire thing


outlierSpeciesTranspose <- function(matrices_list) {
  # Extract the gene names and the list of matrices
  gene_names <- matrices_list[[1]]
  matrices <- matrices_list[[2]]
  
  # Initialize an empty list to store the renamed matrices
  matrices_renamed <- list()
  
  # Loop through the matrices and rename rows based on gene names
  for (i in 1:length(matrices)) {
    # Get the current matrix
    mat <- matrices[[i]]
    
    # Assign the new row names by appending the gene name
    rownames(mat) <- paste0(rownames(mat), '_', names(matrices)[i])
    
    # Add the renamed matrix to the list
    matrices_renamed[[i]] <- mat
  }
  
  # Combine the matrices by rows
  df <- do.call(rbind, matrices_renamed)
  
  return(df)
}

dNdS_df <- outlierSpeciesTranspose(dNdS_outliers)
dN_df <- outlierSpeciesTranspose(dN_outliers)
dS_df <- outlierSpeciesTranspose(dS_outliers)
pure_df <- outlierSpeciesTranspose(pure_outliers)

dNdS_outliers_species <- detect_outliers_and_extract_chisq(dNdS_impute, dNdS_df, conf = 0.95)
dN_outliers_species <- detect_outliers_and_extract_chisq(dN_impute, dN_df, conf = 0.95)
dS_outliers_species <- detect_outliers_and_extract_chisq(dS_impute, dS_df, conf = 0.95)
pure_outliers_species <- detect_outliers_and_extract_chisq(pure_impute, pure_df, conf = 0.95)

text <- dNdS_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(Gene = x[1], Species = x[2])))
write.csv(df,'Results/outliers_species_dnds_q95.csv')

text <- dN_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(Gene = x[1], Species = x[2])))
write.csv(df,'Results/outliers_species_dn_q95.csv')

text <- dS_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(Gene = x[1], Species = x[2])))
write.csv(df,'Results/outliers_species_ds_q95.csv')

text <- pure_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(Gene = x[1], Species = x[2])))
write.csv(df,'Results/outliers_species_pure_q95.csv')



# outlier species - all genes ---------------------------------------------


# keep entire thing


outlierSpeciesTranspose <- function(matrices) {
  # Extract the gene names and the list of matrices

  # Initialize an empty list to store the renamed matrices
  matrices_renamed <- list()
  
  # Loop through the matrices and rename rows based on gene names
  for (i in 1:length(matrices)) {
    # Get the current matrix
    mat <- matrices[[i]]
    
    # Assign the new row names by appending the gene name
    rownames(mat) <- paste0(rownames(mat), '_', names(matrices)[i])
    
    # Add the renamed matrix to the list
    matrices_renamed[[i]] <- mat
  }
  
  # Combine the matrices by rows
  df <- do.call(rbind, matrices_renamed)
  
  return(df)
}

dNdS_df <- outlierSpeciesTranspose(dNdS_impute)
dN_df <- outlierSpeciesTranspose(dN_impute)
dS_df <- outlierSpeciesTranspose(dS_impute)
pure_df <- outlierSpeciesTranspose(pure_impute)

dNdS_outliers_species <- detect_outliers_and_extract_chisq(dNdS_impute, dNdS_df, conf = 0.95)
dN_outliers_species <- detect_outliers_and_extract_chisq(dN_impute, dN_df, conf = 0.95)
dS_outliers_species <- detect_outliers_and_extract_chisq(dS_impute, dS_df, conf = 0.95)
pure_outliers_species <- detect_outliers_and_extract_chisq(pure_impute, pure_df, conf = 0.95)

text <- dNdS_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(Gene = x[1], Species = x[2])))
write.csv(df,'Results/allgenes_species_dnds_q95.csv')

text <- dN_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(Gene = x[1], Species = x[2])))
write.csv(df,'Results/allgenes_species_dn_q95.csv')

text <- dS_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(Gene = x[1], Species = x[2])))
write.csv(df,'Results/allgenes_species_ds_q95.csv')

text <- pure_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(Gene = x[1], Species = x[2])))
write.csv(df,'Results/allgenes_species_pure_q95.csv')






