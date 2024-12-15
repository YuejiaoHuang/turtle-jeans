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

# Optional scaling/normalization


# Rowwise -----------------------------------------------------------------


## OUTLIER SPECIES
# for taxon 1: compare row for taxon1 between all genes
# ...

# for taxon 1: take mean distance to other taxa for each gene
# plot all means
# take outliers from that distribution
# -> Outlier genes for taxa 1

all <- data.frame()

for(gene_name in names){
  # print(gene_name)
  matrix <- matrices[[gene_name]]
  sums <- rowSums(matrix)
  means <- sums/(ncol(matrix)-1)
  frame <- data.frame(name = names(means), 
                      mean = means, 
                      gene = gene_name,
                      row.names = NULL)
  all <- rbind(all,frame)
}


library(reshape2)
reshape <- dcast(all, name ~ gene, value.var = 'mean')
all2 <- reshape[,-1]
rownames(all2) <- reshape[,1]

for(i in tips){
  print(i)
  plot(density(all$mean[all$name == i]),main = i)
}

write.csv(all,'Results/rowwise_species_gene_long.csv')
write.csv(all2,'Results/rowwise_species_gene_wide.csv')


# triangle-wise -----------------------------------------------------------


## OUTLIER GENES
# for build distribution for each gene (based on values for all pairwise distances)
#     plot all pariwise distances from upper triangle of matrix
#     fit parameters for that distribution ISSUE: WHICH DISTRIBUTIONS, HOW TO COMPARE DISTRIBUTIONS
# which distribution is different (based on comparing parameters)?




# mahalabonis -------------------------------------------------------------

# Required Libraries
library(MASS)  # For Mahalanobis distance
library(stats) # For Chi-squared distribution

# calculate upper tri instead of matrix
# get upper triangle
matrices_example <- matrices$`100088at32523`
matrices_example_vector<-as.vector(matrices_example[upper.tri(matrices_example)])
mean(matrices_example_vector)

data_frame(val = matrices_example_vector) %>%
  ggplot(., aes(val)) + 
  geom_density()

# think about stats
# include min
# fit distribution / parameter estimates?  

summary_stats <- do.call(rbind, lapply(matrices, function(mat) {
  tri <- as.vector(mat[upper.tri(mat)])
  c(mean = mean(tri, na.rm = TRUE),
    median = median(tri, na.rm = TRUE),
    sd = sd(tri, na.rm = TRUE),
    max = max(tri, na.rm = TRUE),
    min = min(tri, na.rm = TRUE))
}))

# Compute Mahalanobis distances for matrices
center_matrix <- colMeans(summary_stats, na.rm = TRUE)
cov_matrix <- cov(summary_stats, use = "pairwise.complete.obs")
distances <- mahalanobis(summary_stats, center_matrix, cov_matrix)
length(summary_stats)

detect_outliers_and_extract_quantiles <- function(matrix_list, distances, quantile_threshold = 0.95) {
  # Calculate the quantile threshold
  threshold <- quantile(distances, quantile_threshold)
  # Identify outlier indices
  outlier_indices <- which(distances > threshold)
  # Extract outlier matrices
  outlier_matrices <- matrix_list[outlier_indices]
  # Return a list with indices and matrices
  list(indices = outlier_indices, matrices = outlier_matrices)
}

detect_outliers_and_extract_chisq <- function(matrix_list, distances, conf = 0.95) {
  # Calculate the quantile threshold
  df <- ncol(summary_stats) 
  threshold <- qchisq(conf, df)
  # Identify outlier indices
  outlier_indices <- which(distances > threshold)
  # Extract outlier matrices
  outlier_matrices <- matrix_list[outlier_indices]
  # Return a list with indices and matrices
  list(indices = outlier_indices, matrices = outlier_matrices)
}

outliergene_q90 <- detect_outliers_and_extract_quantiles(matrices,distances,quantile_threshold = 0.90)
outliergene_q95 <- detect_outliers_and_extract_quantiles(matrices,distances,quantile_threshold = 0.95)
outliergene_q99 <- detect_outliers_and_extract_quantiles(matrices,distances,quantile_threshold = 0.99)
outliergene_chi90 <- detect_outliers_and_extract_chisq(matrices,distances,conf = 0.90)
outliergene_chi95 <- detect_outliers_and_extract_chisq(matrices,distances,conf = 0.95)
outliergene_chi99 <- detect_outliers_and_extract_chisq(matrices,distances,conf = 0.99)



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

for(i in seq_along(matrics)){

upper_triangle_list <- apply(which(upper.tri(mat), arr.ind = TRUE), 1, function(idx) {
  # Extract the row name, column name, and value
  row_name <- rownames(mat)[idx[1]]
  col_name <- colnames(mat)[idx[2]]
  value <- mat[idx[1], idx[2]]
  # Return as a list for each cell
  list(row = row_name, column = col_name, value = value)
  })

# Extract the vector of values
upper_values <- sapply(upper_triangle_list, function(x) x$value)

# Create alphabetical combinations of row and column names
upper_combinations <- sapply(upper_triangle_list, function(x) {
  paste(sort(c(x$row, x$column)), collapse = "_")
})

name <- names(matrices)[[1]]
for (i in seq_along(upper_values)) {
  # Fill in the value for the corresponding column
  df[name, upper_combinations[i]] <- upper_values[i]
}

}




create_combined_dataframe <- function(matrices) {
  # Identify all unique upper triangle combinations across all matrices
  
  # Initialize the dataframe
  df <- as.data.frame(matrix(NA, nrow = length(matrices), ncol = length(combinations_vect)))
  colnames(df) <- combinations_vect
  rownames(df) <- names(matrices)
  
  # Populate the dataframe
  for (name in names(matrices)) {
    print(name)
    mat <- matrices[[name]]
    
    # Extract upper triangle information
    upper_triangle_list <- apply(which(upper.tri(mat), arr.ind = TRUE), 1, function(idx) {
      row_name <- rownames(mat)[idx[1]]
      col_name <- colnames(mat)[idx[2]]
      value <- mat[idx[1], idx[2]]
      list(row = row_name, column = col_name, value = value)
    })
    
    # Extract values and combinations
    upper_values <- sapply(upper_triangle_list, function(x) x$value)
    upper_combinations <- sapply(upper_triangle_list, function(x) {
      paste(sort(c(x$row, x$column)), collapse = "_")
    })
    
    # Fill the dataframe for the current matrix
    for (i in seq_along(upper_values)) {
      df[name, upper_combinations[i]] <- upper_values[i]
    }
  }
  return(df)
}

df <- create_combined_dataframe(matrices)


write.csv(df,'Results/gene_speciescombo_uppertri_matrix.csv')

# impute NA's based on gene means 
df <- read.csv('Results/gene_speciescombo_uppertri_matrix.csv',header = T,row.names = 1)
for (i in seq_len(nrow(df))) {
  # Check if the row has NA values
  if (anyNA(df[i, ])) {
    # Replace NA values with the row mean
    row_data <- as.numeric(df[1,])
    row_mean <- mean(row_data, na.rm = TRUE)
    df[i, ][is.na(df[i, ])] <- row_mean
  }
}

head(df,2)

# Compute Mahalanobis distances for matrices
center_matrix <- colMeans(df, na.rm = T)
cov_matrix <- cov(df, use = "pairwise.complete.obs")
distances <- mahalanobis(df, center_matrix, cov_matrix)

detect_outliers_and_extract_quantiles <- function(matrices, distances, quantile_threshold = 0.95) {
  # Calculate the quantile threshold
  threshold <- quantile(distances, quantile_threshold)
  # Identify outlier indices
  outlier_indices <- names(which(distances > threshold))
  # Extract outlier matrices
  outlier_matrices <- matrices[outlier_indices]
  # Return a list with indices and matrices
  list(indices = outlier_indices, matrices = outlier_matrices)
}

detect_outliers_and_extract_chisq <- function(matrix_list, distances, conf = 0.95) {
  # Calculate the quantile threshold
  degrees_freedom <- ncol(df) 
  threshold <- qchisq(conf, degrees_freedom)
  
  # Identify outlier indices
  outlier_indices <- names(which(distances > threshold))
  
  # Extract outlier matrices
  outlier_matrices <- matrices[outlier_indices]
  
  # Return a list with indices and matrices
  list(indices = outlier_indices, matrices = outlier_matrices)
}

outliergene_q90 <- detect_outliers_and_extract_quantiles(matrices,distances,quantile_threshold = 0.90)
outliergene_q95 <- detect_outliers_and_extract_quantiles(matrices,distances,quantile_threshold = 0.95)
outliergene_q99 <- detect_outliers_and_extract_quantiles(matrices,distances,quantile_threshold = 0.99)
outliergene_chi90 <- detect_outliers_and_extract_chisq(matrices,distances,conf = 0.90)
outliergene_chi95 <- detect_outliers_and_extract_chisq(matrices,distances,conf = 0.95)
outliergene_chi99 <- detect_outliers_and_extract_chisq(matrices,distances,conf = 0.99)

write.csv(data.frame(Outlier = outliergene_chi90$indices),'Results/outliergene_chi90.csv',row.names = F)
write.csv(data.frame(Outlier = outliergene_chi95$indices),'Results/outliergene_chi95.csv',row.names = F)
write.csv(data.frame(Outlier = outliergene_chi99$indices),'Results/outliergene_chi99.csv',row.names = F)
write.csv(data.frame(Outlier = outliergene_q90$indices),'Results/outliergene_q90.csv',row.names = F)
write.csv(data.frame(Outlier = outliergene_q95$indices),'Results/outliergene_q95.csv',row.names = F)
write.csv(data.frame(Outlier = outliergene_q99$indices),'Results/outliergene_q99.csv',row.names = F)


