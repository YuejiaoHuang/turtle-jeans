
matrices <- readRDS('Scripts/test_zscores/matrices.rds')
df <- readRDS('Scripts/test_zscores/df.rds')
quantile_threshold = 0.95

detect_outliers_and_extract_quantiles <- function(matrices, df, quantile_threshold = 0.95) {
  # Compute Mahalanobis distances for matrices
  center_matrix <- colMeans(df, na.rm = T)
  cov_matrix <- cov(df, use = "pairwise.complete.obs")
  #return mahalanobis distances
  distances <- mahalanobis(df, center_matrix, cov_matrix)
  
  ### get chisq values 
  df_dim <- ncol(df)
  pchisq <- pchisq(distances, df = df_dim, lower.tail = FALSE)
  
  # Calculate the quantile threshold
  threshold <- quantile(distances, quantile_threshold)
  
  # Identify outlier indices
  outlier_indices <- names(which(distances > threshold))
  # Extract outlier matrices
  outlier_matrices <- matrices[outlier_indices]
  # Return a list with indices and matrices
  list(indices = outlier_indices, matrices = outlier_matrices,
       distances_all = distances, pchisq_all = pchisq)
}

t <- detect_outliers_and_extract_quantiles(matrices,df,quantile_threshold)
str(t)
### Z-scores 
# scale 
hist(distances)
z_scores <- scale(distances)
p_values_normal <- pnorm(z_scores, lower.tail = FALSE)
# get value 
vec <- p_values_normal[,1]
# named vector
names(vec) <- rownames(p_values_normal)
# index 
hist(vec[outlier_indices])


