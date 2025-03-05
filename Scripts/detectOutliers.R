detect_outliers_and_extract_quantiles <- function(matrices, df, quantile_threshold = 0.95) {
  # Compute Mahalanobis distances for matrices
  center_matrix <- colMeans(df, na.rm = T)
  cov_matrix <- cov(df, use = "pairwise.complete.obs")
  
  #return mahalanobis distances
  distances <- mahalanobis(df, center_matrix, cov_matrix)
  
  # Calculate the quantile threshold
  threshold <- quantile(distances, quantile_threshold)
  # Identify outlier indices
  outlier_indices <- names(which(distances > threshold))
  # Extract outlier matrices
  outlier_matrices <- matrices[outlier_indices]
  # Return a list with indices and matrices
  list(indices = outlier_indices, matrices = outlier_matrices)
}

detect_outliers_and_extract_chisq <- function(matrix_list, df, conf = 0.95) {
  # Compute Mahalanobis distances for matrices
  center_matrix <- colMeans(df, na.rm = T)
  cov_matrix <- cov(df, use = "pairwise.complete.obs")
  
  #return mahalanobis distances
  distances <- mahalanobis(df, center_matrix, cov_matrix)
  
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

