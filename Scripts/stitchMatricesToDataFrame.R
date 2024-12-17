stitchMatricesToDataFrame <- function(matrices) {
  # Helper function to flatten the upper triangle of a matrix
  flattenUpperTriangle <- function(matrix) {
    if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix)) {
      stop("All elements in the list must be square matrices.")
    }
    upper_tri_indices <- upper.tri(matrix, diag = FALSE)
    return(matrix[upper_tri_indices])
  }
  
  # Iterate over the list of matrices
  flattened_rows <- lapply(seq_along(matrices), function(i) {
    flattened_row <- flattenUpperTriangle(matrices[[i]]) # Flatten the matrix
    return(flattened_row)       # Return as a numeric vector
  })
  
  # Combine flattened rows into a data frame
  stitched_df <- do.call(rbind, flattened_rows)
  
  # Add gene names as row names
  gene_names <- names(matrices)
  rownames(stitched_df) <- gene_names
  
  # Return the final data frame
  return(as.data.frame(stitched_df))
}
