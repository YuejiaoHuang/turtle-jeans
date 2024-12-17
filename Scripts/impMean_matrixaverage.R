impMean_matrices <- function(matrices) {
  # Add missing rows and columns to match all species across matrices
  AddColAndRow <- function(matrix, allrowsandcol) {
    newmat <- matrix(nrow = length(allrowsandcol), ncol = length(allrowsandcol), 
                     dimnames = list(allrowsandcol, allrowsandcol))
    newmat[rownames(matrix), colnames(matrix)] <- matrix
    newmat
  }
  
  # Reorder rows and columns to a uniform order
  ReorderColAndRow <- function(matrix, allrowsandcol) {
    newmat <- matrix[allrowsandcol, allrowsandcol]
    newmat
  }
  
  # Replace missing values in a matrix with its mean value
  ReplaceMissingValueWithMean <- function(matrix) {
    mean_value <- mean(matrix, na.rm = TRUE) # Compute the mean ignoring NA
    matrix[is.na(matrix)] <- mean_value     # Replace NA with the mean
    matrix
  }
  
  # Step 1: Identify all unique species across matrices
  sp.per.mat <- lapply(matrices, rownames)
  listsp <- unique(unlist(sp.per.mat))
  
  # Step 2: Extend all matrices to the same dimensions
  matrices.extended <- lapply(matrices, AddColAndRow, allrowsandcol = listsp)
  
  # Step 3: Replace missing values in each matrix with the mean value of that matrix
  matrices.imputed <- lapply(matrices.extended, ReplaceMissingValueWithMean)
  
  # Step 4: Ensure the order of rows and columns is consistent across matrices
  matrices.final <- lapply(matrices.imputed, ReorderColAndRow, allrowsandcol = listsp)
  
  names(matrices.final) <- names(matrices) # Preserve matrix names
  return(matrices.final)
}
