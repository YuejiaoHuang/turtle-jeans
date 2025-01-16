# Create a list of test matrices with unequal sizes and specified structure
test_matrices_symmetric <- list(
  # Matrix 1: 3x3, symmetric, filled with non-zero values
  matrix(c(0, 1, 2,
           1, 0, 3,
           2, 3, 0), nrow = 3, byrow = TRUE),
  
  # Matrix 2: 4x4, lower triangle empty (0s and NAs) but upper triangle has non-zero values
  matrix(c(0, 1, 2, 3,
           NA, 0, 4, 5,
           NA, NA, 0, 6,
           NA, NA, NA, 0), nrow = 4, byrow = TRUE),
  
  # Matrix 3: 5x5, upper triangle empty (0s and NAs) but lower triangle has non-zero values
  matrix(c(0, NA, NA, NA, NA,
           1, 0, NA, NA, NA,
           2, 1, 0, NA, NA,
           3, 2, 1, 0, NA,
           4, 3, 2, 1, 0), nrow = 5, byrow = TRUE),
  
  # Matrix 4: 6x6, symmetric, filled with non-zero values
  matrix(c(0, 1, 2, 3, 4, 5,
           1, 0, 6, 7, 8, 9,
           2, 6, 0, 10, 11, 12,
           3, 7, 10, 0, 13, 14,
           4, 8, 11, 13, 0, 15,
           5, 9, 12, 14, 15, 0), nrow = 6, byrow = TRUE),
  
  # Matrix 5: 7x7, lower triangle empty (0s and NAs) but upper triangle has non-zero values
  matrix(c(0, 1, 2, 3, 4, 5, 6,
           NA, 0, 7, 8, 9, 10, 11,
           NA, NA, 0, 12, 13, 14, 15,
           NA, NA, NA, 0, 16, 17, 18,
           NA, NA, NA, NA, 0, 19, 20,
           NA, NA, NA, NA, NA, 0, 21,
           NA, NA, NA, NA, NA, NA, 0), nrow = 7, byrow = TRUE)
)

# You can now use this list for testing
test_matrices_symmetric

mat <- matrices[[1]]

matrixCheck <- function(matrix_list) {
  # Check if input is a non-empty list
  if (!is.list(matrix_list) || length(matrix_list) == 0) {
    stop("Input must be a non-empty list of matrices.")
  }
  
  # Check if all elements in the list are matrices
  if (!all(sapply(matrix_list, is.matrix))) {
    stop("All elements in the list must be matrices.")
  }
  
  # Function to test if upper is empty or not, and transpose 
  transposeLower <- function(mat) {
    if(all(mat[upper.tri(mat)] == 0 | is.na(mat[upper.tri(mat)]))){
      mat_t <- t(mat)
      return(mat_t)
    }else{
      return(mat)
    }
  }
  
  matrix_list_t <- lapply(matrix_list,transposeLower)
  return(matrix_list_t)
}

test <- matrixCheck(test_matrices_symmetric)

