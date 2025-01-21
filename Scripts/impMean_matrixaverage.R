library('gdata')

#matrices <- dN_checked

impMean_matrices <- function(matrices) {

  AddColAndRow <- function(matrix,listsp){
    difference <- setdiff(listsp,rownames(matrix))
        # only att empty columns if any columns are missing
    if (length(difference) > 0){
      for(i in difference){
        
        col <- matrix(NA, ncol = 1, nrow = nrow(matrix))
        colnames(col) <- i
        matrix <- cbind(matrix,col)
        
        row <- matrix(NA, nrow = 1, ncol = ncol(matrix))
        rownames(row) <- i
        matrix <- rbind(matrix,row)
      }
    }
  
    #newmat <- matrix[listsp,listsp]
    sort_order <- order(colnames(matrix))
    sorted_matrix <- matrix[sort_order,sort_order]
  
    return(sorted_matrix)
  }
  
  # Identify all unique species across matrices
  sp.per.mat <- lapply(matrices, rownames)
  listsp <- sort(unique(unlist(sp.per.mat)))

  # Extend all matrices to the same dimensions
  matrices.extended <- lapply(matrices, AddColAndRow,listsp)

  # get means of s*s through all genes 
  mean.list <- matrix(nrow = length(listsp), ncol = length(listsp), 
                      dimnames = list(listsp, listsp))
  
  for (i in listsp) {
    for(j in listsp){
      values <- sapply(matrices.extended, function(mat) mat[i, j])
      mean <- mean(values,na.rm = T)
      mean.list[i,j] <- mean
    }
  }
  
  # Replace missing values in a matrix with its mean value
  ReplaceMissingValueWithMean <- function(matrix,mean.list) {
    
    if(any(is.na(matrix[upper.tri(matrix)]))){
      
      mean_g <- mean(upper.tri(matrix), na.rm = TRUE) # Compute the mean ignoring NA
      na_indices <- upper.tri(matrix) & is.na(matrix)
      na_indices <- which(upper.tri(matrix) & is.na(matrix), arr.ind = T)
  
      # Identify indices where NA exists in the upper triangle
      na_indices <- which(upper.tri(matrix) & is.na(matrix), arr.ind = TRUE)
      
      # Loop over each NA location
      for (i in seq_len(nrow(na_indices))) {
        row_idx <- na_indices[i, 1] # Row index
        col_idx <- na_indices[i, 2] # Column index
        
        row_name <- rownames(matrix)[row_idx] # Row name
        col_name <- colnames(matrix)[col_idx] # Column name
        
        # Retrieve the mean for the specific species combination
        mean_ss <- mean.list[row_name, col_name]
        
        # Calculate the combined mean (global mean and specific mean)
        mean_ss_g <- mean(c(mean_g, mean_ss), na.rm = TRUE)
        
        # Assign the value to the matrix
        matrix[row_idx, col_idx] <- mean_ss_g
      }
    }
      
    # ensure lower and upper triangle are the same 
    lowerTriangle(matrix,byrow=T) <- upperTriangle(matrix)
    
    # fill in the diagonal 
    # Calculate the row-wise mean excluding the diagonal and set it to the diagonal
    for (i in 1:nrow(matrix)) {
      row_values <- matrix[i, -i]  # Exclude the diagonal element
      diag(matrix)[i] <- mean(row_values)  # Assign the mean of the remaining row elements
    }
    
    return(matrix)
  }
  
  matrices.final <- lapply(matrices.extended,ReplaceMissingValueWithMean,mean.list = mean.list)
  return(matrices.final)
  
}

  
