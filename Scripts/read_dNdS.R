
# Read dN and dS files
n_species <- as.integer(readLines("Data/dN_dS_data/100073at32523_dN.dist", n = 1))

dN <- read.table("Data/dN_dS_data/100073at32523_dN.dist", row.names=1, skip=1, fill=TRUE, 
                 col.names = paste0("V", 1:(n_species+1)))
dS <- read.table("Data/dN_dS_data/100073at32523_dS.dist", row.names=1, skip=1, fill=TRUE, 
                 col.names = paste0("V", 1:(n_species+1)))

# adjust colnames to species
colnames(dN) <- rownames(dN)
colnames(dS) <- rownames(dS)

# look at data :)
dNdS <- (as.matrix(dN)) / (as.matrix(dS))
dNdS <- as.vector(as.matrix(dN)) / as.vector(as.matrix(dS))
summary(dNdS)
hist(dNdS)



# loop through all files --------------------------------------------------

# Define the folder containing the data
data_folder <- "Data/dN_dS_data"

# List all _dN.dist and _dS.dist files
dN_files <- list.files(data_folder, pattern = "_dN\\.dist$", full.names = TRUE)
dS_files <- list.files(data_folder, pattern = "_dS\\.dist$", full.names = TRUE)

# Ensure the files are paired and match
if (!all(gsub("_dN\\.dist$", "", basename(dN_files)) == gsub("_dS\\.dist$", "", basename(dS_files)))) {
  stop("Mismatch between dN and dS files!")
}

# Initialize empty lists to store matrices
dN_matrices <- list()
dS_matrices <- list()


# Loop through all files and read data
for (i in seq_along(dN_files)) {
  # Check if the file is empty
  if (file.info(dN_files[i])$size == 0 || file.info(dS_files[i])$size == 0) {
    warning(paste("Skipping empty file pair:", dN_files[i], "and", dS_files[i]))
    next
  }
  
  matrix_name <- gsub("_dN\\.dist$", "", basename(dN_files[i]))
  
  
  
  # Read the number of species from the first line
  n_species <- as.integer(readLines(dN_files[i], n = 1))
  
  # Load dN and dS matrices
  dN <- read.table(dN_files[i], row.names = 1, skip = 1, fill = TRUE, 
                   col.names = paste0("V", 1:(n_species + 1)))
  dS <- read.table(dS_files[i], row.names = 1, skip = 1, fill = TRUE, 
                   col.names = paste0("V", 1:(n_species + 1)))
  
  # Adjust column names to match species names
  colnames(dN) <- rownames(dN)
  colnames(dS) <- rownames(dS)
  
  # Store matrices in lists
  dN_matrices[[matrix_name]] <- as.matrix(dN)
  dS_matrices[[matrix_name]] <- as.matrix(dS)
}


# Optionally convert lists to arrays if dimensions match
# dN_stack <- array(unlist(dN_matrices), dim = c(dim(dN_matrices[[1]]), length(dN_matrices)))
# dS_stack <- array(unlist(dS_matrices), dim = c(dim(dS_matrices[[1]]), length(dS_matrices)))

# Check the data
summary(dN_matrices)
summary(dS_matrices)

### Divide

# Initialize an empty list to store the dN/dS results
dN_dS_matrices <- list()

# Loop through the names of dN_matrices
for (name in names(dN_matrices)) {
  # Check if there is a corresponding dS matrix
  if (!name %in% names(dS_matrices)) {
    warning(paste("No matching dS matrix found for", name))
    next
  }
  
  # Perform element-wise division
  dN_dS <- dN_matrices[[name]] / dS_matrices[[name]]
  
  # Store the result in the list
  dN_dS_matrices[[name]] <- dN_dS
}

# Check the result
names(dN_dS_matrices)  # Names of the resulting matrices
summary(dN_dS_matrices)  # Summary of the results


# Function to make the lower triangle the upper triangle
transpose_lower_to_upper <- function(matrix) {
  # Create a matrix of the same dimensions with zeros
  new_matrix <- matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix))
  
  # Copy the lower triangle to the upper triangle
  new_matrix[upper.tri(new_matrix)] <- matrix[lower.tri(matrix)]
  
  # Copy the diagonal as well (optional, based on your needs)
  diag(new_matrix) <- diag(matrix)
  colnames(new_matrix) <- colnames(matrix)
  rownames(new_matrix) <- rownames(matrix)
  return(new_matrix)
}

# Apply this transformation to all matrices in the dN_dS_matrices list
dN_dS_matrices_transposed <- lapply(dN_dS_matrices, transpose_lower_to_upper)
names(dN_dS_matrices_transposed)[sapply(dN_dS_matrices_transposed, function(mat) any(mat == 0))]
dN_dS_matrices_transposed['152310at32523']

# Check the result
names(dN_dS_matrices_transposed)  # Names of the matrices
summary(dN_dS_matrices_transposed)
