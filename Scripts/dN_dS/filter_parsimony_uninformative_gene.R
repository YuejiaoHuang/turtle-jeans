# Script that takes multiple sequence alignment and copies it to new location
# if it contains at least 10 parsimony informative sites.

library(ape)

# read in current alignment file that has been trimmed based on codons
al <- read.dna('BUSCO_alignments_trimmed/99895at32523.fa', format="fasta")

#' function that identifies parsimony informative sites
#' 
#' @param x object of class DNAbin
#' @param what string that defines output of function, can be "absolute", "fraction", or "index
#' @return number, fraction or indices (depending on value of what param) of parsimony informative sites in x
pis <- function (x, what = "fraction"){
  if (!inherits(x, "DNAbin"))
    stop("'x' is not of class 'DNAbin'")
  what <- match.arg(what, c("absolute", "fraction", "index"))
  pars.inf <- function(x) {
    x <- table(x)
    x <- x[x > 1]
    n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w",
           "y")
    if (length(x[!names(x) %in% n]) > 1)
      x <- TRUE
    else x <- FALSE
  }
  x <- as.character(x)
  out <- apply(x, 2, pars.inf)
  if (what %in% c("absolute", "fraction")) {
    # get absolute number of pis
    out <- length(out[out])
    if (what == "fraction") {
      # get fraction of pis
      out <- round(out/ncol(x) * 100, digits = 2)
    }
  }
  else {
    # get indices of pis
    out <- which(out)
  }
  out
}

# get indices of parsimony informative sites of alignment
parsimony_informative_sites <- pis(al, what = "index")

# print number of parsimony informative sites 
# (used later to get file with n for each gene and look at distribution)
print(length(parsimony_informative_sites))

# filter out all genes that have less 10 parsimony informative sites
if (length(parsimony_informative_sites) > 9) {
  
  # filtering is done by only coping those alignments to the new location that fulfill criterion
  write.dna(al, 'BUSCO_alignments_trimmed_final_pis/99895at32523.fa', 
            format="fasta", colsep = "")
  
}
