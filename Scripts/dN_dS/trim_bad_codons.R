# Script that takes multiple sequence alignment and trims columns in two steps
# 1. Trim columns with missing data 
# 2. Trim columns with high heterozygosity

library(ape)

#' function that function trims alingment columns (or full codons) with less than a given proportion of data.
#' 
#' @param al object of class DNAbin
#' @param prop proportion of data acceptable for column, anything below will lead to the column being removed
#' @param codon boolean whether the alignment is codon-aware, or columns should be treated independently
#' @return trimmed alignment al, with all columns having at least given proportion of data
trimCols <- function(al, prop, codon = T){
	 mat <- as.character(as.matrix(al))
	 ntax <- nrow(mat)
	 propthres <- 1-prop
	 compliantSites <- apply(mat, 2, function(x){
	 		x <- as.character(x)
			compl <- (length(which(x %in% c("N", "n", "?", "-", "O", "o", "X", "x"))) / ntax) < propthres
			return(compl)
			})
			
	 if(codon){
			codIDs <- rep(1:(length(compliantSites)/3), each = 3)
			codsToKeep <- rep(F, length(compliantSites))
			for(i in 1:max(codIDs)){
			      if(all(compliantSites[which(codIDs == i)])) codsToKeep[which(codIDs == i)] <- T
			}
			al <- al[, as.logical(codsToKeep)]
			return(al)
	 } else {
			al <- al[, as.logical(compliantSites)]
			return(al)
	 }
}

#' function that function trims codons with heterozygosity above a given threshold.
#' 
#' @param al object of class DNAbin
#' @param threshet proportion of acceptable heterozygosity for column, anything above will lead to the column being removed
#' @return trimmed alignment al, with all columns having maximum given proportion of heterozygosity
trimCodons <- function(al, threshet = 0.5){
  aaal <- trans(al)
  ncaaal <- ncol(aaal)
  toremove <- which(sapply(1:ncaaal, function(x) max(table(as.character(aaal[,x]))) / nrow(aaal)) < threshet)*3
  if(length(toremove) > 0) al <- al[,-c(toremove, toremove-1, toremove-2)]
  return(al)
}

# read in current alignment file that has been trimmed based on codons
al <- read.dna('BUSCO_alignments_no_stop/1000at32523.fasta', format="fasta")

# get trimmed alignment without gappy columns
al_no_gaps <- trimCols(al, prop = 0.5)

# trimmed alignment alignment even further and remove heterozygous columns
al_no_heterozygous <- trimCodons(al_no_gaps)

# write trimmed alignment into new directory
write.dna(al_no_heterozygous, 'BUSCO_alignments_trimmed/1000at32523.fa', 
          format="fasta", colsep = "")

