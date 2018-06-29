#!/usr/bin/Rscript

suppressWarnings(library(ape, quietly = T))
suppressWarnings(library(readr, quietly = T))

#' Get Pairwise Information
#'
#' @param treeList Generated previously. List containing all trees for a smc file.
#' @param hap1 a string for a haplotype
#' @param hap2 a string for a haplotype
#' @param labels labels provided by smcfile
#'
#' @return data frame, describes for all regions 1. start of region, 2. end of region, 3. number of variant sites in a region, 4. distance between the two haplotypes, and 5. the mrca node id.
#' @export
#'
#' @examples
#' 
#' DRB1_04_03_to_04_05_df <- getPairwiseInformation(treeList,
#'                             "DRB1_04.03.01.01", "DRB1_04.05.01.01", labels)
#' 
getPairwiseInformation <- function(treeList, hap1, hap2, labels) {
  
  if (!(hap1 %in% labels) | !(hap2 %in% labels)) {
    cat("Please enter two valid haplotype names. Options are: \n\n")
    cat(labels, sep="\t")
    cat("\n")
    return(FALSE)
  }
  
  # from "PatrNDRB1_02.09" -> "121"
  hap_num1 <- which(labels == hap1)
  hap_num2 <- which(labels == hap2)
  
  dist_mrca <- sapply(treeList,function(l) {
    
    # from 121 -> 24 (per tree ordering) (0 indexed in labels)
    h1 <- which(l$node_label_ordering == as.character(hap_num1 - 1))
    h2 <- which(l$node_label_ordering == as.character(hap_num2 - 1))
    
    list(start_of_region = l$start_region,
         end_of_region = l$end_region,
         num_sites_in_region = l$num_sites_in_region,
         distance_between = l$node_inf$distances[h1,h2],
         mrca = l$node_inf$mrca_matrix[h1,h2]
    )
  })
  
  x <- data.frame(t(dist_mrca))
  x$start_of_region <- as.numeric(x$start_of_region)
  x$end_of_region <- as.numeric(x$end_of_region)
  x$num_sites_in_region <- as.numeric(x$num_sites_in_region)
  x$distance_between <- as.numeric(x$distance_between)
  x$mrca <- as.numeric(x$mrca)
  
  x
}

###########

suppressWarnings(library(argparse, quietly=TRUE))

# Create a parser
parser <- ArgumentParser(description = "Pairwise Allele Information Parser")

# Add command line arguments
parser$add_argument("haplotype1", nargs = 1, help="Haplotype 1")
parser$add_argument("haplotype2", nargs = 1, help="Haplotype 2")
parser$add_argument("output", help="output file")
parser$add_argument("--rdata", action="store_true",
                    help="path to treeList.RData [defaults to ./treeList.RData]",
                    default='./treeList.RData')

argv <- parser$parse_args()

load(argv$rdata)
df <- getPairwiseInformation(treeList, argv$haplotype1, argv$haplotype2, labels)
if (class(df) == "data.frame")  write_tsv(df, path = argv$output)
