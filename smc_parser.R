#!/usr/bin/Rscript

suppressWarnings(library(ape))
suppressWarnings(library(readr))

args = commandArgs(trailingOnly=TRUE)
readTrees <- function(smcfile) {
  x <- read.delim(smcfile,
                  stringsAsFactors = F)
  treeLines <- x[which(x$NAMES == "TREE"),2:4]
  labels <- colnames(x)[-1]
  
  colnames(treeLines) <- c("seq1", "seq2", "tree")
  treeLines$seq1 <- as.numeric(treeLines$seq1)
  treeLines$seq2 <- as.numeric(treeLines$seq2)
  return(list(treeLines, labels))
}
getSites <- function(sitesfile) {
  
  sites <- suppressMessages(suppressWarnings(read_tsv(sitesfile, col_names = TRUE)))
  sites <- sites[2:nrow(sites),1]
  sites <- as.numeric(sites$NAMES)
  
  sites
}

trees_and_labels <- readTrees(args[1])
sites <- getSites(args[2])

trees = trees_and_labels[[1]]
labels = trees_and_labels[[2]]

nodeList <- list()
treeList <- list()

for (num in 1:nrow(trees)) {
  
  tree <-read.tree(text = trees[num,3] )

  # PER TREE CALCULATIONS
  ordering <- tree$tip
  
  dist_matrix <- dist.nodes(tree)
  mrca_matrix <- mrca(tree)

  nodeList <- list(distances = dist_matrix, mrca_matrix = mrca_matrix, order = tree$tip)
  
  start_region = as.numeric(trees[num,1])
  end_region = as.numeric(trees[num,2])
  
  num_sites_in_region <- sum(sites >= start_region & sites <= end_region)
  
  treeList[[num]] <- list(start_region = start_region,
                          end_region = end_region,
                          num_sites_in_region = num_sites_in_region,
                          node_info = nodeList,
                          node_label_ordering = ordering)
  
}

save(treeList, labels, file = "treeList.RData")
