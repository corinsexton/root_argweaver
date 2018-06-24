library(ape)
library(readr)

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
trees_and_labels <- readTrees(args[1])
trees = trees_and_labels[[1]]
labels = trees_and_labels[[2]]

tree <-read.tree(text = trees[1,3] )
tree$tip.label <- labels[as.numeric(tree$tip.label) + 1]

matrix <- cophenetic(tree)
for (num in 2:nrow(trees)) {

  tree <-read.tree(text = trees[num,3] )
  tree$tip.label <- labels[as.numeric(tree$tip.label) + 1]
 
  new_matrix <- cophenetic(tree)
  matrix <- pmin(matrix, new_matrix)

}

rownames(matrix) <- labels[rownames(matrix)]
colnames(matrix) <- labels[colnames(matrix)]

matrix <- cbind(names = rownames(matrix), matrix)
write.table(matrix,file = args[2],quote = F,row.names = T,col.names = T,sep = '\t')