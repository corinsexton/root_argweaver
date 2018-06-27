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

# trees_and_labels <- readTrees("root_hla/HLA-DRB1_tests/HLA-DRB1_maxtime_2000e3_recomb_100xslower/smc_files_HLA-DRB1_maxtime_2000e3_recomb_100xslower/HLA-DRB1_maxtime_2000e3_recomb_100xslower.1000.smc.gz")
# sites <- getSites("/Users/coripenrod/root_hla/HLA-DRB1_tests/HLA-DRB1_maxtime_2000e3_recomb_100xslower/HLA-DRB1_maxtime_2000e3_recomb_100xslower.sites")

trees = trees_and_labels[[1]]
labels = trees_and_labels[[2]]



# # NONHUMAN
# labels[grepl('Gogo|Patr|Popy|Poab|Papa', labels)]
# nonhuman_ind <- grep('Gogo|Patr|Popy|Poab|Papa', labels, value = F)
# 
# # HUMAN
# labels[!grepl('Gogo|Patr|Popy|Poab|Papa', labels)]
# human_grepl <- !grepl('Gogo|Patr|Popy|Poab|Papa', labels)
# #grep('Gogo|Patr|Popy|Poab|Papa', labels, value = F)

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




# minima
# 
# index <- which.min(dist_mrca[1,])
# min_mrca <- as.numeric(dist_mrca[2, index])
# start <- as.numeric(dist_mrca[3, index])
# end <- as.numeric(dist_mrca[4, index])
# min_dist <- min(unlist(dist_mrca[1,]))

# distance_vector <- sapply(treeList,function(l) {
#   l$node_inf$distances[hap1,hap2]
# })



# for (hap1 in 1:length(ordering)) {
#   for (hap2 in 1:length(ordering)) {
#       if (hap1 != hap2 && !(paste(hap2,hap1, sep = '-') %in% seen)) {
#         
#         name_haplo <- paste(hap1,hap2, sep = '-')
#         #nodeList[[count]] <- list(nodepath = nodepath(tree, hap1, hap2),
#         nodeList[[count]] <- list(
#                               mrca_node = mrca_matrix[ordering[hap1],ordering[hap2]],
#                               distance = dist_matrix[hap1,hap2],
#                               mrca_distance = dist_matrix[hap1, mrca_matrix[hap1,hap2]])
#         names(nodeList)[count] <- name_haplo
#         
#         seen[[count]] <- name_haplo
#         count <- count + 1
#       }  
#   }
# }








  #### ACROSS ALL TREES:
  # across trees::::
  #   
  # - minimum distance between them
  # - what tree it was found in (sequence)
  # - 5 minimum distances
  # - path through the tree was the minimum distance

    # - minimum distance greater than (pick mins) (across all trees)
  
  # 
  # 
  # tree <-read.tree(text = trees[num,3])
  # tree$tip.label <- labels[as.numeric(tree$tip.label) + 1]
  # 
  # 
  # 
  # 
  # 
  # new_matrix <- cophenetic(tree)
  # matrix <- pmin(matrix, new_matrix)
  # 
  # 
  # 
# evaluateTree <- function(tree, human_grepl) {
#   #### BASED ON EACH TREE: 
#   # min distance
#   dist_matrix <- cophenetic(tree)
#   dist_matrix <- dist.topo(tree)
#   z <- dist.nodes(tree)
#   #[1] 384 809 808 806 805 804 802 801 800 799 794 793 740 741 747 749 750 751 752 765 341
#   z[341,765] + z[765,752] + z[753,751] + z[750,749] + z[749,747]  + z[747,741] + z[741,740]
#   
#   # - spec - spec: closer to a non species than to each other (T/F)
#   # - spec - non spec: closer to a species than this non species (T/F)
#   species_matrix <- matrix(nrow = nrow(dist_matrix),
#                            ncol = ncol(dist_matrix),
#                            dimnames = list(rownames(dist_matrix), colnames(dist_matrix)))
#   for (i in 1:ncol(dist_matrix)) {
#     
#     current_node <- as.numeric(colnames(dist_matrix)[i])
#     node_is_human <- human_grepl[current_node]
#     
#     # distances for nonhumans and humans
#     nonhuman_dist_matrix <- dist_matrix[!human_grepl,]
#     human_dist_matrix <- dist_matrix[human_grepl,]
#     
#     min_node_nonhuman <- rownames(dist_matrix[-i,])[which.min(dist_matrix[-i,][,i])]
#     min_dist_nonhuman <- min(dist_matrix[-i,][,i])
#     
#     min_node_human <- rownames(dist_matrix[-i,])[which.min(dist_matrix[-i,][,i])]
#     min_dist_human <- min(dist_matrix[-i,][,i])
#     
#     same_species <- min_node %in% nonhuman_ind == current_node %in% nonhuman_ind
#     if (node_is_human) {
#       ### Closest node is same species
#       
#       
#     } else {
#       ### Closest node is different species
#     }
#     # equal mins?
#     #rownames(dist_matrix[-i,])[which(dist_matrix[-i,][,i] == min(dist_matrix[-i,][,i]))]
#     #min(dist_matrix[-i,][,i])
#   }
#   
#   # find last common ancestor node
#   x <- mrca(tree)
#   
# }
  
# }
# 
# rownames(matrix) <- labels[rownames(matrix)]
# colnames(matrix) <- labels[colnames(matrix)]
# 
# matrix <- cbind(names = rownames(matrix), matrix)
# write.table(matrix,file = args[2],quote = F,row.names = T,col.names = T,sep = '\t')
# 
# 
# tree <-read.tree(text = trees[1,3] )
# tree2 <-read.tree(text = trees[2,3] )
# 
# tree$tip.label <- labels[as.numeric(tree$tip.label) + 1]
# tree2$tip.label <- labels[as.numeric(tree2$tip.label) + 1]
# 
# 
# cophyloplot(tree,tree2)
