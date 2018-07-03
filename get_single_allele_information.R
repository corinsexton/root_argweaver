#!/usr/bin/Rscript

suppressWarnings(library(ape, quietly = T))
suppressWarnings(library(readr, quietly = T))

#' Get Single Allele Information
#'
#' @param treeList Generated previously. List containing all trees for a smc file.
#' @param allele name of allele (as found in smc file)
#' @param labels labels provided by smcfile
#' @param min_split Number specifying minimum allowable time for a MRCA to occur
#' @param max_split Number specifying maximum allowable time for a MRCA to occur
#'
#' @return data frame, describes for all other alleles 1. first allele id, 2. second allele id, 3. first allele deme, 4. second allele deme, 5. shortest distance found to second allele, 6. start of region for tree with the shortest distance, 7. end of region for tree with the shortest distance, 8. number of variant sites in the region for the shortest distance tree, 9. mrca node in shortest distance tree, 10. maximum distance found to that allele, 11. start of region for tree with the maximum distance, 12. end of region for tree with the maximum distance, 13. number of variant sites in the region for the maximum distance tree, 14. mrca node in maximum distance tree
#' @export
#'
#' @examples
#' df <- getSingleAlleleInformation(treeList, "DRB1_12.01.01.05", labels)

getSingleAlleleInformation <- function(treeList, allele, labels, min_split = 0, max_split = Inf) {
  
  if (!(allele %in% labels)) {
    cat("Please enter a valid allele name. Options are: \n\n")
    cat(labels, sep="\t")
    cat("\n")
    return(FALSE)
  }
  
  # from "PatrNDRB1_02.09" -> "121"
  allele_num1 <- which(labels == allele)
  other_alleles <- labels[which(labels != allele)]
  
  allele1_deme <- getDemeIdentifier(min_split, allele)
  
  total_df <- tibble::tibble(allele_1_id = character(), allele_2_id = character(),
                             deme_allele_1 = character(), deme_allele_2 = character(),
                             shortest_dist = numeric(),
                             shortest_dist_region_start = numeric(),
                             shortest_dist_region_end = numeric(),
                             shortest_dist_num_sites_in_region = numeric(),
                             shortest_mrca = numeric(),
                             longest_dist = numeric(),
                             longest_dist_region_start = numeric(),
                             longest_dist_region_end = numeric(),
                             longest_dist_num_sites_in_region = numeric(),
                             longest_mrca = numeric()
                             )

  for (i in 1:length(other_alleles)) {
    
    allele2 <- other_alleles[i]
    allele_num2 <- which(labels == allele2)
    
    # pull data for all trees
    dist_mrca <- sapply(treeList,function(l) {
      
      # from 121 -> 24 (per tree ordering) (0 indexed in labels)
      h1 <- which(l$node_label_ordering == as.character(allele_num1 - 1))
      h2 <- which(l$node_label_ordering == as.character(allele_num2 - 1))
      
      list(start_of_region = l$start_region,
           end_of_region = l$end_region,
           num_sites_in_region = l$num_sites_in_region,
           distance_between = l$node_inf$distances[h1,h2],
           mrca = l$node_inf$mrca_matrix[h1,h2]
           
      )
    })
    
    # format all trees data with proper types
    dist_mrca_df <- data.frame(t(dist_mrca))
    dist_mrca_df$start_of_region <- as.numeric(dist_mrca_df$start_of_region)
    dist_mrca_df$end_of_region <- as.numeric(dist_mrca_df$end_of_region)
    dist_mrca_df$num_sites_in_region <- as.numeric(dist_mrca_df$num_sites_in_region)
    dist_mrca_df$distance_between <- as.numeric(dist_mrca_df$distance_between)
    dist_mrca_df$mrca <- as.numeric(dist_mrca_df$mrca)
<<<<<<< HEAD
    
    allele2_deme <- getDemeIdentifier(min_split, allele2)
    
    # get min and max distances according to deme
    min_spl <- getSplitTime(min_split, allele1_deme, allele2_deme)
    max_spl <- getSplitTime(max_split, allele1_deme, allele2_deme)

    # subset by min and max distances
    minimum_dist <- subset(dist_mrca_df, distance_between >= min_spl)
    maximum_dist <- subset(dist_mrca_df, distance_between <= max_spl)

    # determine mins and maxs 
    if (nrow(minimum_dist) == 0) {
      cat("WARNING: Minimum split time specified is too ancient (too high), no MRCA found for: ", allele, " and ", other_alleles[i], '\n')
    }
    min_list <- getDistancePerAllele(minimum_dist, min)
    
    if (nrow(maximum_dist) == 0) {
      cat("WARNING: Maximum split time specified is too recent (too low), no MRCA found for: ", allele, " and ", other_alleles[i], '\n')
    }
    max_list <- getDistancePerAllele(maximum_dist, max)
    
    total_df[i,] <- c(allele, allele2, allele1_deme, allele2_deme,
                    min_list[[1]],
                    paste(min_list[[2]], collapse = ','),
                    paste(min_list[[3]], collapse = ','),
                    paste(min_list[[4]], collapse = ','),
                    paste(min_list[[5]], collapse = ','),
                    max_list[[1]],
                    paste(max_list[[2]], collapse = ','),
                    paste(max_list[[3]], collapse = ','),
                    paste(max_list[[4]], collapse = ','),
                    paste(max_list[[5]], collapse = ',')
                    )
=======
    
    allele2_deme <- getDemeIdentifier(min_split, allele2)
    
    # get min and max distances according to deme
    min_spl <- getSplitTime(min_split, allele1_deme, allele2_deme)
    max_spl <- getSplitTime(max_split, allele1_deme, allele2_deme)

    # subset by min and max distances
    minimum_dist <- subset(dist_mrca_df, distance_between >= min_spl)
    maximum_dist <- subset(dist_mrca_df, distance_between <= max_spl)

    # determine mins and maxs 
    if (nrow(minimum_dist) == 0) {
      cat("WARNING: Minimum split time specified is too ancient (too high),
           no MRCA found for: ", allele, " and ", other_alleles[i], '\n')
    }
    min_list <- getDistancePerAllele(minimum_dist, min)
    
    if (nrow(maximum_dist) == 0) {
      cat("WARNING: Maximum split time specified is too recent (too low),
           no MRCA found for: ", allele, " and ", other_alleles[i], '\n')
    }
    max_list <- getDistancePerAllele(maximum_dist, max)
    
    total_df[i,] <- c(allele, allele2, allele1_deme, allele2_deme,
                    min_list[[1]],
                    paste(min_list[[2]], collapse = ','),
                    paste(min_list[[3]], collapse = ','),
                    paste(min_list[[4]], collapse = ','),
                    paste(min_list[[5]], collapse = ','),
                    max_list[[1]],
                    paste(max_list[[2]], collapse = ','),
                    paste(max_list[[3]], collapse = ','),
                    paste(max_list[[4]], collapse = ','),
                    paste(max_list[[5]], collapse = ',')
                    )
  }
  as.data.frame(total_df)
}


#' Get Deme Indetifier (internal)
#'
#' @param min_split either 1) data frame containing identical row and column names for deme identifiers OR 2) an integer to specify min split limits for all species
#' @param allele allele name
#'
#' @return the deme identifier for that allele
#' @export
#'
#' @examples
getDemeIdentifier <- function(min_split, allele) {
  if (class(min_split) == 'data.frame') {
    demes_available <- colnames(min_split)[2:length(colnames(min_split))]
    allele_deme <- 'human'
    
    for (deme in demes_available) {
      if (grepl(tolower(deme), allele, ignore.case = T)) {
        allele_deme <- tolower(deme)
        break
      }
    }
  } else {
    if (grepl('gogo', allele, ignore.case = T)) {
      allele_deme <- 'gogo'
    } else if (grepl('popy', allele, ignore.case = T)) {
      allele_deme <- 'popy'
    } else if (grepl('patr', allele, ignore.case = T)) {
      allele_deme <- 'patr'
    } else if (grepl('poab', allele, ignore.case = T)) {
      allele_deme <- 'poab'
    } else if (grepl('papa', allele, ignore.case = T)) {
      allele_deme <- 'papa'
    } else {
      allele_deme <- 'human'
    }
  }
  return(allele_deme)
}


#' Get Distance Per Allele (internal)
#'
#' @param subsetted_matrix matrix for all the trees for one allele pair
#' @param limit_func either max or min
#'
#' @return a list with either all max or all min statistics
#' @export
#'
#' @examples
getDistancePerAllele <- function(subsetted_matrix, limit_func) {
  if (nrow(subsetted_matrix) == 0) {
    data_dist[[j]] <- 'NA'
    data_start[[j]] <- 'NA'
    data_end[[j]] <- 'NA'
    data_sites[[j]] <- 'NA'
    data_mrca[[j]] <- 'NA'
  } else {
    dist<- limit_func(subsetted_matrix$distance_between)
    indices <- which(subsetted_matrix$distance_between == dist)
    
    data_start <- list()
    data_end <- list()
    data_sites <- list()
    data_mrca <- list()
    for (j in 1:length(indices)) {
      index <- indices[j]
      
      data_dist <- subsetted_matrix$distance_between[index]
      data_start[[j]] <- subsetted_matrix$start_of_region[index]
      data_end[[j]] <- subsetted_matrix$end_of_region[index]
      data_sites[[j]] <- subsetted_matrix$num_sites_in_region[index]
      data_mrca[[j]] <- subsetted_matrix$mrca[index]
    }
>>>>>>> 3c78d63f31b1bdc0b5ae5ed9dca572a4da72f36f
  }
  return(list(data_dist,data_start,data_end,data_sites,data_mrca))
}


<<<<<<< HEAD
#' Get Deme Indetifier (internal)
#'
#' @param min_split either 1) data frame containing identical row and column names for deme identifiers OR 2) an integer to specify min split limits for all species
#' @param allele allele name
#'
#' @return the deme identifier for that allele
#' @export
#'
#' @examples
getDemeIdentifier <- function(min_split, allele) {
  if (class(min_split) == 'data.frame') {
    demes_available <- colnames(min_split)[2:length(colnames(min_split))]
    allele_deme <- 'human'
    
    for (deme in demes_available) {
      if (grepl(tolower(deme), allele, ignore.case = T)) {
        allele_deme <- tolower(deme)
        break
      }
    }
  } else {
    if (grepl('gogo', allele, ignore.case = T)) {
      allele_deme <- 'gogo'
    } else if (grepl('popy', allele, ignore.case = T)) {
      allele_deme <- 'popy'
    } else if (grepl('patr', allele, ignore.case = T)) {
      allele_deme <- 'patr'
    } else if (grepl('poab', allele, ignore.case = T)) {
      allele_deme <- 'poab'
    } else if (grepl('papa', allele, ignore.case = T)) {
      allele_deme <- 'papa'
    } else {
      allele_deme <- 'human'
    }
  }
  return(allele_deme)
}


#' Get Distance Per Allele (internal)
#'
#' @param subsetted_matrix matrix for all the trees for one allele pair
#' @param limit_func either max or min
#'
#' @return a list with either all max or all min statistics
#' @export
#'
#' @examples
getDistancePerAllele <- function(subsetted_matrix, limit_func) {
  if (nrow(subsetted_matrix) == 0) {
    data_dist[[j]] <- 'NA'
    data_start[[j]] <- 'NA'
    data_end[[j]] <- 'NA'
    data_sites[[j]] <- 'NA'
    data_mrca[[j]] <- 'NA'
  } else {
    dist<- limit_func(subsetted_matrix$distance_between)
    indices <- which(subsetted_matrix$distance_between == dist)
    
    data_start <- list()
    data_end <- list()
    data_sites <- list()
    data_mrca <- list()
    for (j in 1:length(indices)) {
      index <- indices[j]
      
      data_dist <- subsetted_matrix$distance_between[index]
      data_start[[j]] <- subsetted_matrix$start_of_region[index]
      data_end[[j]] <- subsetted_matrix$end_of_region[index]
      data_sites[[j]] <- subsetted_matrix$num_sites_in_region[index]
      data_mrca[[j]] <- subsetted_matrix$mrca[index]
    }
  }
  return(list(data_dist,data_start,data_end,data_sites,data_mrca))
}


=======
>>>>>>> 3c78d63f31b1bdc0b5ae5ed9dca572a4da72f36f
#' Get Split Time
#'
#' @param split_data matrix containing split times or a numeric
#' @param deme1 allele 1 deme
#' @param deme2 allele 2 deme
#'
#' @return
#' @export
#'
#' @examples
getSplitTime <- function(split_data, deme1, deme2) {
  if (class(split_data) == 'data-frame') {
    colnames(split_data) <- lapply(colnames(split_data),tolower)
    split_data[,1] <- sapply(split_data[,1],tolower)
    
    max_spl <- split_data[which(split_data[,1] == deme1), which(colnames(split_data) == deme2)]
  } else {
    max_spl <- split_data
  }
}
###########

suppressWarnings(library(argparse, quietly=TRUE))

# Create a parser
parser <- ArgumentParser(description = "Single Allele Information Parser")

# Add command line arguments
parser$add_argument("allele", nargs = 1, help="Allele Name")
parser$add_argument("output", help="output file")

parser$add_argument("--min_split", help="Matrix specifying minimum allowable time for a MRCA to occur (tab separated)")
parser$add_argument("--max_split", help="Matrix specifying maximum allowable time for a MRCA to occur (tab separated)")
parser$add_argument("--rdata", action="store_true",
                    help="path to treeList.RData [defaults to ./treeList.RData]",
                    default='./treeList.RData')

argv <- parser$parse_args()

load(argv$rdata)

if (!is.null(argv$min_split)) {
  min_split <- read.table(argv$min_split, header = T, stringsAsFactors = F)
} else {
  min_split = 0
}

if (!is.null(argv$max_split)) {
  max_split <- read.table(argv$max_split, header = T, stringsAsFactors = F)
} else {
  max_split = Inf
}

df <- getSingleAlleleInformation(treeList, argv$allele, labels, min_split, max_split)
if (class(df) == "data.frame")  write_tsv(df, path = argv$output)

# CAVEATS: A max split matrix can only be specified with a min split matrix. Min split can be specified
# without a max split matrix. Defaults to looking for 'popy', 'gogo', 'patr', 'poab','papa', and 'human' and
# is not case sensitive. (IMPORTANT) If only some of these are specified in the min split matrix, the taxa
# not specified will be considered human.
