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
  
  if (class(min_split) == 'data.frame') {
    demes_available <- colnames(min_split)[2:length(colnames(min_split))]
    allele1_deme <- 'human'
    for (deme in demes_available) {
      if (grepl(tolower(deme), allele, ignore.case = T)) {
        allele1_deme <- tolower(deme)
        break
      }
    }
  } else {
    if (grepl('gogo', allele, ignore.case = T)) {
      allele1_deme <- 'gogo'
    } else if (grepl('popy', allele, ignore.case = T)) {
      allele1_deme <- 'popy'
    } else if (grepl('patr', allele, ignore.case = T)) {
      allele1_deme <- 'patr'
    } else if (grepl('poab', allele, ignore.case = T)) {
      allele1_deme <- 'poab'
    } else if (grepl('papa', allele, ignore.case = T)) {
      allele1_deme <- 'papa'
    } else {
      allele1_deme <- 'human'
    }
  }
  
  total_df <- tibble::tibble(allele_1_id = character(),
                             allele_2_id = character(),
                             deme_allele_1 = character(),
                             deme_allele_2 = character(),
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
    
    allele_num2 <- which(labels == other_alleles[i])
    
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
    
    x <- data.frame(t(dist_mrca))
    x$start_of_region <- as.numeric(x$start_of_region)
    x$end_of_region <- as.numeric(x$end_of_region)
    x$num_sites_in_region <- as.numeric(x$num_sites_in_region)
    x$distance_between <- as.numeric(x$distance_between)
    x$mrca <- as.numeric(x$mrca)
    
    allele2 <- other_alleles[i]
    
    if (class(min_split) == 'data.frame') {
      allele2_deme <- 'human'
      
      for (deme in demes_available) {
        if (grepl(tolower(deme), allele2, ignore.case = T)) {
          allele2_deme <- tolower(deme)
          break
        }
      }
    } else {
        if (grepl('gogo', allele2, ignore.case = T)) {
          allele2_deme <- 'gogo'
        } else if (grepl('popy', allele2, ignore.case = T)) {
          allele2_deme <- 'popy'
        } else if (grepl('patr', allele2, ignore.case = T)) {
          allele2_deme <- 'patr'
        } else if (grepl('poab', allele2, ignore.case = T)) {
          allele2_deme <- 'poab'
        } else if (grepl('papa', allele2, ignore.case = T)) {
          allele2_deme <- 'papa'
        } else {
          allele2_deme <- 'human'
        }
    }
    
    if (class(min_split) == 'data.frame') {
      colnames(min_split) <- lapply(colnames(min_split),tolower)
      min_split[,1] <- sapply(min_split[,1],tolower)
      
      min_spl <- min_split[which(min_split[,1] == allele1_deme), which(colnames(min_split) == allele2_deme)]
    } else {
      min_spl <- min_split
    }
    
    if (class(max_split) == 'data-frame') {
      colnames(max_split) <- lapply(colnames(max_split),tolower)
      max_split[,1] <- sapply(max_split[,1],tolower)
      
      max_spl <- max_split[which(max_split[,1] == allele1_deme), which(colnames(max_split) == allele2_deme)]
    } else {
      max_spl <- max_split
    }
    
    minimum_dist <- subset(x, distance_between >= min_spl)
    maximum_dist <- subset(x, distance_between <= max_spl)

    if (nrow(minimum_dist) == 0) {
      cat("WARNING: Minimum split time specified is too ancient (too high), no MRCA found for: ", allele, " and ", other_alleles[i], '\n')
      min_data_dist[[j]] <- 'NA'
      min_data_start[[j]] <- 'NA'
      min_data_end[[j]] <- 'NA'
      min_data_sites[[j]] <- 'NA'
      min_data_mrca[[j]] <- 'NA'
    } else {
      min_dist <- min(minimum_dist$distance_between)
      min_indices <- which(minimum_dist$distance_between == min_dist)
      
      min_data_start <- list()
      min_data_end <- list()
      min_data_sites <- list()
      min_data_mrca <- list()
      for (j in 1:length(min_indices)) {
        min_index <- min_indices[j]
        
        min_data_dist <- minimum_dist$distance_between[min_index]
        min_data_start[[j]] <- minimum_dist$start_of_region[min_index]
        min_data_end[[j]] <- minimum_dist$end_of_region[min_index]
        min_data_sites[[j]] <- minimum_dist$num_sites_in_region[min_index]
        min_data_mrca[[j]] <- minimum_dist$mrca[min_index]
      }
    }
    
    
    if (nrow(maximum_dist) == 0) {
      cat("WARNING: Maximum split time specified is too recent (too low), no MRCA found for: ", allele, " and ", other_alleles[i], '\n')
      max_data_dist[[j]] <- 'NA'
      max_data_start[[j]] <- 'NA'
      max_data_end[[j]] <- 'NA'
      max_data_sites[[j]] <- 'NA'
      max_data_mrca[[j]] <- 'NA'
    } else {
      max_dist<- max(maximum_dist$distance_between)
      max_indices <- which(maximum_dist$distance_between == max_dist)
      
      max_data_start <- list()
      max_data_end <- list()
      max_data_sites <- list()
      max_data_mrca <- list()
      for (j in 1:length(max_indices)) {
        max_index <- max_indices[j]
        
        max_data_dist <- maximum_dist$distance_between[max_index]
        max_data_start[[j]] <- maximum_dist$start_of_region[max_index]
        max_data_end[[j]] <- maximum_dist$end_of_region[max_index]
        max_data_sites[[j]] <- maximum_dist$num_sites_in_region[max_index]
        max_data_mrca[[j]] <- maximum_dist$mrca[max_index]
      }
    }
    
    total_df[i,] <- c(allele, other_alleles[i], allele1_deme, allele2_deme,
                    min_data_dist,
                    paste(min_data_start,collapse = ','),
                    paste(min_data_end,collapse = ','),
                    paste(min_data_sites,collapse = ','),
                    paste(min_data_mrca,collapse = ','),
                    max_data_dist,
                    paste(max_data_start,collapse = ','),
                    paste(max_data_end,collapse = ','),
                    paste(max_data_sites,collapse = ','),
                    paste(max_data_mrca,collapse = ',')
                    )
  }
  as.data.frame(total_df)
}

###########

suppressWarnings(library(argparse, quietly=TRUE))

# Create a parser
parser <- ArgumentParser(description = "Single Allele Information Parser")

# Add command line arguments
parser$add_argument("allele", nargs = 1, help="Allele Name")
parser$add_argument("output", help="output file")

parser$add_argument("--min_split", help="Matrix specifying minimum allowable
                                         time for a MRCA to occur (tab separated)")
parser$add_argument("--max_split", help="Matrix specifying maximum allowable
                                         time for a MRCA to occur (tab separated)")
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
