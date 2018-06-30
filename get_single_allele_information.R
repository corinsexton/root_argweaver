#!/usr/bin/Rscript

suppressWarnings(library(ape, quietly = T))
suppressWarnings(library(readr, quietly = T))


# A haplotype's-eye-view script, which for an input-specified haplotype i of an input-specified locus, outputs a file comprising a row for each haplotype j that includes
# . j's haplotype id
# . j's taxon
# . the shortest* distance separating haplotypes i and j, as found among all trees in the run
# . the segment(s) of the locus where that shortest* distance was inferred,including (as subfields in a delimited list) 
#     . the start_end of each such segment in question
#     . the node id of the oldest node on that shortest-found path between the two haplotypes in that given segment;
#     . the sites count within that segment
# . the longest* distance separating haplotypes i and j among all trees in the run
# . the segment(s) of the locus where that longest* distance was inferred, including (as subfields in a delimited list)
#     . the start_end of each such segment in question
#     . the node id of the oldest node on that longest-found path between the two haplotypes in that given segment
#     . the count of variable (or otherwise 'informative', if defined better by ARGweaver) sites within that segment

# *: where 'shortest' can be restricted to >= an optional --min_split parameter, and 'longest' can be restricted to <= an optional --max_split parameter




#' Get Single Allele Information
#'
#' @param treeList Generated previously. List containing all trees for a smc file.
#' @param hap name of allele (as found in smc file)
#' @param labels labels provided by smcfile
#' @param min_split Number specifying minimum allowable time for a MRCA to occur
#' @param max_split Number specifying maximum allowable time for a MRCA to occur
#'
#' @return data frame, describes for all other alleles 1. allele id, 2. taxon of allele, 3. shortest distance found to that allele, 4. start of region for tree with the shortest distance, 5. end of region for tree with the shortest distance, 6. number of variant sites in the region for the shortest distance tree, 7. mrca node in shortest distance tree, 8. maximum distance found to that allele, 9. start of region for tree with the maximum distance, 10. end of region for tree with the maximum distance, 11. number of variant sites in the region for the maximum distance tree, 12. mrca node in maximum distance tree
#' @export
#'
#' @examples
#' df <- getSingleAlleleInformation(treeList, "DRB1_12.01.01.05", labels)
getSingleAlleleInformation <- function(treeList, hap, labels, min_split = 0, max_split = Inf, all_species = FALSE) {
  
  if (!(hap %in% labels)) {
    cat("Please enter a valid haplotype name. Options are: \n\n")
    cat(labels, sep="\t")
    cat("\n")
    return(FALSE)
  }
  
  # from "PatrNDRB1_02.09" -> "121"
  hap_num1 <- which(labels == hap)
  other_haps <- labels[which(labels != hap)]
  
  if (grepl('gogo', hap, ignore.case = T)) {
    hap1_taxon <- 'Gorilla_gorilla'
  } else if (grepl('popy', hap, ignore.case = T)) {
    hap1_taxon <- 'Pongo_pygmaeus'
  } else if (grepl('patr', hap, ignore.case = T)) {
    hap1_taxon <- 'Pan_troglodytes'
  } else if (grepl('poab', hap, ignore.case = T)) {
    hap1_taxon <- 'Pongo_abelii'
  } else if (grepl('papa', hap, ignore.case = T)) {
    hap1_taxon <- 'Pan_paniscus'
  } else {
    hap1_taxon <- 'Homo_sapiens'
  }
  
  total_df <- tibble::tibble(allele_1_id = character(),
                             allele_2_id = character(),
                             taxon_allele_2 = character(),
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

  for (i in 1:length(other_haps)) {
    
    hap_num2 <- which(labels == other_haps[i])
    
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
    
    
    if (grepl('gogo', other_haps[i], ignore.case = T)) {
      taxon <- 'Gorilla_gorilla'
    } else if (grepl('popy', other_haps[i], ignore.case = T)) {
      taxon <- 'Pongo_pygmaeus'
    } else if (grepl('patr', other_haps[i], ignore.case = T)) {
      taxon <- 'Pan_troglodytes'
    } else if (grepl('poab', other_haps[i], ignore.case = T)) {
      taxon <- 'Pongo_abelii'
    } else if (grepl('papa', other_haps[i], ignore.case = T)) {
      taxon <- 'Pan_paniscus'
    } else {
      taxon <- 'Homo_sapiens'
    }
    
    if (class(min_split) == 'data-frame') {
      min_split <- min_split[which(min_split[,1] == hap1_taxon), which(colnames(min_split) == hap1_taxon)]
    }
    
    if (class(max_split) == 'data-frame') {
      max_split <- max_split[which(max_split[,1] == hap1_taxon), which(colnames(max_split) == hap1_taxon)]
    }
    
    min_index <- which.min(subset(x, distance_between >= min_split)$distance_between) # TIE BREAKERS? (only 20 unique)
    max_index <- which.max(subset(x, distance_between <= max_split)$distance_between) # TIE BREAKERS?
    
    minimum_dist <- subset(x, distance_between >= min_split)
    maximum_dist <- subset(x, distance_between <= max_split)
    
    if (nrow(minimum_dist) == 0) {
      cat("WARNING: Minimum split time specified is too ancient (too high), no MRCA found for: ", hap, " and ", other_haps[i], '\n')
      min_data <- c(NA, NA, NA, NA, NA)
    } else {
      min_index <- which.min(minimum_dist$distance_between) # TIE BREAKERS? (only 20 unique)
      
      min_data <- c(minimum_dist$distance_between[min_index],
                    minimum_dist$start_of_region[min_index],
                    minimum_dist$end_of_region[min_index],
                    minimum_dist$num_sites_in_region[min_index],
                    minimum_dist$mrca[min_index])
    }
    
    if (nrow(maximum_dist) == 0) {
      cat("WARNING: Maximum split time specified is too recent (too low), no MRCA found for: ", hap, " and ", other_haps[i], '\n')
      max_data <- c(NA, NA, NA, NA, NA)
    } else {
      max_index <- which.max(maximum_dist$distance_between) # TIE BREAKERS?
      
      max_data <- c(maximum_dist$distance_between[max_index],
                    maximum_dist$start_of_region[max_index],
                    maximum_dist$end_of_region[max_index],
                    maximum_dist$num_sites_in_region[max_index],
                    maximum_dist$mrca[max_index])
    }

    total_df[i,] <- c(hap, other_haps[i],
                     taxon, min_data, max_data)
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
parser$add_argument("-a", "--all_species", default=FALSE,
                    help="enforce split time parameters for all species (including the same species as the allele in question) [default: False]")
parser$add_argument("--rdata", action="store_true",
                    help="path to treeList.RData [defaults to ./treeList.RData]",
                    default='./treeList.RData')

argv <- parser$parse_args()

load(argv$rdata)

if (!is.null(argv$min_split)) {
  min_split <- read_table(min_split)
} else {
  min_split = 0
}

if (!is.null(argv$max_split)) {
  max_split <- read_table(max_split)
} else {
  max_split = Inf
}

df <- getSingleAlleleInformation(treeList, argv$allele,
                                 labels, min_split, max_split,
                                 as.logical(argv$all_species))
# df <- getSingleAlleleInformation(treeList, hap, labels, min_split, max_split)

if (class(df) == "data.frame")  write_tsv(df, path = argv$output)



