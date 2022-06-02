library(tibble)
library(tidyr)
library(dplyr)
library(data.table)

#' Load OTU abundance file from mothur (shared format)
#'
#' @param file Path to shared file
#' @return Abundance dataframe (sample x OTU)
#' @examples
#' load_shared("abundance.100.shared")
load_shared <- function(file) {
    abundance <- fread(file, drop=c("label", "numOtus"), header=T, blank.lines.skip=T) %>%
        tibble %>%
        column_to_rownames("Group")
    return(abundance)
}

#' Load OTU abundance file from mothur (count table format)
#'
#' @param file Path to count_table file
#' @return Abundance dataframe (OTU x sample)
#' @examples
#' load_counttable("OTUs.100.rep.count_table")
load_counttable <- function(file) {
    abundance <- fread(file, drop=c("total"), header=T, blank.lines.skip=T) %>%
        tibble %>%
        column_to_rownames("Representative_Sequence")
    return(abundance)
}

#' Load consensus taxonomy from mothur
#' 
#' @param file Path to cons.taxonomy file
#' @return Taxonomy dataframe (OTU x taxonomic ranks)
#' @examples
#' load_constaxonomy("annotations.100.cons.taxonomy")
load_constaxonomy <- function(file) {
    tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    taxonomy <- read.table(file, header=T, row.names=1, sep="\t") %>%
        separate(Taxonomy, tax_ranks, sep=";") %>%
        select(-c(Size, Species))
    return(taxonomy)
}
