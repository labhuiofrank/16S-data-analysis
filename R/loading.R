library(tibble)
library(dplyr)
library(data.table)

#' Load OTU abundance file from mothur
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
