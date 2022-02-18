library(tibble)
library(dplyr)
library(data.table)


load_shared <- function(file) {
    abundance <- fread(file, drop=c("label", "numOtus"), header=T, blank.lines.skip=T) %>%
        tibble %>%
        column_to_rownames("Group")
    return(abundance)
}

load_constaxonomy <- function(file) {
    tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    taxonomy <- read.table(file, header=T, row.names=1, sep="\t") %>%
        separate(Taxonomy, tax_ranks, sep=";") %>%
        select(-c(Size, Species))
    return(taxonomy)
}
