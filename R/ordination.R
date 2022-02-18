library(phyloseq)
library(vegan)

#' Wrapper around phyloseq's ordinate function to compute NMDS and include metadata
#' This function needs to be updated with more specific behaviour for other types
#' of ordination methods (PCoA, CCA, ...)
#'
#' @param ps Phyloseq object. Needs to contain a phylogenetic tree if the distance matrix is unifrac-based
#' @param method Ordinaton method. See `ordinate()` function to see all the available methods
#' @param ... Extra arguments to provide to the underlying ordination function
#' @return Dataframe with first 2 components and metadata information
#' @examples
#' run_nmds(ps, method="NMDS", distance="bray", trymax=300)
#' run_nmds(ps, method="PCoA")
run_nmds <- function(ps, method="NMDS", ...) {
    
    ord.mod <- ordinate(ps, method=method, ...)

    if (method == "NMDS") {
        print(sprintf("stress=%s, converged=%s", ord.mod$stress, ord.mod$converged))
        components <- scores(ord.mod)
    }

    df <- sample_data(ps)[rownames(components),]
    df[, sprintf("%s%s", method, 1:2)] <- components[, 1:2]
    
    return(df)
}
