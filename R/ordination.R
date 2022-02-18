library(phyloseq)
library(vegan)

#' Wrapper around phyloseq's ordinate function
#' Adds the metadata to the component scores
#'
#' @param ps Phyloseq object. Needs to contain a phylogenetic tree if the distance matrix is unifrac-based
#' @param method Ordinaton method. See `ordinate()` function to see all the available methods
#' @param ... Extra arguments to provide to the underlying ordination function
#' @return List with ordination model, dataframe with first 2 components and metadata information, and biplot data if supervised
#' @examples
#' run_ordination(ps, method="NMDS", distance="bray", trymax=300)
#' run_ordination(ps, method="PCoA")
#' run_ordination(ps, "CCA", formula=~ DO + pH + SPC)
run_ordination <- function(ps, method="NMDS", ...) {

    ord.mod <- ordinate(ps, method=method, ...)

    if (method == "PCoA") {
        components <- ord.mod$vectors
    } else {
        components <- scores(ord.mod, display="sites")
    }
    components <- as.data.frame(components)
    n_components <- dim(components)[2]

    df <- sample_data(ps)[rownames(components),]
    df[, sprintf("%s%s", method, 1:2)] <- components[, 1:n_components]

    results <- list(data=df, model=ord.mod)
    
    if (method %in% c("CCA", "RDA")) {
        ## Supervised --> extract the vectors for biplot
        results$arrows <- scores(ord.mod, display="bp")
    }
    
    return(results)
}
