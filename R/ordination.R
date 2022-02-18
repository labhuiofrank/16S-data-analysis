library(phyloseq)
library(vegan)
library(dplyr)

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

    ## Extract components
    if (method == "PCoA") {
        components <- ord.mod$vectors
    } else {
        components <- scores(ord.mod, display="sites")
    }

    ## Add metadata to components
    metadata <- sample_data(ps)[rownames(components),] %>% data.frame
    components <- components %>% data.frame %>%
        bind_cols(metadata)

    results <- list(data=components, model=ord.mod)
    
    if (method %in% c("CCA", "RDA")) {
        ## Supervised --> we also extract the arrows for the biplot
        results$arrows <- scores(ord.mod, display="bp")
    }
    
    return(results)
}
