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
    if (tolower(method) == "pcoa") {
        components <- ord.mod$vectors
    } else {
        components <- scores(ord.mod, display="sites")
    }

    ## Add metadata to components
    metadata <- sample_data(ps)[rownames(components),] %>% data.frame
    components <- components %>% data.frame %>%
        bind_cols(metadata)

    results <- list(data=components, model=ord.mod)
    
    if (tolower(method) %in% c("cca", "rda")) {
        ## Supervised --> we also extract the arrows for the biplot
        results$arrows <- scores(ord.mod, display="bp")
    }
    
    return(results)
}

#' Wrapper around vegan::adonis
#'
#' @param ps Phyloseq object
#' @param meta_vars Variable(s) to test
#' @param ... Extra parameters to pass to the adonis function
#' @return dataframe with r2 and pvalue for PERMANOVA test
#' @examples
#' run_adonis(ps, "Group+Season", permutations=999)
run_adonis <- function(ps, meta_vars, ...) {
    dists <- phyloseq::distance(ps, method="bray")
    meta <- data.frame(sample_data(ps))
    formula <- as.formula(sprintf("dists ~ %s", meta_vars))
    adonis.tab <- adonis(formula, data=meta, ...)$aov.tab
    adonis_res <- data.frame(r2=adonis.tab$R2[1], pval=adonis.tab$`Pr(>F)`[1])

    return(adonis_res)
}
