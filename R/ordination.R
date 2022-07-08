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

    ord.mod <- phyloseq::ordinate(ps, method=method, ...)

    ## Extract components
    if (tolower(method) == "pcoa") {
        components <- ord.mod$vectors
    } else {
        components <- vegan::scores(ord.mod, display="sites")
    }

    ## Add metadata to components
    metadata <- phyloseq::sample_data(ps)[rownames(components),] %>% data.frame
    components <- components %>% data.frame %>%
        dplyr::bind_cols(metadata)

    results <- list(data=components, model=ord.mod)
    
    if (tolower(method) %in% c("cca", "rda")) {
        ## Supervised --> we also extract the arrows for the biplot
        results$arrows <- vegan::scores(ord.mod, display="bp", scale="species") %>%
            data.frame
    }
    
    return(results)
}

#' Plot of constrained ordination plot with environment variable effect shown as arrows
#' @param data Ordination dataframe with metadata (long format for components)
#' @param arrows Environment factor (arrows) to overlay
#' @param x Variable name for x-axis
#' @param y Variable name for y-axis
#' @param fill Variable name for points fill color
#' @param size Point size
#' @param arrow_rescale Rescale all arrows with by a factor
#' @param ellipse Whether to draw ellipses around fill groups
#' @param ... Additional params to pass to geom_point()
#' @return ggplot2 object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 stat_ellipse
#' @export
#' @examples
#' plot_ordination_with_arrows(df, arrows, "NMDS1", "NMDS2", fill="condition")
plot_ordination_with_arrows <- function(data, arrows, x="CCA1", y="CCA2", fill="black", size=3,
                                        arrow_rescale=0.5, ellipse=TRUE, ...) {
    scale_factor <- arrow_rescale * max(abs(data[, c(x, y)])) / max(abs(arrows))
    arrows <- arrows * scale_factor
    label_nudge <- max(arrows) * 0.2

    g <- ggplot(data) +
        geom_point(aes_string(x=x, y=y, fill=fill), size=size, shape=21, ...) +
        geom_segment(data=arrows, aes_string(x=0, y=0, xend=x, yend=y), 
                     arrow=arrow(length=unit(0.2, "cm"))) + 
        geom_text(data=arrows, aes_string(x=x, y=y), label=rownames(arrows),
                  cex=6, nudge_x=label_nudge, nudge_y=label_nudge, fontface='bold')

    if(ellipse) {
        g <- g + stat_ellipse(geom="polygon", aes_string(x=x, y=y, fill=fill), alpha=0.2, show.legend=FALSE, level=0.75)
    }

    return(g)
}

#' Wrapper around vegan::adonis
#'
#' @param ps Phyloseq object
#' @param meta_vars Variable(s) to test
#' @param ... Extra parameters to pass to the adonis function
#' @return dataframe with r2 and pvalue for PERMANOVA test
#' @export
#' @examples
#' run_adonis(ps, "Group+Season", permutations=999)
run_adonis <- function(ps, meta_vars, ...) {
    dists <- phyloseq::distance(ps, method="bray")
    meta <- data.frame(phyloseq::sample_data(ps))
    formula <- as.formula(sprintf("dists ~ %s", meta_vars))
    adonis.tab <- vegan::adonis(formula, data=meta, ...)$aov.tab
    adonis_res <- data.frame(r2=adonis.tab$R2[1], pval=adonis.tab$`Pr(>F)`[1])

    return(adonis_res)
}
