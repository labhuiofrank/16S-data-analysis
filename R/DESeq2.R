#' Wrapper around DESeq functions
#'
#' @param ps Phyloseq object
#' @param design Experimental design. The last variable needs to be the experimental condition to test for
#' @param auto_gm Whether to automatically use modified geometric mean for size factors when too many zeros are present
#' @param gm Whether or not to use geometric mean for size factor estimation. Might be necessary if there are a lot of zero (or add pseudocount
#' @param pseudocount Pseudocount to add to abundance table if there are too many zeros
#' @param ... Extra arguments to pass the DESeq function
#' @return DESeq2 object
#' @export
#' @examples
#' run_deseq2(ps, ~ Season + Group, gm=TRUE)
#' run_deseq2(ps, ~ Season + Group, pseudocount=1)
run_deseq2 <- function(ps, design, auto_gm=TRUE, gm=FALSE, pseudocount=0, ...) {
    if (auto_gm && (gm || pseudocount > 0)) {
        warning("Setting auto_gm to FALSE since one of {gm, pseudocount} has been set")
        auto_gm = FALSE
    }
    if (auto_gm) {
        ## Normalize if all OTUs have at least one zero
        taxa_with_zeros <- filter_taxa(ps, function(x) sum(x == 0) > 0, TRUE)
        if (ntaxa(taxa_with_zeros) == ntaxa(ps)) {
            message(
                "All OTUs have at least one zero. ",
                "Using geometric mean on positive values for size factor estimation"
            )
            gm = TRUE
        }
    }
    ## Add pseudocount if needed
    phyloseq::otu_table(ps) <- phyloseq::otu_table(ps) + pseudocount

    ## Convert to DESeqDataSet object
    diagdds <- phyloseq::phyloseq_to_deseq2(ps, design)

    ## Use geometric mean for size factors if needed
    if (gm) {
        geom_mean <- function(x) exp(mean(log(x[x>0]), na.rm=TRUE))
        geoMeans <- apply(DESeq2::counts(diagdds), 1, geom_mean)
        diagdds <- DESeq2::estimateSizeFactors(diagdds, geoMeans=geoMeans)
    }

    ## Run DESeq2
    DESeq2::DESeq(diagdds, ...)
}

#' Get result table from DESeq2 object between given variable level vs reference
#'
#' @param diagdds DESeq2 object
#' @param factor_name condition variable to test
#' @param ref reference level for variable
#' @param level other level for variable
#' @return result dataframe
#' @examples
#' get_results_for_single_contrast(diagdds, "Group", "A", "B")
get_results_for_single_contrast <- function(diagdds, factor_name, ref, level) {
  DESeq2::results(diagdds, contrast=c(factor_name, level, ref))  %>% data.frame %>%
    dplyr::mutate(enriched_in=ifelse(log2FoldChange < 0, ref, level), 
                  depleted_in=ifelse(log2FoldChange > 0, ref, level)) %>%
    dplyr::mutate(log2FoldChange=abs(log2FoldChange)) %>%
    dplyr::na.omit(padj) %>%
    tibble::rownames_to_column("OTU")
}

#' Get result table from DESeq2 object between all levels vs reference
#' If the same OTU is selected for multiple comparisons, the padj and log2FoldChange
#' retained are the ones for the most significant comparison (lowest padj)
#'
#' @param diagdds DESeq2 object
#' @param variable condition variable to test
#' @param comparisons n x 2 dataframe, each row is a comparison to perform
#' @return result dataframe
#' @export
#' @examples
#' get_results_for_contrasts(diagdds, "Group", c("B", "C"))
get_results_for_contrasts <- function(diagdds, variable, comparisons) {
    var_order <- unique(c(comparisons))

    all_results <- list()
    for (i in 1:nrow(comparisons)) {
        all_results[[i]] <- get_results_for_contrast(diagdds, variable, comparisons[i,2], comparisons[i,1])
    }

    do.call(rbind, all_results) %>%
      dplyr::mutate(enriched_in=factor(enriched_in, var_order), 
                    depleted_in=factor(depleted_in, var_order))
}
