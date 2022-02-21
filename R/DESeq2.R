library(DESeq2)
library(phyloseq)
library(tibble)
library(dplyr)


#' Geometric mean function for vectors with zeros
geom_mean <- function(x) exp(mean(log(x[x>0]), na.rm=TRUE))

#' Wrapper around DESeq functions
#'
#' @param ps Phyloseq object
#' @param design Experimental design. The last variable needs to be the experimental condition to test for
#' @param gm Whether or not to use geometric mean for size factor estimation. Might be necessary if there are a lot of zero (or add pseudocount
#' @param pseudocount Pseudocount to add to abundance table if there are too many zeros
#' @return DESeq2 object
#' @examples
#' run_deseq2(ps, ~ Season + Group, gm=TRUE)
#' run_deseq2(ps, ~ Season + Group, pseudocount=1)
run_deseq2 <- function(ps, design, gm=FALSE, pseudocount=0) {
    ## Add pseudocount if needed
    otu_table(ps) <- otu_table(ps) + pseudocount

    ## Convert to DESeqDataSet object
    diagdds <- phyloseq_to_deseq2(ps, design)

    ## Use geometric mean for size factors if needed
    if (gm) {
        geoMeans <- apply(counts(diagdds), 1, geom_mean)
        diagdds <- estimateSizeFactors(diagdds, geoMeans=geoMeans)
    }

    ## Run DESeq2
    DESeq(diagdds, fitType="local", test="Wald", parallel=TRUE)
}

#' Get result table from DESeq2 object between given variable level vs reference
#'
#' @param diagdds DESeq2 object
#' @param variable condition variable to test
#' @param ref reference level for variable
#' @param level other level for variable
#' @return result dataframe
#' @examples
#' get_results_for_contrast(diagdds, "Group", "A", "B")
get_results_for_contrast <- function(diagdds, variable, ref, level) {
  results(diagdds, contrast=c(variable, level, ref))  %>%
    data.frame %>%
    select(baseMean, log2FoldChange, padj) %>% 
    mutate(which=level) %>%
    rownames_to_column("OTU")
}

#' Get result table from DESeq2 object between all levels vs reference
#' If the same OTU is selected for multiple comparisons, the padj and log2FoldChange
#' retained are the ones for the most significant comparison (lowest padj)
#'
#' @param diagdds DESeq2 object
#' @param variable condition variable to test
#' @param ref reference level for variable
#' @param levels all levels for comparison with ref
#' @return result dataframe
#' @examples
#' get_all_results(diagdds, "Group", "A", c("B", "C"))
get_all_results <- function(diagdds, variable, ref, levels) {
    all_results <- list()

    for (level in levels) {
        all_results[[level]] <- get_results_for_contrast(diagdds, variable, ref, level)
    }

    all_results <- do.call(rbind, all_results) %>%
        group_by(OTU) %>%
        slice_min(padj, with_ties=FALSE) %>%
        na.omit(padj) %>%
        mutate(signif=ifelse(log2FoldChange < 0,
                             sprintf("overexpressed in %s", ref),
                             sprintf("overexpressed in %s", which))) %>%
        mutate(signif=factor(signif, sprintf("overexpressed in %s", c(ref, levels)))) %>%
        arrange(signif)

    return(all_results)
}