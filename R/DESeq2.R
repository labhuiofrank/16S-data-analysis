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
#' @param factor_name condition variable to test
#' @param ref reference level for variable
#' @param level other level for variable
#' @return result dataframe
#' @examples
#' get_results_for_contrast(diagdds, "Group", "A", "B")
get_results_for_contrast <- function(diagdds, factor_name, ref, level) {
  results(diagdds, contrast=c(factor_name, level, ref))  %>%
    data.frame %>%
    mutate(enriched_in=ifelse(log2FoldChange < 0, ref, level), 
           depleted_in=ifelse(log2FoldChange > 0, ref, level)) %>%
    mutate(log2FoldChange=abs(log2FoldChange)) %>%
    na.omit(padj) %>%
    rownames_to_column("OTU")
}

#' Get result table from DESeq2 object between all levels vs reference
#' If the same OTU is selected for multiple comparisons, the padj and log2FoldChange
#' retained are the ones for the most significant comparison (lowest padj)
#'
#' @param diagdds DESeq2 object
#' @param variable condition variable to test
#' @param comparisons n x 2 dataframe, each row is a comparison to perform
#' @return result dataframe
#' @examples
#' get_all_results(diagdds, "Group", "A", c("B", "C"))
get_all_results <- function(diagdds, variable, comparisons) {
    var_order <- unique(c(comparisons))

    all_results <- list()
    for (i in 1:nrow(comparisons)) {
        all_results[[i]] <- get_results_for_contrast(diagdds, variable, comparisons[i,2], comparisons[i,1])
    }

    do.call(rbind, all_results) %>%
      mutate(enriched_in=factor(enriched_in, var_order), 
             depleted_in=factor(depleted_in, var_order))
}
