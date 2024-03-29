#' Wrapper around speedyseq::merge_samples2 for combination of factors
#'
#' @param ps Phyloseq object
#' @param groups variables to define sample groups
#' @return Phyloseq object with sample merged
#' @export
#' @examples
#' merge_samples3(ps, c("Site", "Season"))
merge_samples3 <- function(ps, groups) {
    all_cols <- sapply(groups, function(v) phyloseq::get_variable(ps, v))
    phyloseq::sample_data(ps)$combined <- apply(all_cols, 1, paste0, collapse="_")
    ps <- speedyseq::merge_samples2(ps, "combined")
    sample_data(ps)$combined <- NULL

    return(ps)
}

#' Palette generator for taxonomic data
#' Uses database in github database repository
#'
#' @param taxa_names Vector of taxa names (factor)
#' @return named vector, values are colors in HEX and names are taxa names
#' @export
#' @examples
#' get_taxonomic_palette(c("Alphaproteobacteria", "Zetaproteobacteria", "Other"))
get_taxonomic_palette <- function(taxa_names) {
    stopifnot(is.factor(taxa_names))
    taxa_names <- levels(taxa_names)
    db_path <- system.file("extdata", "taxa_palette.csv", package="labhuiofrank.16S")
    ## tibble::deframe converts a 2 column dataframe to named vector
    custom_pal <- read.csv(db_path) %>% tibble::deframe()

    is_in_palette <- taxa_names %in% names(custom_pal)
    is_in_taxa <- names(custom_pal) %in% taxa_names
    # is_in_palette <- taxa_names %in% names(custom_pal)
    ## In case when we have unknown taxa names, we pick among the remaining colors
    n_colors_remaining <- sum(!is_in_taxa)
    n_taxa_remaining <- sum(!is_in_palette)
    replace <- ifelse(n_taxa_remaining > n_colors_remaining, TRUE, FALSE)

    ## Assign remaining taxa to colors
    remaining_colors <- sample(custom_pal[!is_in_palette], n_taxa_remaining, replace=replace)
    names(remaining_colors) <- taxa_names[!is_in_palette]
    
    ## Special case of "Other" group
    filler_idx <- grep("[Oo]ther.*", names(remaining_colors))
    if (any(filler_idx)) {
        remaining_colors[filler_idx] <- "gray"
    }

    colors <- c(custom_pal[is_in_taxa], remaining_colors)
    return(colors[taxa_names])
}

#' Taxonomic stacked barplot of phyloseq object.
#'
#' @param ps Phyloseq object to plot
#' @param x Categorical variable in metadata, displayed in x-axis
#' @param y Variable for the height of bars. Either "Abundance" or
#' "relabund" (the latter normalize the bars to the same height)
#' @param taxrank Taxonomic rank to gather taxa (stacked for each level of x)
#' @param rows Factor for row facets
#' @param cols Factor for col facets
#' @param min_relabund Below this threshold, replace taxa name to "< x%" to
#' make the plot for readable and reduce the number of colors
#' @param return_df Whether to return the formatted dataframe or plot
#' @return ggplot object (geom_bar) or a dataframe if return_df is TRUE
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' 
#' @export
#' @examples
#' taxa_barplot(ps, x="sample_type", taxrank="Phylum", min_relabund=0.01)
#' taxa_barplot(ps, x="sample_type", y="relabund", taxrank="Phylum", rows="Season")
taxa_barplot <- function(ps, x=NULL, y="relabund", taxrank="Class",
                         rows=NULL, cols=NULL, min_relabund=0.01,
                         return_df=FALSE, nrows_legend=20) {
    # Group taxa by {tax_rank}
    ps <- speedyseq::tax_glom(ps, taxrank)
    # Group samples for each combination of factors in [x, rows, cols]
    ps <- labhuiofrank::merge_samples3(ps, c(x, rows, cols))
    # Transform phyloseq object to data frame and compute relative abundances
    data <- speedyseq::psmelt(ps) %>% dplyr::group_by(combined) %>%
        dplyr::mutate(relabund=Abundance/sum(Abundance))
    # Replace low abundance species with filler "< xx %"
    filler <- sprintf('Other(<%.f%%)', 100*min_relabund)
    data[data$relabund < min_relabund, taxrank] <- filler
    
    # Get the taxa ordering (by total relabund)
    sorted_taxa <- data %>% dplyr::group_by_at(taxrank) %>%
        dplyr::summarize(score=sum(relabund)) %>%
        dplyr::arrange(score) %>% dplyr::pull(1)

    # Order the taxa and put the filler at the end
    data <- data %>% dplyr::mutate_at(taxrank, ~ factor(., sorted_taxa))
    if (filler %in% sorted_taxa) {
        data <- data %>% dplyr::mutate_at(taxrank, ~ relevel(., filler))
    }

    # Return
    if(return_df) {
        return(data)
    }

    # Choose color palette from database folder
    palette <- get_taxonomic_palette(data %>% dplyr::pull(taxrank))
    
    p <- ggplot(data=data, aes_string(x=x, y=y, fill=taxrank)) +
        geom_bar(stat="identity", position="stack", color="black", size=0.5, width=0.7) +
        scale_fill_manual(values=palette) +
        guides(fill=guide_legend(reverse=T, nrow=nrows_legend))

    # Facet arguments, need to handle null values
    facets <- sprintf(
        "%s ~ %s", 
        ifelse(is.null(rows), ".", rows), 
        ifelse(is.null(cols), ".", cols)
    )

    p <- p + facet_grid(as.formula(facets))
    return(p)
}
