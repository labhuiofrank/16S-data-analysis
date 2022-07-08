library(speedyseq)
library(ggplot2)

#' Custom palette for large categorical groups
#' Cycles through the values every 50 colors
#'
#' @param n Number of colors to generate
#' @return list of n colors (hex color codes)
#' @examples
#' get_palette(10)
get_palette <- function(n) {
    palette <- c("#45cc65","#8e3ebc","#62ba3f","#9164dc","#a75355","#91ba33","#cd74e5",
                 "#33a150","#b8319c","#4dcc90","#e75dc1","#538d2a","#586cd8","#dca831",
                 "#9049a2","#b2b33a","#d57ccb","#80bc6c","#e54486","#4fcfc3","#df4b40",
                 "#4db9df","#de632b","#5e91d3","#d68330","#4363a5","#9c8026","#ae99e4",
                 "#60751c","#7961a8","#3d864a","#af3a77","#66b78b","#c3364e","#2ba198",
                 "#aa4b28","#35896b","#e2798a","#1a6447","#e390bf","#306a3c","#a56395",
                 "#4d6d30","#924869","#aab56e","#e38f6b","#666020","#d1a965","#966433",
                 "#868949")
  while (n > length(palette)) {
    palette <- rep(palette, 1+floor(n/50))
  }
  return(palette[1:n])
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
#' @examples
#' taxa_barplot(ps, x="sample_type", taxrank="Phylum", min_relabund=0.01)
#' taxa_barplot(ps, x="sample_type", y="relabund", taxrank="Phylum", rows="Season")
taxa_barplot <- function(ps, x=NULL, y="Abundance", taxrank="Class",
                         rows=NULL, cols=NULL,
                         min_relabund=0.01, return_df=FALSE) {
    # Group taxa by {tax_rank}
    ps <- speedyseq::tax_glom(ps, taxrank)
    # Group samples for each combination of factors in [x, rows, cols]
    all_cols <- sapply(c(x, rows, cols), function(v) phyloseq::get_variable(ps, v))
    sample_data(ps)$combined <- apply(all_cols, 1, paste0, collapse="_")
    ps <- speedyseq::merge_samples2(ps, "combined")
    # Transform phyloseq object to data frame and compute relative abundances
    data <- speedyseq::psmelt(ps) %>% group_by(combined) %>%
        mutate(relabund=Abundance/sum(Abundance))
    # Replace low abundance species with filler "< xx %"
    filler <- sprintf('Other(<%.f%%)', 100*min_relabund)
    data[data$relabund < min_relabund, taxrank] <- filler
    
    # Get the taxa ordering (by total relabund)
    sorted_taxa <- data %>% group_by_at(taxrank) %>%
        summarize(score=sum(relabund)) %>%
        arrange(score) %>% pull(1)

    # Order the taxa and put the filler at the end
    data <- data %>% mutate_at(taxrank, ~ factor(., sorted_taxa))
    if (filler %in% sorted_taxa) {
        data <- data %>% mutate_at(taxrank, ~ relevel(., filler))
    }

    # Return
    if(return_df) {
        return(data)
    } 
    p <- ggplot(data=data, aes_string(x=x, y=y, fill=taxrank)) +
        geom_bar(stat="identity", position="stack", color="black", size=0.5, width=0.7) +
        scale_fill_manual(values=get_palette(nlevels(data %>% pull(taxrank)))) +
        guides(fill=guide_legend(reverse=T))

    # Facet arguments, need to handle null values
    facets <- sprintf(
        "%s ~ %s", 
        ifelse(is.null(rows), ".", rows), 
        ifelse(is.null(cols), ".", cols)
    )

    p <- p + facet_grid(as.formula(facets))
    return(p)
}
