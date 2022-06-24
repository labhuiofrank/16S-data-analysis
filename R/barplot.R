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
  if (n > 50) {
    palette <- rep(palette, 1+floor(n/50))[1:n]
  }
  return(palette)
}

#' Taxonomic stacked barplot of phyloseq object.
#'
#' @param ps Phyloseq object to plot
#' @param x Categorical variable in metadata, displayed in x-axis
#' @param taxrank Taxonomic rank to gather taxa (stacked for each level of x)
#' @param min_relabund Below this threshold, replace taxa name to "< x%" to
#' make the plot for readable and reduce the number of colors
#' @return ggplot object (geom_bar)
#' @examples
#' taxa_barplot(ps, x="sample_type", taxrank="Phylum", min_relabund=0.01)
taxa_barplot <- function(ps, x=NULL, taxrank=NULL, min_relabund=0.01) {
    # Sum taxa by {tax_rank}
    ps <- speedyseq::tax_glom(ps, taxrank)
    # Sum samples for each level of {x}
    ps <- speedyseq::merge_samples2(ps, x)
    # Convert to relative abundance
    ps <- transform_sample_counts(ps, function(x) x/sum(x))
    
    # Transform phyloseq object to data frame
    data <- speedyseq::psmelt(ps)[, c(taxrank, x, "Abundance")]

    # Replace low abundance species with filler "< xx %"
    filler <- sprintf('Other(<%.f%%)', 100*min_relabund)
    data[data$Abundance < min_relabund, taxrank] <- filler
    
    # Get the taxa ordering (by total relabund)
    sorted_taxa <- data %>% group_by_at(taxrank) %>%
        summarize(score=sum(Abundance)) %>%
        arrange(score) %>% pull(1)

    # Order the taxa and put the filler at the end
    data <- data %>% mutate_at(taxrank, ~ factor(., sorted_taxa))
    if (filler %in% sorted_taxa) {
        data <- data %>% mutate_at(taxrank, ~ relevel(., filler))
    }

    # Display
    ggplot(data=data, aes_string(x=x, y="Abundance", fill=taxrank)) +
        geom_bar(stat="identity", position="stack", color="black", size=0.5, width=0.7) +
        scale_fill_manual(values=get_palette(nlevels(data[, taxrank]))) +
        ylab('Relative abundance') + 
        guides(fill=guide_legend(reverse=T))
}
