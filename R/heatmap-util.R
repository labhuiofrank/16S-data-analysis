#' Format sample metadata from a phyloseq object for the ComplexHeatmap package.
#' @param ps Phyloseq object
#' @param columns Columns from sample metadata to display in the heatmap
#' @param extra (list) Extra columns (not in the metadata) to include. 
#' The number of values must match the number of samples in the phyloseq object.
#' @param colors (list) Custom colors to use for any of the metadata in the
#' selected columns or extra metadata. Fields are variable names and values are
#' named vectors (names: levels, values: colors).
#' @return List object with each variable as factor and an extra "col" field
#' with custom colors for each of those variable, as provided in the colors
#' parameter.
#' @export
#' 
#' @examples
#' get_sample_annotations(ps, c("Site_id", "Season"))
#' get_sample_annotations(ps, "Site_id", extra=list(pH=rep(7, nsamples(ps))))
#' get_sample_annotations(
#'   ps, "Site_id", 
#'   colors=list(Site_id=c(S1="red", S2="blue", S3="green"))
#' )
get_sample_annotations <- function(ps, columns=NULL, extra=NULL, colors=NULL) {
  annot <- list()
  
  if(!is.null(extra)) { annot <- extra } # initialize with custom annotations
  if(!is.null(colors)) { annot$col <- colors } # initialize with custom colors
  if(is.null(columns)) {
      warning("No sample metadata provided")
      return(annot)
  }
  
  # fill annotations, colors automatically assigned if not in `col` field
  for (column in columns) {
    if (!column %in% names(annot)) {
      values <- as.factor(phyloseq::get_variable(ps, column))
      annot[[column]] <- values
    }
  }
  return(annot)
}

#' Format OTU metadata from a phyloseq object for the ComplexHeatmap package.
#' @param ps Phyloseq object
#' @param levels Taxonomic levels from tax table to display in the heatmap
#' @param extra (list) Extra columns (not in the tax table) to include. 
#' The number of values must match the number of OTUs in the phyloseq object.
#' @param colors (list) Custom colors to use for any of the annotation variables.
#' Fields are variable names and values are named vectors 
#' (names: factor levels, values: colors).
#' @return List object with each variable as factor and an extra "col" field
#' with custom colors for each of those variable, as provided in the colors
#' parameter, plus custom colors for all of the taxonomic levels provided
#' @export
#' 
#' @examples
#' get_taxa_annotations(ps, c("Phylum", "Class"), extra=list(Nox=rep(F, ntaxa(ps))))
#' get_taxa_annotations(ps, "Phylum", extra=list(Nox=rep(F, ntaxa(ps))), 
#'                      cols=list(Nox=c(T="green", F="red")))
get_taxa_annotations <- function(ps, levels=NULL, extra=NULL, colors=NULL) {
  annot <- list()
  if(!is.null(extra)) { annot <- extra } # initialize with custom annotations
  annot$col <- list()
  if(!is.null(colors)) { annot$col <- colors } # initialize with custom colors
  if(is.null(levels)) {
    warning("No taxonomic levels provided")
    return(annot)
  }
  
  # fill annotations
  for(lvl in levels) {
    if(!lvl %in% names(annot)) {
      values <- as.factor(tax_table(ps)[, lvl])
      annot[[lvl]] <- values
    }
    # set colors
    if(!lvl %in% names(annot$col)) {
      tax_colors <- labhuiofrank.16S::get_taxonomic_palette(values)
      names(tax_colors) <- levels(values)
      annot$col[[lvl]] <- tax_colors
    }
  }
  return(annot)
}
