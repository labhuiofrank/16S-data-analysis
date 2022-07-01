## Script for loading the microbiome data into phyloseq
## Includes path to all of the necessary files
## Need to download first the labhuiofrank helper scripts in: https://github.com/labhuiofrank/16S-data-analysis.git

source("labhuiofrank-16S-util/R/loading.R")
library(phyloseq)

# path to data
abundance_file <- "data/OTUs.100.rep.count_table"
taxonomy_file <- "data/OTUs.100.cons.taxonomy"
metadata_file <- "data/metadata.csv"

# Loader function to import for each analyses
load_microbiome <- function() {
    abundance <- load_counttable(abundance_file)
    taxonomy <- load_constaxonomy(as.matrix(taxonomy_file))
    metadata <- read.table(metadata_file, sep=",", row.names=1, header=T) 
    metadata <- metadata %>% 
      select(-c(Water_Level_cm, Turbidity)) %>%
      na.omit() %>%
      mutate(Ahupuaa_Location=factor(Ahupuaa_Location, c("Before_Kawainui", "Kawainui", "After_Kawainui", "Kailua_Bay")))

    ps <- phyloseq(
        otu_table(abundance, taxa_are_rows=T),
        tax_table(as.matrix(taxonomy)),
        sample_data(metadata)
    )
    return(ps)
}

#' Custom palette function for large categorical groups
#' Repeats the palette if more than 50 groups
get_palette <- function(n) {
  palette <- c("#45cc65","#8e3ebc","#62ba3f","#9164dc","#a75355","#91ba33","#cd74e5","#33a150","#b8319c","#4dcc90","#e75dc1","#538d2a","#586cd8","#dca831","#9049a2","#b2b33a","#d57ccb","#80bc6c","#e54486","#4fcfc3","#df4b40","#4db9df","#de632b","#5e91d3","#d68330","#4363a5","#9c8026","#ae99e4","#60751c","#7961a8","#3d864a","#af3a77","#66b78b","#c3364e","#2ba198","#aa4b28","#35896b","#e2798a","#1a6447","#e390bf","#306a3c","#a56395","#4d6d30","#924869","#aab56e","#e38f6b","#666020","#d1a965","#966433","#868949")
  if (n > 50) {
    palette <- rep(palette, 1+floor(n/50))[1:n]
  }
  return(palette)
}

options(ggplot2.discrete.colour=RColorBrewer::brewer.pal(9, "Set1"))
options(ggplot2.discrete.fill=RColorBrewer::brewer.pal(9, "Set1"))


