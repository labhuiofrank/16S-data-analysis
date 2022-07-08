source("util.R")
library("ggplot2")
library("speedyseq") # phyloseq extension

# Global variable
# What are we focusing on?
target <- "Site_Type"

## =====================##
##     Load our data    ##
## =====================##

mb <- load_microbiome()
mb <- merge_samples2(mb, "Sample_no")

## =====================##
##    Preprocess data   ##
## =====================##

mb_per_type <- merge_samples2(mb, target)
mb_per_type <- transform_sample_counts(mb_per_type, function(x) x/sum(x))

# Keep ASVs that make at least 1% of the total abundance
mb_core <- filter_taxa(mb_per_type, function(x) any(x > 0.01), TRUE)
mb_core <- transform_sample_counts(mb_core, function(x) x/sum(x))

## =================================##
## Display in barplot with phyloseq ##
## =================================##

plot_bar(mb_core, target, fill="Class") +
  geom_bar(aes(fill=Order), stat="identity", position="stack") +
  scale_fill_manual(values=get_palette(ntaxa(mb_core)))

## ========================##
## Alpha diversity metrics ##
## ========================##

indices <- c("Chao1", "Shannon", "InvSimpson")

# One dot for each sample in group
plot_richness(mb, color=target, measures=indices)
# Or as a boxplot per site
plot_richness(mb, x=target, color=target, measures=indices) + 
  geom_boxplot()

# A bit prettier
plot_richness(mb, x=target, color=target, measures=indices) + 
  geom_boxplot(aes_string(fill=target), color="black") +
  scale_fill_brewer(palette="Set1") + scale_color_brewer(palette="Set1")

