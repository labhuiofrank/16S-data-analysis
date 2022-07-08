source("util.R")
library("ggplot2")
library("speedyseq") # phyloseq extension
library("vegan")

# Global variable
# What are we focusing on?
target <- "Site_Type"

## =====================##
##     Load our data    ##
## =====================##
mb <- load_microbiome()

## ======================##
##  Normalize samples    ##
## ======================##

# Compute sample counts
sample_counts <- sample_sums(mb)

# quickly
hist(log10(sample_counts), n=30)
threshold <- 4000

# Same, but a bit prettier
sample_sizes <- data.frame(x=sample_counts)

ggplot(sample_sizes, aes(x=x)) +
    geom_histogram(aes(y=..density..), alpha=0.5) + # bars
    geom_density(fill="blue", alpha=0.3) +          # curve
    geom_vline(xintercept=threshold, color="red", linetype=2) + # red dotted line to see the threshold
    scale_x_log10()

# Subsampling
mb_subsampled <- rarefy_even_depth(mb, threshold, rngseed=123, replace=FALSE)

## ==============================##
## Unsupervised ordination: NMDS ##
## ==============================##

nmds <- ordinate(mb_subsampled, method="NMDS", trymax=300, parallel=8)
plot_ordination(mb_subsampled, nmds, color=target) + 
  geom_point(size=3) +
  stat_ellipse(aes_string(x="NMDS1", y="NMDS2", fill=target), 
               alpha=0.1, show.legend=FALSE, level=0.5, geom="polygon") 

## ============================##
## Supervised ordination: CCA  ##
## ============================##

cca <- ordinate(mb_subsampled, method="CCA", formula=~ DO_percent + Temp_C + Sal_ppt + pH)

# A bit more complex: add the environmental variables as arrows
arrow_coords <- scores(cca, display = "bp") %>% data.frame %>% rownames_to_column("labels")
arrowhead = arrow(length = unit(0.05, "npc"))

plot_ordination(mb_subsampled, cca, type="Sample", color=target, shape="Ahupuaa_Location") + geom_point(size=3) +  
  geom_segment(aes(x=0, xend=CCA1, y=0, yend=CCA2, color=NULL, shape=NULL), arrow = arrowhead, data=arrow_coords) +
  geom_text(aes(x=1.2 * CCA1, y=1.2 * CCA2, label=labels, color=NULL, shape=NULL), data=arrow_coords)

# Other option: custom function from labhuiofrank github
source("labhuiofrank-16S-util/R/ordination.R")
cca2 <- run_ordination(mb_subsampled, "CCA", formula=~ DO_percent + Temp_C + Sal_ppt + pH)
plot_ordination_with_arrows(cca2$data, cca2$arrows, "CCA1", "CCA2", fill=target) 

