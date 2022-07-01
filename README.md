# 16S data analysis helper scripts

This repository is a collection of useful scripts and function for the analysis 16S sequencing data.

## **R --- R scripts**

### For helping to format sequencing data into phyloseq objects
- **loading.R**: Helper script to appropriately load and read all the files needed into a phyloseq 
	- load\_shared(file) 
	- load\_counttable(file)
	- load\_constaxonomy(file)
- util.R

### Analyzing alpha diversity
- alpha\_diversity.R
- barplot.R
- taxa\_barplot\_kiana.R

### Analyzing beta diversity
- beta\_diversity.R
- ordination.R
- run\_nmds.R
- DESeq2.R

### Subsetting data
- subset\_barplot.R
- single\_taxa\_barplot.R
- single\_taxa_boxplot.R

### Presentation and Reporting Aesthetics
- theme\_presentation.R
- consistent\_color\_pallete.R