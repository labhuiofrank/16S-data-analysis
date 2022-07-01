# 16S data analysis helper scripts

This repository is a collection of useful scripts and function for the analysis 16S sequencing data.

## R 
**R scripts**

### For helping to format sequencing data into phyloseq objects
- **loading.R**: Helper script to appropriately load and read all the files needed into a phyloseq 
	- load\_shared(file) : Load OTU abundance shared file from mothur (abundance.100.shared)
	- load\_counttable(file) : Load OTU abundance count table file from mothur (OTUs.100.rep.count_table)
	- load\_constaxonomy(file) : Load consensus taxonomy from mothur (annotations.100.cons.taxonomy)
- **util.R**: Script for loading the microbiome data into phyloseq
	- Need to set path to all necessary files and adjust metadata levels
	- load\_microbiome() : Loader function to import for each analyses
	- get\_palette(n) : Custom palette function for large categorical groups,  repeats the palette if more than 50 groups

### Analyzing alpha diversity
- **alpha\_diversity.R** : See the Rmd\_examples\_explanaitions folder for more detailed explanation
	- Load data sourcing the loading.R and util.R scripts
	- Preprocess data via subsetting, merging or applying a prevalence filter
	- Generate taxanomic barplots by target/variable of interest
	- Generate alpha diversity metrics and display in different forms
		- Chao1: Estimator of the richness that looks at rare species to get closer to the true sample richness
		- Shannon: Provides information about both richness and evenness.
		- Inverse Simpson: “This simply equals true diversity of order 2, i.e. the effective number of types that is obtained when the weighted arithmetic mean is used to quantify average proportional abundance of types in the dataset of interest. The index is also used as a measure of the effective number of parties.”
- **barplot.R** : cedrics pretty version will plot taxonomic stacked barplot at different levels of phyloseq object
- **taxa\_barplot\_kiana.R** : kianas version will plot taxonomic stacked barplot at different levels of phyloseq object

### Analyzing beta diversity
- **beta\_diversity.R** : See the Rmd\_examples\_explanaitions folder for more detailed explanation
	- Load data sourcing the loading.R and util.R scripts
	- Normalize samples by rarefying to set threshold
	- Ordnation figures (NMDS, CCA) by sourcing ordination.R
- **ordination.R** : 
	- run\_ordnation(ps, method = "NMDS", ...) : Wrapper around phyloseq's ordinate function, adds the metadata to the component scores
	- plot\_ordination\_with\_arrows(data, arrows, x="CCA1", y="CCA2", fill="black", size=3, arrow_rescale=0.5, ellipse=TRUE, ...)
	- run\_adonis(ps, meta_vars, ...): Wrapper around vegan::adonis
- **DESeq2.R** : 
	- geom_mean(x) : Geometric mean function for vectors with zeros
	- run_deseq2(ps, design, gm=FALSE, pseudocount=0) : Wrapper around DESeq functions
	- get_results_for_contrast(diagdds, factor_name, ref, level)
	- get_all_results(diagdds, variable, comparisons) : Get result table from DESeq2 object between all levels vs reference
		

### Visualizing Subset data as part of the whole
- **subset\_barplot.R** : will plot a taxonomic stacked bar plot of only the subsetted taxa
	-subset\_barplot(ps, tax\_lvl="Genus", groups, thresh=1e-4, normalize=TRUE,nrows\_legend=30, subset=pathogens, palette=pathogen\_palette)
- **single\_taxa\_barplot.R** : will plot a stacked bar plot at ASV level of a single genus relative abundance
	- single\_taxa\_barplot(ps, tax\_lvl="Genus", groups, thresh=1e-4, normalize=TRUE, nrows\_legend=15, single\_tax="Staphylococcus")
- **single\_taxa\_boxplot.R** : will plot a boxplot of relative abundance of single genus from all samples 
	- single\_taxa\_boxplot(ps, tax\_lvl="Genus", groups, normalize=TRUE, nrows\_legend=2, single_tax="Staphylococcus" )

### Presentation and Reporting Aesthetics 
- **theme\_presentation.R**: functions that will help with presentation and report writing to easily modify legend and text sizes.
	- theme\_presentation(tsize=24): add + to any ggplot to format for black background presentation
	- theme\_report(tsize=18, xaxis=90, xhjust=1) : add + to any ggplot to format for a written report
	- addSmallLegend(myPlot, pointSize = 3, textSize = 6, spaceLegend = 0.1) : will shrink the legend which is especially helpful for taxa barplots
- **consistent\_color\_pallete.R** : formats custom palettes when you want to compare different sets of data but have a consistent coloring of taxa
	- Need to do this to set palettes for subset\_barplot
	- get_palette(n): Custom color palette for large groups, repeats palette when n > 195

## Rmd\_examples\_explained 
**Annotated R-markdown scripts and examples**

## databases 
**.csv databases for subsetting or setting color palettes**
- **taxa\_barfile\_color\_palette.csv** : List of all taxa at every level, based on Silva v 138 taxa annotations
- **FAPROTAX1.2.4_to_Silvav138.csv** : Mapping file to Silva v 138 with all of FAPROTAX annotations and metabolisms, can use this file to get at all types of metabolisms
- **function_N_bygenus.csv** : : List of putative Nitrogen cyclers from FAPROTAX version 1.2.4 at genus level
- **nfix.csv** : List of putative nitrogen ficerss from FAPROTAX version 1.2.4 at genus level
- **DRNA.csv** : List of putative organisms capable of DNRA from FAPROTAX version 1.2.4 at genus level
- **function_N_S_bygenus.csv**: List of putative Nitrogen and Sulfur cyclers from FAPROTAX version 1.2.4 at genus level
- **function_S_bygenus.csv** : List of putative Sulfur cyclers from FAPROTAX version 1.2.4 at genus level
- **methanecycling.csv** : List of putative methane cyclers from FAPROTAX version 1.2.4 at genus level
- **PATRIC\_human\_pathogens\_lineage.csv** : PATRIC database to species/strain level
- **pathogens_genus.csv** : List of putative pathogens from PATRIC database at genus level
- **plant_pathogens.csv** : List of putative plant pathogens from FAPROTAX version 1.2.4 at genus level
