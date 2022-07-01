#function will plot the relative abundances of only the subseted taxa that you are interested in.  You must provide a .csv of a list of Genus to subset, and a palette specific for that subset.  You can make a specific palette with the consistent color palette script.
library(ggpubr)
library(ggplot2)

subset_barplot <- function(mg, tax_lvl="Genus", groups, thresh=1e-4, normalize=TRUE,nrows_legend=30, 
                           subset=pathogens, palette=pathogen_palette) {
  # Format data to long format
  data_full <- speedyseq::psmelt(
    speedyseq::tax_glom(mg, tax_lvl)
  )[, c(tax_lvl, groups, "Sample", "Abundance")]
  
  data_full1 <- data_full %>% 
    group_by_at(groups) %>% 
    mutate(Total_Abundance=sum(Abundance))
  
  data <- subset(data_full1, Genus %in% rownames(subset))
  
  filler <- sprintf('Other(<%.g%%)', 100*thresh)
  
  data1 <- data %>% 
    group_by_at(groups) %>% 
    mutate(relabund=Abundance/(Total_Abundance)) %>%
    group_by_at(c(groups, tax_lvl)) %>%
    mutate(score=sum(relabund)) %>% 
    arrange(score)%>%
    mutate_at(tax_lvl, ~ ifelse(score<thresh, filler, .))
  
  library(data.table)  
  setDT(data1)
  data2 <- data1[, .(relabund = sum(relabund)), by = c(tax_lvl, groups)]
  
  data3 <- data1[, .(relabund = sum(relabund)), by = c(tax_lvl)]
  
  # Set palette
  
  taxa1 <- palette %>% add_row(Genus=filler, value="grey")
  taxa1 <- as.tibble(taxa1)
  pal <- left_join(data3, taxa1, by= tax_lvl)
  pal <- pal[order(relabund)]
  palette <- pal$value
  names <- pal[,1]
  names(names)<- c("taxa")
  names(palette) <- names$taxa
  palette
  
  y <- ifelse(normalize, 'relabund', 'count')
  title <- paste("Abundance by", groups, sep=" ") 
  p <- data2 %>%
    ggplot(aes_string(x=groups[1], y=y, fill=tax_lvl)) +
    geom_bar(stat="identity", position="stack", color='black', size=0.3, width=0.7) +
    scale_fill_manual(values=palette) +
    scale_y_continuous(labels=scales::comma) +
    ylab('Relative Abundance') + xlab('') +
    theme_linedraw() +
    ggtitle(title)+
    guides(fill=guide_legend(nrow=nrows_legend, reverse=T)) +
    theme(axis.text=element_text(size=8),
          axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
          axis.title=element_text(size=(10)),                          
          legend.title = element_text(size = (8)), 
          axis.line=element_line("gray25"),
          panel.grid.major.x=element_blank(),
          panel.grid=element_line("gray25"),
          legend.text.align=0,
          legend.text=element_text(size=6),
          legend.key.size=unit(0.5, 'lines'))
  
  if(length(groups) > 1) {
    if(length(groups) == 2) {
      form_str <- sprintf("~ %s", groups[2])
    } else {
      form_str <- sprintf("%s ~ %s", groups[2], groups[3])
    }
    p <- p + facet_grid(as.formula(form_str), scales="free_x")
  }
  
  return(p)
}