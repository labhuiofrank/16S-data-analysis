#function to plot one single taxa as a boxplot

single_taxa_boxplot <- function(mg, tax_lvl="Genus", groups, normalize=TRUE, 
                                nrows_legend=2, single_tax="Staphylococcus" ) {
  
  # Format data to long format
  data_full <- speedyseq::psmelt(
    speedyseq::tax_glom(mg, tax_lvl)
  )[, c(tax_lvl, groups, "Sample", "OTU", "Abundance")]
  
  data_full1 <- data_full %>% 
    group_by_at(groups) %>% 
    mutate(Total_Abundance=sum(Abundance))
  
  data4 <- subset(data_full1, Genus %in% single_tax) 
  
  data1 <- data4 %>% 
    group_by_at(groups) %>% 
    mutate(relabund=Abundance/(Total_Abundance)) %>%
    group_by_at(c(groups, tax_lvl)) %>%
    mutate(score=sum(relabund)) %>% 
    arrange(score)
  
  library(data.table)  
  setDT(data1)
  data4 <- data1[, .(relabund = mean(relabund)), by = c(tax_lvl,"Sample", groups)]
  
  y <- ifelse(normalize, 'relabund', 'count')
  title2 <- paste("Average Sample Abundance of", single_tax, "by", groups,  sep=" ")
  
  s <- data4 %>%
    ggplot(aes_string(x=groups[1], y='relabund', fill=tax_lvl))+
    geom_boxplot() +
    geom_point(aes_string(x=groups[1], y='relabund'), size=1) +
    scale_y_continuous(labels=comma) +
    ylab('Relative Abundance') + xlab('') +
    theme_linedraw() +
    ggtitle(title2) +  
    guides(fill=guide_legend(nrow=nrows_legend, reverse=T)) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
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
    s <- s + facet_grid(as.formula(form_str))
  }
  s
  return(s)
}