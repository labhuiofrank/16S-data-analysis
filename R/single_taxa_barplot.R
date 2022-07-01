#function to plot one single taxa as a barplot

single_taxa_barplot <- function(mg, tax_lvl="Genus", groups, thresh=1e-4, normalize=TRUE, 
                                nrows_legend=15, single_tax="Staphylococcus") {
  
  # Format data to long format
  data_full7 <-speedyseq::psmelt(mg)[, c(tax_lvl, groups, "Sample", "OTU", "Abundance")]
  
  data_full1 <- data_full7 %>% 
    group_by_at(groups) %>% 
    mutate(Total_Abundance=sum(Abundance))
  
  data4 <- subset(data_full1, Genus %in% single_tax) 
  
  filler <- sprintf('Other(<%.g%%)', 100*thresh)
  
  data1 <- data4 %>% 
    group_by_at(groups) %>% 
    mutate(relabund=Abundance/(Total_Abundance)) %>%
    group_by_at(c(groups, "OTU")) %>%
    mutate(score=sum(relabund)) %>% 
    arrange(score)%>%
    mutate_at(tax_lvl, ~ ifelse(score<thresh, filler, .))
  
  library(data.table)  
  setDT(data1)
  data2 <- data1[, .(relabund = sum(relabund)), by = c(tax_lvl,"OTU", groups)]
  data3 <- data1[, .(relabund = sum(relabund)), by = c("OTU")]
  
  #Set palette
  ntaxa <- nrow(data3)
  palette <- rainbow(ntaxa)
  y <- ifelse(normalize, 'relabund', 'count')
  
  title1 <- paste("Total Abundance of", single_tax, "by", groups,  sep=" ")
  b <- data2 %>%
    ggplot(aes_string(x=groups[1], y=y, fill="OTU")) +
    geom_bar(stat="identity", position="stack", color='black', size=0.3, width=0.7) +
    scale_fill_manual(values="forestgreen") +
    scale_y_continuous(labels=comma) +
    ylab('Relative Abundance') + xlab('') +
    theme_linedraw() +
    ggtitle(title1)+
    guides(fill=guide_legend(nrow=nrows_legend, reverse=T)) +
    theme(axis.text.x=element_text(size=6, angle=90, vjust=0.5, hjust=1),
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
    b <- b + facet_grid(as.formula(form_str))
  }
  
  b
  list <- list(b, ntaxa)
  return(b)
}