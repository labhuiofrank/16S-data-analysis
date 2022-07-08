#functions that will help with presentation and report writing to easily modify legend and text sizes.

#' Theme for black background presentations
#'
#' @param tsize Base text size
#' @return ggplot2 theme object
#' @export
#' @examples
#' theme_presentation(20)
theme_presentation<- function(tsize=24){
   theme_bw() %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(colour="white"),
      strip.text.y = element_text(colour="white"),
      axis.text.x = element_text(size=tsize-2, angle = 90, colour="white"),
      axis.text.y = element_text(size=tsize-2, colour="white",hjust=1),
      axis.ticks =  element_line(colour = "white"), 
      axis.title.x= element_text(colour="white"),
      axis.title.y= element_text(size=tsize+2, angle = 90,colour="white"),
      panel.background = element_rect(fill="black"), 
      panel.border =element_blank(),  
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      #panel.margin = unit(1.0, "lines"), 
      plot.background = element_rect(fill="black"), 
      plot.title =element_text(colour="white"), 
      #plot.margin = unit(c(1,  1, 1, 1), "lines"),
      legend.background=element_rect(fill='black'),
      legend.title=element_text(colour="white"),
      legend.text=element_text(size=6, colour="white"),
      legend.text.align=0,
      legend.key.size=unit(0.4, 'lines'),
      legend.key = element_rect( fill = 'black'),
      #legend.key.size = unit(c(1, 1), "lines"),
      axis.line.x = element_line(color="white"),
      axis.line.y = element_line(color="white")
    )
}

#' Theme for reports. Starts with theme_grey and then modify some parts
#'
#' @param tsize Base text size
#' @param xaxis Rotation for x-axis text
#' @param xhjust Horizontal text justification (from 0: left to 1: right)
#' @return ggplot2 theme object
#' @export
#' @examples
#' theme_report(15, 45, 1)
theme_report<- function(tsize=18, xaxis=90, xhjust=1) {
  theme_bw() %+replace%
    theme(
      axis.text.x = element_text(size=tsize-6, angle = xaxis,hjust=xhjust, colour="black"),
      axis.text.y = element_text(size=tsize-6,hjust=1, colour="black"),
      axis.title.y= element_text(size=tsize-4, angle = 90, colour="black"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      legend.text=element_text(size=tsize-10),
      legend.text.align=0,
      legend.title.align=0,
      legend.title=element_text(size=tsize-8),
      legend.key.size=unit(0.5, 'lines'))
}

#' Change plot legend to make it smaller
#'
#' @param myPlot ggplot2 object to edit
#' @param pointSize Point size
#' @param textSize Text size
#' @param textSize Spacing between legend elements
#' @export
#' @return ggplot2 object
#' @examples
#' addSmallLegend(p)
addSmallLegend <- function(myPlot, pointSize = 3, textSize = 6, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
