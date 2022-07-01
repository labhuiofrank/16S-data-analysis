get_palette <- function(n) {
  palette <- c( "#45cc65","#8e3ebc","#62ba3f","#9164dc","#a75355","#91ba33","#cd74e5","#33a150","#b8319c","#4dcc90","#e75dc1","#538d2a","#586cd8","#dca831","#9049a2","#b2b33a","#d57ccb","#80bc6c","#e54486","#4fcfc3","#df4b40","#4db9df","#de632b","#5e91d3","#d68330","#4363a5","#9c8026","#ae99e4","#60751c","#7961a8","#3d864a","#af3a77","#66b78b","#c3364e","#2ba198","#aa4b28","#35896b","#e2798a","#1a6447","#e390bf","#306a3c","#a56395","#4d6d30","#924869","#aab56e","#e38f6b","#666020","#d1a965","#966433","#868949","#069668", "#c7eea8", "#5d99aa", "#81f4fe", "#788ce0", "#955ccd", "#f75ef0", "#b95581", "#f89ade", "#dad2de", "#7efc9d", "#4ba40b", "#6f7d43", "#cfdf34", "#e13224", "#faa566", "#bb5c2a", "#88fe0e", "#e21c7a", "#ab8791", "#399283", "#3eeaef", "#154e56", "#40e18c", "#39970e", "#b8e27d", "#683c00", "#ebb9c7", "#f4327e", "#8d102b", "#f76015", "#bf711e", "#eac328", "#1f3e9e", "#7b7ec0", "#aedbf0", "#8d30ba", "#d179f8", "#3f16f9", "#fa2beb", "#5a3e4f", "#c7628e", "#50e316", "#3f4c08", "#938073","#FF0000", "#FF0F00", "#FF1F00", "#FF2E00", "#FF3D00", "#FF4D00", "#FF5C00", "#FF6B00", "#FF7A00", "#FF8A00", "#FF9900", "#FFA800", "#FFB800", "#FFC700", "#FFD600", "#FFE500" , "#FFF500" , "#FAFF00" , "#EBFF00" , "#DBFF00" , "#CCFF00" , "#BDFF00", "#ADFF00" , "#9EFF00" , "#8FFF00" , "#80FF00" , "#70FF00" , "#61FF00" , "#52FF00" , "#42FF00" , "#33FF00" , "#24FF00" , "#14FF00", "#05FF00" , "#00FF0A" , "#00FF1A" , "#00FF29" , "#00FF38" , "#00FF47" , "#00FF57" , "#00FF66" , "#00FF75" , "#00FF85" , "#00FF94",  "#00FFA3" , "#00FFB3" , "#00FFC2" , "#00FFD1" , "#00FFE0" , "#00FFF0" , "#00FFFF" , "#00F0FF" , "#00E0FF" , "#00D1FF" , "#00C2FF",  "#00B2FF" , "#00A3FF" , "#0094FF" , "#0085FF" , "#0075FF" , "#0066FF" , "#0057FF" , "#0047FF" , "#0038FF" , "#0029FF" , "#0019FF",  "#000AFF" , "#0500FF" , "#1400FF" , "#2400FF" , "#3300FF" , "#4200FF" , "#5200FF" , "#6100FF" , "#7000FF" , "#8000FF" , "#8F00FF",  "#9E00FF" , "#AD00FF" , "#BD00FF" , "#CC00FF" , "#DB00FF" , "#EB00FF" , "#FA00FF" , "#FF00F5" , "#FF00E6" , "#FF00D6" , "#FF00C7",  "#FF00B8" , "#FF00A8" , "#FF0099" , "#FF008A" , "#FF007A" , "#FF006B" , "#FF005C" , "#FF004C" , "#FF003D" , "#FF002E" , "#FF001F","#FF000F")
  if (n > 195) {
    palette <- rep(palette, 1+floor(n/195))[1:n]
  } else {palette <- palette[1:n]}
  return(palette)
}


taxa <- read.table("/Users/klfrank/Dropbox/Post-Grad/Bioinformatics/R-Microbiome-Statistics/Functions/Function_Scripts/Consistent_Pallet/taxa_barfile_color_palette.csv", header=T, row.names=1, sep=",")

#PALETTE FOR TAXA BARPLOT  
n <- nrow(taxa)
col <- get_palette(n)
taxa <- as_tibble(taxa)
col_pal <- as_tibble(col)
taxa_palette <- taxa %>% add_column(col_pal)

#PALETTE FOR PATHOGENS
pathogens <- read.table("/Users/klfrank/Dropbox/Post-Grad/Bioinformatics/R-Microbiome-Statistics/Functions/Function_Scripts/Consistent_Pallet/pathogens_genus.csv", header=T, row.names=1, sep=",")
n <- nrow(pathogens)
col <- get_palette(n)
path <- rownames(pathogens)
path <- as.tibble(path)
names(path) <- "Genus"
col_pal <- as.tibble(col)
pathogen_palette <- path %>% add_column(col_pal)

#PALETTE FOR PLANT PATHOGENS
plant_pathogens <- read.table("/Users/klfrank/Dropbox/Post-Grad/Bioinformatics/R-Microbiome-Statistics/Functions/Function_Scripts/Consistent_Pallet/Plant_Pathogens.csv", header=T, row.names=1, sep=",")
n <- nrow(plant_pathogens)
col <- get_palette(n)
path <- rownames(plant_pathogens)
path <- as.tibble(path)
names(path) <- "Genus"
col_pal <- as.tibble(col)
plant_pathogen_palette <- path %>% add_column(col_pal)