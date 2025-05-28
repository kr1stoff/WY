#!/Bio/Bin/pipeline/System_Programs/Rscript.3.6.3

args <- commandArgs(trailingOnly = TRUE)

file         <- args[1]
outpre       <- args[2]
color        <- args[3]
color_set    <- args[4]
transparency <- args[5]
band_width   <- args[6]
band_knotpos <- args[7]
title        <- args[8]
xlab         <- args[9]
ylab         <- args[10]
title_size   <- args[11]
axis_size    <- args[12]
label_size   <- args[13]
label_color  <- args[14]
stratum_color <- args[15]
alluvium_color <- args[16]
fill_column  <- args[17]

if( is.na(file) || is.na(outpre) ) {
	stop("[Error]: <file> <outpre> are needed. ")
}

if( is.na(color) )       	 color = "none"
if( is.na(color_set) )   	 color_set = "set1"
if( is.na(transparency) )    transparency = 0.7
if( is.na(band_width) )   	 band_width = 0.33
if( is.na(band_knotpos) )    band_knotpos = 0.25
if( is.na(title) )   		 title = "none"
if( is.na(xlab) )   		 xlab = "none"
if( is.na(ylab) )   		 ylab = "none"
if( is.na(title_size) )   	 title_size = 16
if( is.na(axis_size) )       axis_size = 12
if( is.na(label_size) )   	 label_size = 4
if( is.na(label_color) )     label_color = "black"
if( is.na(stratum_color) )   stratum_color = "black"
if( is.na(alluvium_color) )  alluvium_color = NA
if( is.na(fill_column) )     fill_column = "default"

transparency  <- as.numeric(transparency)
band_width    <- as.numeric(band_width)
band_knotpos  <- as.numeric(band_knotpos)
title_size    <- as.numeric(title_size)
axis_size     <- as.numeric(axis_size)
label_size    <- as.numeric(label_size)

library(ggalluvial)
library(ggplot2)
library(scales)
library(RColorBrewer)

if( "showtext" %in% installed.packages() ) {
    library(showtext)
    showtext_auto(enable = TRUE)
}

mytheme <- theme_bw() + 
theme(
	panel.grid = element_blank(),
	panel.border = element_blank(), 

	axis.title = element_text(color = "#000000", size = title_size - 1),
	axis.text = element_text(color = "#000000", size = axis_size),
	axis.text.y = element_blank(),
	axis.line = element_blank(),
	axis.ticks = element_blank(),

	legend.position = "none",

	plot.title = element_text(color = "#000000", size = title_size, hjust = 0.5),
	plot.margin = unit(c(5,5,5,5), "mm")
)

colors <- NULL
if( color != "none" && color != "" ) {
	colors <- unlist(strsplit(color, ","))
} else if( color_set == "set1" ) {
	colors <- c("yellow", "red", "green")
} else if( color_set == "set2" ) {
	colors <- c("blue", "green", "pink")
} else if( color_set == "set3" ) {
	colors <- c("yellow", "green", "blue")
} else if( color_set == "set4" ) {
	colors <- c("steelblue", "green", "lightgreen")
} else {
	colors <- c("yellow", "red", "green")
}

data <- read.table(file, sep = "\t", head = T, check.names = F)
data_raw <- data

if(fill_column == "default") {
	data <- cbind(
		data[1:(length(data) - 1)], 
		#.tmp_class = data[, length(data) - 1], 
		.tmp_class = data[, 1], 
		data[length(data)]
	)
} else if(is.na(as.numeric(fill_column)) == FALSE) {
	tmp_col = as.numeric(fill_column)
	data <- cbind(
		data[1:(length(data) - 1)],
		.tmp_class = data[, tmp_col],
		data[length(data)]
	)
} else {
	stop("[Error]: <fill_column> is not right format. ")
}

data_long <- to_lodes_form(
	data, 
	key = "x", 
	value = "stratum", 
	id = "alluvium", 
	axes = 1:(length(data) - 2)
)
is_ok <- is_lodes_form(data_long, key = "x", value = "stratum", id = "alluvium", weight = "Freq")
if( is_ok != TRUE ) {
	stop("[Error]: the data cannot be converted, check it.")
}

if( length(colors) != length(unique(data_long$.temp_class)) ) {
	colors <- colorRampPalette(colors)(length(unique(data_long$.tmp_class)))
}

width  <- 10
height <- 8
row_len <- length(unique(data_raw[,1]))
for (i in 1:(length(data_raw) - 1)) {
	if(row_len < length(unique(data_raw[,i]))) {
		row_len <- length(unique(data_raw[,i]))
	}
}
if( length(data_raw) > 3 ) {
	width <- width + 2 * (length(data_raw) - 3)
}
if( row_len > 20 ) {
	height <- height + 1.25 * ceiling((row_len - 20) / 5)
}

p <- ggplot(data = data_long, 
  aes(x = x, y = Freq, stratum = stratum, alluvium = alluvium, fill = .tmp_class, label = stratum)) + 
  geom_alluvium(alpha = transparency, width = band_width, knot.pos = band_knotpos, 
  color = alluvium_color) +
  geom_stratum(alpha = transparency, width = band_width, color = stratum_color) + 
  geom_text(stat = "stratum", size = label_size, color = label_color) +
  labs(x = NULL, y = NULL) + 
  scale_fill_manual(values = colors) + 
  scale_x_discrete(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) + 
  mytheme

if( title != "none" && title != "" )    { p <- p + ggtitle(title) }
if( xlab  != "none" && xlab  != "" )    { p <- p + xlab(xlab) }
if( ylab  != "none" && ylab  != "" )    { p <- p + ylab(ylab) }

ggsave(paste0(outpre, ".pdf"), p, width = width, height = height)
#ggsave(paste0(outpre, ".png"), p, width = width, height = height)
#system(paste0("convert -density 300 ", outpre, ".pdf", " ", outpre, ".png"))

