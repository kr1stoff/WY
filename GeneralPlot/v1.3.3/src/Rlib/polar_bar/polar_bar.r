#!/Bio/bin/Rscript-3.6.0
args <- commandArgs(T)

file        <- args[1]
outpre      <- args[2]
color       <- args[3]              # none | red,yellow,green
color_set   <- args[4]              # rainbow | viridis | greyred
title       <- args[5]              # none
format      <- args[6]              # none | log2 log10
header      <- args[7]              # T
show_border <- args[8]              # F
show_label  <- args[9]              # T
show_count  <- args[10]             # T
direction   <- args[11]             # clockwise | anticlockwise
display     <- args[12]				# type1 | type2
width       <- args[13]             # 8
height      <- args[14]             # 8

if( is.na(file) || is.na(outpre) ) {
	stop("[Error]: <file> <outpre> are needed. ")
}
if( is.na(color) )       color = "none"
if( is.na(color_set) )   color_set = "rainbow"
if( is.na(format) )      format = "none"
if( is.na(title) )       title = "none"
if( is.na(header) )      header = T
if( is.na(show_border) ) show_border = F
if( is.na(show_label) )  show_label = T
if( is.na(show_count) )  show_count = T
if( is.na(direction) )   direction = "clockwise"
if( is.na(display) )     display = "type1"
if( is.na(width) )       width = 8
if( is.na(height) )      height = 8
header      = as.logical(header)
show_border = as.logical(show_border)
show_label  = as.logical(show_label)
show_count  = as.logical(show_count)
width       = as.numeric(width)
height      = as.numeric(height)

library(ggplot2)
library(dplyr)

if( "showtext" %in% installed.packages() ) {
    library(showtext)
    showtext_auto(enable = TRUE)
}

mytheme <- theme_void() +
theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
	plot.margin = unit(c(5,5,5,5), "mm")
)

colors <- NULL
if(color != "none") {
	colors <- unlist(strsplit(color, ","))
} else if(color_set == "rainbow") {
	colors <- c("#54778f", "#4EB043", "#E69D2A", "#DD4714", "#A61650")    # rainbow
} else if(color_set == "viridis") {
	colors <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")    # viridis
} else if(color_set == "greyred") {
	colors <- c("#806E6E", "#B7BBC4", "#D8A69F")                          # greyred
} else {
	colors <- c("#54778f", "#4EB043", "#E69D2A", "#DD4714", "#A61650")
}

data <- read.table(file, sep = "\t", head = header)
colnames(data) <- c("name", "count") 
data.len <- length(data$name)

data <- data %>% arrange(count) %>% mutate(
	id = 1:data.len, 
	label = case_when(
		show_label == T && show_count == T ~ paste0(name, "  ", count),
		show_label == T ~ paste0(name),
		show_count == T ~ paste0(count),
		T ~ ""
	), 
	quantile = case_when(
		360 / data.len * id - 360 / data.len / 2 < 90  ~ 1,
		360 / data.len * id - 360 / data.len / 2 < 180 ~ 2,
		360 / data.len * id - 360 / data.len / 2 < 270 ~ 3,
		T ~ 4
	),
	angles = case_when(
		direction == "clockwise" ~ case_when(
			quantile <= 2 ~ -(360 / data.len * id - 180 / data.len) + 90,
			T ~ -(360 / data.len * id - 180 / data.len) + 270
		),
		T ~ case_when(
			quantile <= 2 ~ (360 / data.len * id - 180 / data.len) - 90,
			T ~ (360 / data.len * id - 180 / data.len) - 270
		)
	)
)

if(format == "log2") {
	data$count <- log2(data$count)
	data$count <- as.numeric(sprintf("%0.3f", data$count))
} else if(format == "log10") {
	data$count <- log10(data$count)
	data$count <- as.numeric(sprintf("%0.3f", data$count))
}

textsize <- 4
bordersize <- 0.1
if(data.len <= 10) {
	textsize <- 4
	bordersize <- 0.2
} else if(data.len <= 20) {
	textsize <- 2.8
	bordersize <- 0.1
} else if(data.len <= 30) {
	textsize <- 2.4
	bordersize <- 0.05
} else if(data.len <= 50) {
	textsize <- 2
	bordersize <- 0.02
} else if(data.len <= 100) {
	textsize <- 2 - 1 * (data.len / 100)
	bordersize <- 0.01
} else {
	textsize <- 0.8
	bordersize <- 0.01
}

ymin <- min(data$count) * -0.25
ymax <- max(data$count) * 1.12
if(display != 'type1') {
	# a simple patch for type2 display, nothing could do to deal with toolong labels in type1
	label.str <- as.character(data[data$quantile == 4,]$label)
	label.str.len <- max(nchar(label.str))
	if(label.str.len <= 10) {
		ymax <- max(data$count) * 1.12
	} else if(label.str.len <= 20) {
		ymax <- max(data$count) * 1.32
	} else if(label.str.len <= 30) {
		ymax <- max(data$count) * 1.52
	} else if(label.str.len <= 40) {
		ymax <- max(data$count) * 1.72
	} else {
		ymax <- max(data$count) * 1.82
	}
}

p <- ggplot(data = data, aes(x = id, y = count, label = label)) 
if(show_border == T) {
	p <- p + geom_col(aes(fill = id), width = 1, size = bordersize, color = "#FFFFFF")
} else {
	p <- p + geom_col(aes(fill = id), width = 1, size = 0)
}
p <- p + 
geom_col(aes(y = min(data$count) * 0.5 ), fill = "#FFFFFF", width = 1, alpha = 0.2, size = 0) + 
geom_col(aes(y = min(data$count) * 0.3 ), fill = "#FFFFFF", width = 1, alpha = 0.2, size = 0) + 
scale_fill_gradientn(colors = colors, guide = F) +
scale_y_continuous(limits = c(ymin, ymax)) +
mytheme

if(direction == "clockwise") {
	p <- p + geom_text(
		data = . %>% filter(quantile <= 2), 
		angle = data[data$quantile <= 2,]$angles,
		size = textsize,
		nudge_y = min(data$count) * 0.05,
		vjust = 0.5,
		hjust = 0,
	) + 
	geom_text(
		data = . %>% filter(quantile == 3),
		angle = data[data$quantile == 3,]$angles,
		size = textsize,
		nudge_y = min(data$count) * 0.05,
		vjust = 0.5,
		hjust = 1,
	)
	if(display == "type1") {
		p <- p + geom_text(
			data = . %>% filter(quantile == 4), 
			angle = data[data$quantile == 4,]$angles,
			size = textsize,
			nudge_y = min(data$count) * -0.1,
			vjust = 0.5,
			hjust = 0,
			color = "white",
			alpha = 0.8,
		)
	} else {
		p <- p + geom_text(
			data = . %>% filter(quantile == 4), 
			angle = data[data$quantile == 4,]$angles,
			size = textsize,
			nudge_y = min(data$count) * 0.05,
			vjust = 0.5,
			hjust = 1,
		)
	}
	p <- p + coord_polar("x", direction = 1) 
} else {
	p <- p + geom_text(
		data = . %>% filter(quantile <= 2), 
		angle = data[data$quantile <= 2,]$angles,
		size = textsize,
		nudge_y = min(data$count) * 0.05,
		vjust = 0.5,
		hjust = 1,
	) + 
	geom_text(
		data = . %>% filter(quantile == 3),
		angle = data[data$quantile == 3,]$angles,
		size = textsize,
		nudge_y = min(data$count) * 0.05,
		vjust = 0.5,
		hjust = 0,
	)
	if(display == "type1") {
		p <- p + geom_text(
			data = . %>% filter(quantile == 4), 
			angle = data[data$quantile == 4,]$angles,
			size = textsize,
			nudge_y = min(data$count) * -0.1,
			vjust = 0.5,
			hjust = 1,
			color = "white",
			alpha = 0.8,
		)
	} else {
		p <- p + geom_text(
			data = . %>% filter(quantile == 4), 
			angle = data[data$quantile == 4,]$angles,
			size = textsize,
			nudge_y = min(data$count) * 0.05,
			vjust = 0.5,
			hjust = 0,
		)
	}
	p <- p + coord_polar("x", direction = -1) 
}

if(title != "" && title != "none") p <- p + labs(title = title) 

ggsave(paste0(outpre, ".pdf"), p, width = width, height = height)
#system(paste0("convert -density 1000 ", outpre, ".pdf", " ", outpre, ".png"))

