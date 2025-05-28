#!/Bio/bin/Rscript-3.6.0
#-----------------------------------------------------------------------------
# Program:      lollipop.R
# Author:       dev-group@genedenovo
# Date:         2021-03-23

library(ggplot2, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(reshape2, warn.conflicts = F)
library(RColorBrewer)

if( "showtext" %in% installed.packages() ) {
	library(showtext, warn.conflicts = F)
	showtext_auto(enable = TRUE)
}

#-----------------------------------------------------------------------------
args <- commandArgs(T)
file         <- args[1]
outpref      <- args[2]
group_file   <- args[3]			# none
scale        <- args[4]			# none | log2 | log10 | -log2 | -log10
dot_size     <- args[5]			# 4
line_size    <- args[6]         # 1
dot_color    <- args[7]			# none
line_color   <- args[8]			# grey50
dot_alpha    <- args[9]			# 1
is_label     <- args[10]		# F
is_baseline  <- args[11]		# F
is_flip      <- args[12]		# F
is_sort      <- args[13]		# F
x_lab        <- args[14]		# none
y_lab        <- args[15]		# none
title        <- args[16]		# none
lgd_lab      <- args[17]		# none

if( is.na(file) | is.na(outpref) ) {
	stop("[Error]: <file> <outprefix> are needed. ")
}

if( is.na(group_file) )  group_file <- "none"
if( is.na(scale) )       scale <- "none"
if( is.na(dot_size) )    dot_size <- 4
if( is.na(line_size) )   line_size <- 1
if( is.na(dot_color) )   dot_color <- "none"
if( is.na(line_color) )  line_color <- "grey50"
if( is.na(dot_alpha) )   dot_alpha <- 1
if( is.na(is_label) )    is_label <- FALSE
if( is.na(is_baseline) ) is_baseline <- FALSE
if( is.na(is_flip) )	 is_flip <- FALSE
if( is.na(is_sort) )     is_sort <- FALSE
if( is.na(x_lab) )       x_lab <- "none"
if( is.na(y_lab) )       y_lab <- "none"
if( is.na(title) )       title <- "none"
if( is.na(lgd_lab) )	 lgd_lab <- "none"

dot_size    <- as.numeric(dot_size)
line_size   <- as.numeric(line_size)
dot_alpha   <- as.numeric(dot_alpha)
is_label    <- as.logical(is_label)
is_baseline <- as.logical(is_baseline)
is_flip     <- as.logical(is_flip)
is_sort     <- as.logical(is_sort)

#-----------------------------------------------------------------------------
runMessage <- function(tag = "R", msgs = NULL) {
	if(tag == "R" | tag == "Running") {
		message(paste0("[", as.character(Sys.time()) ,"] [Running]: ", msgs))
	} 
	else if(tag == "W" | tag == "Warning") {
		message(paste0("[", as.character(Sys.time()) ,"] [Warning]: ", msgs))
	}
	else if(tag == "E" | tag == "Error") {
		message(paste0("[", as.character(Sys.time()) ,"] [Error]: ", msgs))
	}
}

FigWidth = function(sample.name, width = 8) {
    sample.num = length(sample.name)
    if(sample.num >= 800) {
        width = 30
    } else if(sample.num >= 500) {
        width = 25
    } else if(sample.num >= 100) {
        width = switch(as.integer(sample.num / 100), 11, 14, 17, 20)
    } else {
        if(sample.num >= 50) {
            width = 10
        } else if(sample.num >= 20) {
            width = 9
        }
    }

    width
}

AxisTextSize <- function(value, width = 8) {
    name <- as.character(unique(value))
    name.len <- sum(nchar(name))
    name.meanlen <- name.len / length(name)

    if(name.meanlen <= 10) {
        size = 12
    } else if(name.meanlen <= 20) {
        size = 11
    } else {
        size = 10
    }

    name.num <- length(name)
    if(name.num >= 30) {
        size = width * 0.6 * 72 / name.num
    }

    size
}

GeomTextSize = function(value, sample.num = 10, width = 8, fontsize = 4, legend = T) {
    name = as.character(value)
    name.len = sum(nchar(name))
    name.maxlen = max(nchar(name))
    name.meanlen = name.len / length(name)

    name.font.len = fontsize * 1 / 72 * name.len * 1
    name.max.font.len = fontsize * 1 / 72 * name.maxlen * 1
    draw.width = ifelse(legend, width - 3, width - 2)

    if(name.max.font.len >= draw.width / sample.num * 0.8) {
        fontsize = draw.width * 72 / sample.num / name.maxlen * 0.5
    } else if(name.max.font.len >= draw.width / (sample.num + 2) * 0.8) {    # recheck
        fontsize = draw.width * 72 / sample.num / name.maxlen * 0.5
    }

    fontsize
}

xlab_angle = function(sample_name, width = 7, fontsize = 14) {
    name = as.character(unique(sample_name))
    name_len = sum(nchar(name))
    name_mean = name_len / length(name)

    name_font_len = fontsize * 1/72 * name_len * 0.6
    cut = c(0.6, 2) * width

    if(name_font_len < cut[1]) {
        angle = 0
        hjust = 0.5
        vjust = 0
    } else if(name_font_len < cut[2]) {
        if(name_mean > 15){
            angle = 45
            hjust = 1
            vjust = 1
        } else {
            angle = 45
            hjust = 1
            vjust = 1
        }
    } else {
        angle = 90
        hjust = 1
        vjust = 0.5 
    }

    theme(axis.text.x = element_text(angle=angle, hjust = hjust, vjust = vjust))
}

SaveFig = function(fig = fig, outpfx = outpfx, width = 7, height = 7) {
    ggsave(file = paste0(outpfx, ".pdf"), fig, width = width, height = height, limitsize = FALSE)
	#ggsave(file = paste0(outpfx, ".png"), fig, width = width, height = height, limitsize = FALSE)
	#system(paste0("convert -density 300 ", outpfx, ".pdf", " ", outpfx, ".png"))
}

default_theme <- function() {
	theme1 <- theme_classic() + 
	theme(
		#panel.grid = element_blank(),
		#panel.grid.minor = element_blank(),
		panel.grid.major = element_line(linetype = "twodash", size = 0.5, color = "grey75"),
		panel.grid.major.x = element_blank(),
		#panel.border = element_rect(color = "#000000", size = 0.8),
		axis.ticks = element_line(color = "#000000", size = 0.5),
		axis.ticks.length = unit(0.1, 'cm'),
		axis.text = element_text(color = "#000000", size = 10),
		axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
		axis.text.y = element_text(hjust = 1, vjust = 0.5),
		axis.title = element_text(color = "#000000", size = 14),
		axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
		axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
		plot.title = element_text(hjust = 0.5, size = 16),
		plot.margin = unit(c(5,5,5,5), "mm")
	)
	theme1
}

default_color <- function() {
	color1 <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", 
	"#91D1C2", "#DC0000", "#7E6148", "#B09C85")
}

lollipop <- function(
file = NULL, outpref = NULL, group_file = "none", scale = "none", dot_size = 4, line_size = 1, 
dot_color = "none", line_color = "grey50", dot_alpha = 1, is_label = F, is_baseline = F, 
is_flip = F, is_sort = F, x_lab = "none", y_lab = "none", title = "none", lgd_lab = "none") {

	data <- read.table(file, head = T, sep = "\t", check.names = F, quote = "")
	coln <- colnames(data)
	colnames(data) <- c("index", "value")

	if( group_file != "none" ) {
		if( file.exists(group_file) ) {
			group <- read.table(group_file, head = F, sep = "\t", check.names = F, quote = "")
			colnames(group) <- c("index", "group")
			data <- left_join(data, group, by = "index")
			data$group <- as.character(data$group)
		} else {
			runMessage("W", "the group file isnot exists.")
			data <- data %>% mutate(group = "1")
			data$group <- as.character(data$group)
		}
	} else {
		data <- data %>% mutate(group = "1")
		data$group <- as.character(data$group)
	}
	data$index <- factor(data$index, levels = unique(data$index))

	if( is_sort == T ) {
		data <- data %>% arrange(desc(group), value)
		data$index <- factor(data$index, levels = unique(data$index))
	}

	if( scale == "log2" ) {
		data$value <- log2(data$value)
	} else if( scale == "log10" ) {
		data$value <- log10(data$value)
	} else if( scale == "-log2" ) {
		data$value <- -log2(data$value)
	} else if( scale == "-log10" ) {
		data$value <- -log10(data$value)
	}

	if( dot_color != "none" ) {
		dot_color <- unlist(strsplit(dot_color, ","))
	} else {
		dot_color <- default_color()
	}
	if( line_color == "none" ) {
		line_color = NA
	} else {
		line_color = line_color
	}
	if( length(unique(data$group)) > 1 ) {
		if( length(unique(data$group)) <= length(dot_color) ) {
			dot_color <- dot_color[1:length(unique(data$group))]
		} else {
			runMessage("W", "the group num is greater than the color num. ")
			dot_color <- colorRampPalette(dot_color)(length(unique(data$group)))
		}
	} else {
		dot_color <- dot_color[1]
	}

	y_min <- min(data$value) * 1.12
	y_max <- max(data$value) * 1.12
	if( y_min > 0 ) y_min = 0

	width <- FigWidth(data$index)
	height <- 7
	axis.size <- AxisTextSize(data$index, width = width)
	label.size <- GeomTextSize(data$value, length(data$value), width = width)

	p <- ggplot(data, aes(x = index, y = value, group = group))
	p <- p + geom_segment(aes(x = index, xend = index, y = 0, yend = value), 
		size = line_size, color = line_color)

	if(length(dot_color) > 1) {
		p <- p + geom_point(aes(color = group), size = dot_size, alpha = dot_alpha) + 
		scale_color_manual(values = dot_color)
	} else {
		p <- p + geom_point(color = dot_color, size = dot_size, alpha = dot_alpha)
	}

	if( is_label == T ) {
		p <- p + geom_text(aes(label = value), size = label.size)
	}
	if(is_baseline == T & min(data$value) < 0 & max(data$value) > 0) {
		p <- p + geom_hline(yintercept = 0, size = line_size, color = line_color, linetype = 4)
	}

	p <- p + scale_y_continuous(expand = c(0, 0), limits = c(y_min, y_max)) 
	
	if( lgd_lab    == "none" )  lgd_lab <- NULL
	if( x_lab      == "none" )  x_lab <- NULL
	if( y_lab      == "none" )  y_lab <- NULL
	if( title      == "none" )  title <- NULL
	p <- p + labs(color = lgd_lab, x = x_lab, y = y_lab, title = title) + 
	default_theme() 

	if( is_flip == T ) {
		p <- p + coord_flip() + 
		theme(axis.text.y = element_text(size = axis.size), panel.grid.major.y = element_blank())
		tmp <- height
		height <- width
		width <- tmp
	} else {
		height <- height - 1
		p <- p + theme(axis.text.x = element_text(size = axis.size)) +
		xlab_angle(unique(data$index), width = width)
	}

	outdir <- dirname(outpref)
	outname <- basename(outpref)
	if(dir.exists(outdir) != TRUE) dir.create(outdir, showWarnings = T, recursive = T)
	SaveFig(p, outpref, width, height)
}

lollipop(file, outpref, group_file, scale, dot_size, line_size, dot_color, line_color, 
dot_alpha, is_label, is_baseline, is_flip, is_sort, x_lab, y_lab, title, lgd_lab)

