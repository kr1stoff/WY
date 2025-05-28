#!/Bio/bin/Rscript-3.6.0
#-----------------------------------------------------------------------------
# Program:      scatter_bubble.R
# Author:       dev-group@genedenovo
# Date:         2021-03-23

library(ggplot2, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(reshape2, warn.conflicts = F)

if( "showtext" %in% installed.packages() ) {
	library(showtext, warn.conflicts = F)
	showtext_auto(enable = TRUE)
}

#-----------------------------------------------------------------------------
args <- commandArgs(T)
exp_file     <- args[1]
pv_file      <- args[2]
outpref      <- args[3]
rowname_file <- args[4]			# none
colname_file <- args[5]			# none
exp_scale    <- args[6]			# none | log2 | log10 | -log2 | -log10
pv_scale     <- args[7]			# none | -log2 | -log10 | log2 | log10
dot_minsize  <- args[8]			# 1
dot_maxsize  <- args[9]			# 6
colors       <- args[10]		# none | "black,yellow,blue,red"
transparency <- args[11]		# 1
exp_lab      <- args[12]		# none
pv_lab       <- args[13]		# none
x_lab        <- args[14]		# none
y_lab        <- args[15]		# none
main_title   <- args[16]		# none

if( is.na(exp_file) | is.na(pv_file) | is.na(outpref) ) {
	stop("[Error]: <exp_file> <pv_file> <outprefix> are needed. ")
}

if( is.na(rowname_file) ) rowname_file <- "none"
if( is.na(colname_file) ) colname_file <- "none"
if( is.na(exp_scale) )    exp_scale <- "none"
if( is.na(pv_scale) )     pv_scale <- "none"
if( is.na(dot_minsize) )  dot_minsize <- 1
if( is.na(dot_maxsize) )  dot_maxsize <- 6
if( is.na(colors) )       colors <- "none"
if( is.na(transparency) ) colors <- 1
if( is.na(exp_lab) )      exp_lab <- "none"
if( is.na(pv_lab) )       pv_lab <- "none"
if( is.na(x_lab) )        x_lab <- "none"
if( is.na(y_lab) )        y_lab <- "none"
if( is.na(main_title) )   main_title <- "none"

transparency <- as.numeric(transparency)
dot_minsize  <- as.numeric(dot_minsize)
dot_maxsize  <- as.numeric(dot_maxsize)

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
	theme1 <- theme_bw() + 
	theme(
		panel.grid = element_blank(),
		panel.border = element_rect(color = "#000000", size = 0.8),
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
	color1 <- c("black", "blue", "yellow", "red")
}

scatter_bubble <- function(
exp_file = NULL, pv_file = NULL, outpref = NULL, rowname_file = "none", colname_file = "none", 
exp_scale = "none", pv_scale = "none", dot_minsize = 1, dot_maxsize = 6, colors = "none", 
exp_lab = "none", pv_lab = "none", x_lab = "none", y_lab = "none", main_title = "none") {

	exp_data <- read.table(exp_file, head = T, sep = "\t", check.names = F, quote = "")
	pv_data  <- read.table(pv_file,  head = T, sep = "\t", check.names = F, quote = "")
	colnames(exp_data)[1] <- "col_idx"
	colnames(pv_data)[1]  <- "col_idx"
	rownames(exp_data) <- exp_data$col_idx
	rownames(pv_data)  <- pv_data$col_idx

	if(rowname_file != "none") {
		if(file.exists(rowname_file)) {
			select_row <- readLines(rowname_file)
			exp_data   <- exp_data %>% filter(col_idx %in% select_row)
			pv_data    <- pv_data  %>% filter(col_idx %in% select_row)
		} else {
			runMessage("W", "the row name file isnot exists.")
		}
	}
	if(colname_file != "none") {
		if(file.exists(colname_file)) {
			select_col <- readLines(colname_file)
			exp_data   <- exp_data %>% select(c("col_idx", select_col))
			pv_data    <- pv_data  %>% select(c("col_idx", select_col))
		} else {
			runMessage("W", "the column name file isnot exists.")
		}
	}

	exp_rowns <- rownames(exp_data)
	exp_colns <- colnames(exp_data)
	pv_rowns  <- rownames(pv_data)
	pv_colns  <- colnames(pv_data)
	inter_rowns <- intersect(exp_rowns, pv_rowns)
	inter_colns <- intersect(exp_colns, pv_colns)

	if(length(exp_rowns) != length(pv_rowns)) {
		runMessage("W", "the row length isnot equal.")
	}
	if(length(exp_colns) != length(pv_colns)) {
		runMessage("W", "the row length isnot equal.")
	}
	if(length(inter_rowns) == 0 || length(inter_colns) == 0) {
		runMessage("E", "no cross id in row/column, check the files.")
	}

	exp_data <- exp_data[inter_rowns, inter_colns]
	pv_data  <- pv_data [inter_rowns, inter_colns]
	exp_data$col_idx <- factor(exp_data$col_idx, levels = rev(unique(exp_data$col_idx)))
	pv_data$col_idx  <- factor(pv_data$col_idx, levels = rev(unique(pv_data$col_idx)))
	data     <- left_join(melt(exp_data), melt(pv_data), by = c("col_idx", "variable"))
	data$value.x[data$value.x == 0] <- mean(data$value.x) * 0.01
	data$value.x[data$value.x == 0] <- 0.001
	data$value.y[data$value.y == 0] <- 0.00009

	if(exp_scale == "log2") {
		data$value.x <- log2(data$value.x)
	} else if(exp_scale == "log10") {
		data$value.x <- log10(data$value.x)
	} else if(exp_scale == "-log2") {
		data$value.x <- -log2(data$value.x)
	} else if(exp_scale == "-log10") {
		data$value.x <- -log10(data$value.x)
	}
	if(pv_scale == "-log2") {
		data$value.y <- -log2(data$value.y)
	} else if(pv_scale == "-log10") {
		data$value.y <- -log10(data$value.y)
	} else if(pv_scale == "log2") {
		data$value.y <- log2(data$value.y)
	} else if(pv_scale == "log10") {
		data$value.y <- log10(data$value.y)
	}
	if(dot_minsize > dot_maxsize) {
		runMessage("W", "the scatter min size is greater than the scatter max size.")
		dot_minsize = dot_maxsize
	}
	if(colors != "none") {
		colors <- unlist(strsplit(colors, ","))
	} else {
		colors <- default_color()
	}

	width  <- FigWidth(inter_colns)
	height <- FigWidth(inter_rowns)
	axis.size.x <- AxisTextSize(data$variable, width = width)
	axis.size.y <- AxisTextSize(data$col_idx, width = height)


	p <- ggplot(data, aes(x = variable, y = col_idx, color = value.x, size = value.y)) + 
	geom_point(alpha = transparency) + 
	scale_color_gradientn(colors = colors) + 
	scale_size_continuous(range = c(dot_minsize, dot_maxsize)) + 
	default_theme() + 
	theme(
	  axis.text.x = element_text(size = axis.size.x), 
	  axis.text.y = element_text(size = axis.size.y)) + 
	xlab_angle(unique(data$variable), width = width)

	if( exp_lab    == "none" )  exp_lab    <- "value.1"
	if( pv_lab     == "none" )  pv_lab     <- "value.2"
	if( x_lab      == "none" )  x_lab      <- NULL
	if( y_lab      == "none" )  y_lab      <- NULL
	if( main_title == "none" )  main_title <- NULL
	p <- p + labs(color = exp_lab, size = pv_lab, x = x_lab, y = y_lab, title = main_title)

	outdir <- dirname(outpref)
	outname <- basename(outpref)
	if(dir.exists(outdir) != TRUE) dir.create(outdir, showWarnings = T, recursive = T)
	SaveFig(p, outpref, width, height)
}

scatter_bubble(exp_file, pv_file, outpref, rowname_file, colname_file, exp_scale, pv_scale, 
dot_minsize, dot_maxsize, colors, exp_lab, pv_lab, x_lab, y_lab, main_title)

