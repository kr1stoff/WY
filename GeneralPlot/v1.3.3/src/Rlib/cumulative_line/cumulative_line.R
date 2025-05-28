#!/Bio/bin/Rscript-3.6.0
#-----------------------------------------------------------------------------
# Program:      cumulative_line.R
# Author:       dev-group@genedenovo
# Date:         2021-05-25

library(ggplot2, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(reshape2, warn.conflicts = F)
library(scales)
library(RColorBrewer)

if( "showtext" %in% installed.packages() ) {
    library(showtext, warn.conflicts = F)
    showtext_auto(enable = TRUE)
}

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
		quit(status = 1)
    }
}

SaveFig = function(fig = fig, outpfx = outpfx, width = 8, height = 8) {
    ggsave(file = paste0(outpfx, ".pdf"), fig, width = width, height = height, limitsize = FALSE)
    #ggsave(file = paste0(outpfx, ".png"), fig, width = width, height = height, limitsize = FALSE)
    #system(paste0("convert -density 300 ", outpfx, ".pdf", " ", outpfx, ".png"))
}

default_color <- function() {
    color1 <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
    "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
}

default_theme <- function() {
    theme1 <- theme_classic() +
    theme(
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "#000000", size = 0.5),
        axis.ticks.length = unit(0.1, 'cm'),
        axis.text = element_text(color = "#000000", size = 10),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 1, vjust = 0.5),
        axis.title = element_text(color = "#000000", size = 14),
        axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
        axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
		legend.title = element_blank(),
		legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 16),
        plot.margin = unit(c(5,5,5,5), "mm")
    )
    theme1
}

#-----------------------------------------------------------------------------
cumulative_exp <- function ( exp_file = NULL, outpref = NULL, group_file = "none", 
filter0 = FALSE, scale_data = "none", namefix = "none", color_set = "none", 
title = "none", x_lab = "none", y_lab = "none", line_size = 1, line_alpha = 0.9, 
width = 8, height = 7, out_format = "pic_data" ) {

	if( is.null(outpref) ) {
		runMessage( "E", "the outpref isnot set." )
	}
	outdir  <- dirname(outpref)
	outname <- basename(outpref)

	if( out_format != "pic" & out_format != "data" & out_format != "pic_data" ) {
		runMessage( "E", "the out_format should be [ pic | data | pic_data ]." )
	}

	## read exp file
	if( is.null(exp_file) == FALSE & file.exists(exp_file) == TRUE ) {
		data <- read.table( exp_file, head = T, sep = "\t", check.names = F, quote = "" )
	} else {
		runMessage( "E", "the exp matrix file isnot exists." )
	}
	if( namefix != "none" ) {
		getcol <- grep( namefix, colnames(data) )
		if(length(getcol) > 0) {
			data <- data[, c(1, getcol)]
			colnames(data) <- gsub( namefix, "", colnames(data) )
		} else {
			runMessage( "E", paste0( "no data found with namefix : [", namefix, "]" ) )
		}
	}

	colnames(data)[1] <- "id"
	data.tb <- melt(data)
	colnames(data.tb) <- c("id", "sample", "value")
	data.tb$group <- data.tb$sample

	## read group file
	if ( is.null(group_file) == FALSE & group_file != "none" ) {
		if( file.exists(group_file) == TRUE ) {
			group <- read.table( group_file, head = F, sep = "\t", check.names = F, quote = "" )
			colnames(group) <- c( "sample", "group" )
			rownames(group) <- group[, 1]
			group$group <- factor( group$group, levels = unique(group$group) )
			group <- group[ as.vector(data.tb$sample), ]
			data.tb$group <- group$group
		} else {
			runMessage( "W", "no group file." )
		}
	}

	# group data
	data.tb2 <- data.tb %>% group_by(group, id) %>% summarize(value = mean(value))

	## filter zero data
	if( as.logical(filter0) ) {
		data.tb2 <- data.tb2[data.tb2$value > 0, ]
	}

	## format the data
	data.tb2$fmt_value <- data.tb2$value
	if( scale_data == "log2+1" ) {
		data.tb2$fmt_value <- log2(data.tb2$fmt_value + 1)
	} else if( scale_data == "log10+1" ) {
		data.tb2$fmt_value <- log10(data.tb2$fmt_value + 1)
	} else if( scale_data == "log2" ) {
		data.tb2$fmt_value[data.tb2$fmt_value == 0] = 0.001
		data.tb2$fmt_value <- log2(data.tb2$fmt_value)
	} else if( scale_data == "log10" ) {
		data.tb2$fmt_value[data.tb2$fmt_value == 0] = 0.001
		data.tb2$fmt_value <- log10(data.tb2$fmt_value)
	}

	## stat
	data.tb2.st <- data.tb2 %>% group_by(group) %>% summarize(num = n(), sum = sum(value), 
	mean = mean(value), fmt_sum = sum(fmt_value), fmt_mean = mean(fmt_value))

	## calculate cumulative frequence
	data.tb2$freq <- 0
	for ( i in data.tb2.st$group ) {
		data.tb2[data.tb2$group == i, ] <- data.tb2[data.tb2$group == i, ] %>% mutate(
			freq = case_when(
				group == i ~ value / data.tb2.st[data.tb2.st$group == i, ]$sum
			)
		)
	}
	data.tb2 <- data.tb2 %>% group_by(group) %>% arrange(freq) %>% 
	mutate(cumsum_freq = cumsum(freq))
	data.tb3 <- data.tb2 %>% group_by(group, fmt_value) %>% 
	summarize(mean_cumsum_freq = mean(cumsum_freq)) %>% arrange(fmt_value) 

	## plot settings
	if( color_set != "none" ) {
		color_set <- unlist(strsplit(color_set, ","))
	} else {
		color_set <- default_color()
	}

	if( length(unique(data.tb2$group)) > 1 ) {
		if( length(unique(data.tb2$group)) <= length(color_set) ) {
			color_set <- color_set[1:length(unique(data.tb2$group))]
		} else {
			runMessage("W", "the group num is greater than the color num. ")
			color_set <- colorRampPalette(color_set)(length(unique(data.tb2$group)))
		}
	} else {
		color_set <- color_set[1]
	}

	if( x_lab == "none" )  x_lab <- NULL
	if( y_lab == "none" )  y_lab <- NULL
	if( title == "none" )  title <- NULL

	if( out_format == "pic" | out_format == "pic_data" ) {
		p <- ggplot( data.tb3, aes(x = fmt_value, y = mean_cumsum_freq, color = group) ) + 
		geom_line( size = line_size, alpha = line_alpha ) + 
		scale_color_manual( values = color_set ) + 
		labs( x = x_lab, y = y_lab, title = title ) + 
		default_theme()

		if(dir.exists(outdir) != TRUE) dir.create(outdir, showWarnings = T, recursive = T)
		SaveFig(p, outpref, width, height)
	}

	if( out_format == "data" | out_format == "pic_data" ) {
		write.table(data.tb3, file = paste0(outpref, ".data.xls"), quote = F, sep = "\t", 
		row.names = F, col.names = T)
	}

}

#-----------------------------------------------------------------------------
args <- commandArgs(T)
exp_file    <- args[1]
outpref     <- args[2]
group_file  <- args[3]        # none
filter0     <- args[4]        # FALSE | TRUE
scale_data  <- args[5]        # none | log2 | log10 
namefix     <- args[6]        # none | _fpkm, _count ...
color_set   <- args[7]        # none
title       <- args[8]        # none
x_lab       <- args[9]        # none
y_lab       <- args[10]       # none
line_size   <- args[11]       # 1
line_alpha  <- args[12]       # 0.9
width       <- args[13]       # 7
height      <- args[14]       # 7
out_format  <- args[15]       # pic_data | pic | data

if( is.na(exp_file) | is.na(outpref) ) {
    stop("[Error]: <file> <outprefix> are needed. ")
}

if( is.na(group_file) )  group_file <- "none"
if( is.na(filter0) )     filter0 <- FALSE
if( is.na(scale_data) )  scale_data <- "none"
if( is.na(namefix) )     namefix <- "none"
if( is.na(color_set) )   color_set <- "none"
if( is.na(title) )       title <- "none"
if( is.na(x_lab) )       x_lab <- "none"
if( is.na(y_lab) )       y_lab <- "none"
if( is.na(line_size) )   line_size <- 1
if( is.na(line_alpha) )  line_alpha <- 0.9
if( is.na(width) )       width <- 8
if( is.na(height) )      height <- 7
if( is.na(out_format) )  out_format <- "pic_data"

filter0    <- as.logical(filter0)
line_size  <- as.numeric(line_size)
line_alpha <- as.numeric(line_alpha)
width      <- as.numeric(width)
height     <- as.numeric(height)

#-----------------------------------------------------------------------------
cumulative_exp( exp_file, outpref, group_file, filter0, scale_data, namefix, 
color_set, title, x_lab, y_lab, line_size, line_alpha, width, height, out_format )


