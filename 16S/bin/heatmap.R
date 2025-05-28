#!/usr/bin/env Rscript
options(warn = -1)
suppressPackageStartupMessages({
  library(funr)
  library(stringr)
  library(yaml)
  BIN = dirname(sys.script())
  f_software <- str_c(BIN, "../config/software.yml", sep = '/')
  software <- yaml.load_file(f_software)
  library(optparse)
  library(RColorBrewer)
  library(pheatmap)
  library(logging)
})
basicConfig()


option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "input [default %default]"),
  make_option(c("-l", "--list"), type = "character", default = NULL,
              help = "gene list [default %default]"),
  make_option(c("-m", "--map"), type = "character", default = NULL,
              help = "groups data [default %default]"),
  make_option(c("-a", "--annrow"), type = "character", default = NULL,
              help = "annotation data for row [default %default]"),
  make_option(c("-p", "--palette"), type = "character", default = 'Spectral',
              help = "color palette : red,black,green RdYlBu Blues Spectral ... [default %default]"),
  make_option(c("-c", "--cluster"), type = "character", default = 'row',
              help = "cluster [default %default]"),
  make_option(c("-s", "--scale"), type = "character", default = 'z-score',
              help = "scale : z-score|min-max|log|none [default %default]"),
  make_option(c("-w", "--cw"), type = "numeric", default = NULL,
              help = "cell width [default 16.18]"),
  make_option(c("-o", "--outpre"), type = "character", default = NULL,
              help = "output directory [default %default]"),
  make_option(c("-t", "--title"), type = "character", default = NA,
              help = "graph title [default %default]")
)
opts <- parse_args(OptionParser(usage = "%prog [options]", option_list = option_list), positional_arguments = F)
logdebug(opts)

args <- commandArgs(T)
if (is.null(opts$input) | is.null(opts$outpre)) stop("Rscript heatmap.v1.pl -i <file> -c row -o heatmap\n")

## read data
data <- read.table(opts$input, sep = "\t", header = T, row.names = 1, check.names = F, quote = "", na = '')
if (!is.null(opts$list)) {
  glist <- read.table(opts$list, head = F, comment = '', stringsAsFactors = F)[, 1]
  data <- data[glist,]
}

## read group data
ann_col <- NA
if (!is.null(opts$map)) {
  ann_col <- read.table(opts$map, sep = "\t", check = F, row = 1, head = T, comment = '')
  ann_col[1] <- factor(ann_col[, 1], levels = unique(ann_col[, 1]))
  logdebug(str(ann_col))
  data <- data[, rownames(ann_col)]
}
ann_row <- NA
if (!is.null(opts$annrow)) {
  ann_row <- read.table(opts$annrow, sep = "\t", check = F, row = 1, head = T, comment = '')
  ann_row[1] <- factor(ann_row[, 1], levels = unique(ann_row[, 1]))
  logdebug(str(ann_row))
  data <- data[rownames(ann_row),]
}

## filter
mt <- data[apply(data, 1, sd) != 0,]

## scale
if (opts$scale == 'z-score') {
  mt <- t(scale(t(mt)))
} else if (opts$scale == 'min-max') {
  min <- apply(mt, 1, min)
  max <- apply(mt, 1, max)
  mt <- (mt - min) / (max - min)
} else if (opts$scale == 'log') {
  mt <- log(mt)
} else if (opts$scale == 'none') {
  logdebug("No scale use")
} else {
  logerror('Scale function not support!(choose : z-score,min-max,log,none)\n')
  q()
}
## set colors
#  red,black,green RdYlBu Blues Spectral
cols <- NULL
if (grepl(',', opts$palette)) {
  cols <- colorRampPalette(rev(strsplit(opts$palette, ',')[[1]]))(50)
} else {
  cols <- colorRampPalette(rev(brewer.pal(7, opts$palette)))(50)
}

## set cluster
cluster <- c(T, F)
if (opts$cluster == 'none') {
  cluster <- c(F, F)
} else if (opts$cluster == 'col') {
  cluster <- c(F, T)
} else if (opts$cluster == 'all') {
  cluster <- c(T, T)
}

## cell width/height
cw <- 10 * 1.618
ch <- 10

if (!is.null(opts$cw)) {
  cw <- opts$cw
} else if (ncol(mt) == nrow(mt)) {
  cw <- ifelse(ncol(mt) > 20, 8, 160 / ncol(mt))
  ch <- cw
} else if (nrow(mt) / ncol(mt) > 6) {
  cw <- nrow(mt) / 6 / ncol(mt) * 16.18
  logdebug('cw :', cw, '\n')
}

png <- paste0(opts$outpre, '.png')
pdf <- paste0(opts$outpre, '.pdf')

if (cw > 50) {
  pdf(file = pdf, width = 12, height = 20)
  pheatmap(mt, scale = 'none', cluster_row = cluster[1], cluster_col = cluster[2], border_color = NA,
           color = cols, annotation_col = ann_col, annotation_row = ann_row, main = opts$title)
  invisible(dev.off())
  png(file = png, width = 12, height = 20, units = 'in', res = 300)
  pheatmap(mt, scale = 'none', cluster_row = cluster[1], cluster_col = cluster[2], border_color = NA,
           color = cols, annotation_col = ann_col, annotation_row = ann_row, main = opts$title)
  invisible(dev.off())
  unlink('Rplots.pdf')

} else {
  ## draw heatmap
  pdf(file = pdf, width = 12, height = 20)
  pheatmap(mt, scale = 'none', cluster_row = cluster[1], cluster_col = cluster[2],
           color = cols, filename = pdf, cellwidth = cw, cellheight = ch,
           annotation_col = ann_col, annotation_row = ann_row, border_color = NA,
           main = opts$title)
  invisible(dev.off())
  png(file = png, width = 12, height = 20, units = 'in', res = 300)
  pheatmap(mt, scale = 'none', cluster_row = cluster[1], cluster_col = cluster[2],
           color = cols, filename = png, cellwidth = cw, cellheight = ch,
           annotation_col = ann_col, annotation_row = ann_row, border_color = NA,
           main = opts$title)
  invisible(dev.off())
  unlink('Rplots.pdf')
}

