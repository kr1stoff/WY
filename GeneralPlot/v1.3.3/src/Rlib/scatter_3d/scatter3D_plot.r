#!/Bio/bin/Rscript
args <- commandArgs(trailingOnly = FALSE)
Bin <- dirname(sub("--file=", "", args[grep("--file=", args)]))
source(paste(Bin, "addgrids3d.r", sep = "/"))

library(RColorBrewer)
library(scatterplot3d)
library(getopt)

command <- matrix(c(
	'file'    , 'f', 1, 'character', 'set the input file path', 
	'outdir'  , 'o', 1, 'character', 'set the output dir',
	'outname' , 'n', 1, 'character', 'set the output name',
	'xcol'    , 'x', 2, 'integer'  , 'set the x column, default is 2',
	'ycol'    , 'y', 2, 'integer'  , 'set the y column, default is 3', 
	'zcol'    , 'z', 2, 'integer'  , 'set the z column, default is 4',
	'sample'  , 's', 2, 'integer'  , 'set the sample column, default is 1',
	'group'   , 'u', 2, 'integer'  , 'set the group column, default is "none"',
	'title'   , 'm', 2, 'character', 'set the title, default is "none"',
	'xlab'    , 'a', 2, 'character', 'set the xlab, default is the column name of x',
	'ylab'    , 'b', 2, 'character', 'set the ylab, default is the column name of y',
	'zlab'    , 'c', 2, 'character', 'set the zlab, default is the column name of z',
	'pcex'    , 'p', 2, 'integer'  , 'point size [1]',
	'tcex'    , 't', 2, 'integer'  , 'point text size [0.8]',
	'vjust'   , 'v', 2, 'double'   , 'point label position adjsut [0.5]',
	'color'   , 'r', 2, 'character', 'set the color string by ","',
	'isgrid'  , 'd', 2, 'logical'  , 'add grid background or not, default is "T"',
	'islabel' , 'l', 2, 'logical'  , 'add label to point, default is "F"',
    'help'    , 'h', 0, "logical"  , 'show this message.'

), byrow = TRUE, ncol = 5)
opt = getopt(command)

if ( !is.null(opt$help) || is.null(opt$file) || is.null(opt$outdir) || is.null(opt$outname) ) {
	cat(paste(getopt(command, usage = T), "\n"))
	q(status = 1)
}

if ( is.null(opt$x      ) ) { opt$x       = 2      }
if ( is.null(opt$y      ) ) { opt$y       = 3      }
if ( is.null(opt$z      ) ) { opt$z       = 4      }
if ( is.null(opt$sample ) ) { opt$sample  = 1      }
if ( is.null(opt$group  ) ) { opt$group   = 'none' }
if ( is.null(opt$color  ) ) { opt$color   = 'none' }
if ( is.null(opt$pcex   ) ) { opt$pcex    = 1      }
if ( is.null(opt$tcex   ) ) { opt$tcex    = 0.8    }
if ( is.null(opt$vjust  ) ) { opt$vjust   = 0.5    }
if ( is.null(opt$isgrid ) ) { opt$isgrid  = T      }
if ( is.null(opt$islabel) ) { opt$islabel = F      }
if ( is.null(opt$width  ) ) { opt$width   = 8      }
if ( is.null(opt$height ) ) { opt$height  = 8      }

data <- read.table(opt$file, header = T, check.names = F, sep = "\t", quote = "")
colname <- colnames(data)

if(dir.exists(opt$outdir) == FALSE) {
	dir.create(opt$outdir, recursive = T)
}

data.p <- NULL
legend = "none"
if(opt$group != "none" && is.numeric(opt$group) == T) {
	data.p <- data[, c(opt$x, opt$y, opt$z, opt$sample, opt$group)]
	colnames(data.p)[4] <- "sample"
	colnames(data.p)[5] <- "group"
	legend = "right"
} else {
	data.p <- data[, c(opt$x, opt$y, opt$z)]
	data.p$sample <- data[, opt$sample]
	data.p$group <- data.p$sample
	legend = "none"
}
samples <- unique(data.p$sample)
groups  <- unique(data.p$group)
data.p$group <- factor(data.p$group, levels = unique(data.p$group), order = T)

MaxMin <- function(x, ratio = 1, ref_min_max = "none") {
    values = unlist(strsplit(ref_min_max, ","))
    if(length(values) == 3) ratio = values[3]
    max_v = ifelse(max(x) > 0 , max(x) * ratio, max(x) / ratio)
    min_v = ifelse(min(x) > 0 , min(x) / ratio, min(x) * ratio)
    if(ref_min_max != "none") {
        if(values[1]!="none") min_v = as.numeric(values[1])
        if(values[2]!="none") max_v = as.numeric(values[2])
    }
    c(min_v,max_v)
}
xlim <- MaxMin(data.p[1], ratio = 1.1)
ylim <- MaxMin(data.p[2], ratio = 1.1)
zlim <- MaxMin(data.p[3], ratio = 1.1)

if ( is.null(opt$xlab   ) ) { opt$xlab    = colnames(data.p)[1] }
if ( is.null(opt$ylab   ) ) { opt$ylab    = colnames(data.p)[2] }
if ( is.null(opt$zlab   ) ) { opt$zlab    = colnames(data.p)[3] }

colors <- NULL
if(opt$color != "none") {
	colors <- unlist(strsplit(opt$color, ","))
} else {
	colors <- c("blue", "red", "green", "cyan", "yellow", "mediumpurple", "orange", "purple", 
	"pink", "gray", "wheat", "brown", "darkgreen", "greenyellow", "black", "chocolate")
}
if(length(colors) > length(unique(data.p$group))) {
	colors <- colors[1:length(unique(data.p$group))]
} else {
	colors <- colorRampPalette(colors)(length(unique(data.p$group)))
}
colors <- colors[as.numeric(data.p$group)]

shapes <- NULL
if(length(unique(data.p$group)) < 6) {
	shapes <- c(16, 17, 18, 15, 5)
	shapes <- shapes[1:length(unique(data.p$group))]
	shapes <- shapes[as.numeric(data.p$group)]
} else {
	shapes <- 16
}

group.num <- length(unique(data.p$group))
col.num <- ceiling(group.num / 30)
opt$width <- 10 + col.num
mar.right <- 5 + col.num * 8
inset <- -0.25 * col.num + 0.05
angle = 55
if(legend == "none") {
	opt$width = 8
	mar.right = 5.01
}

outfig <- paste0(opt$outdir, "/", opt$outname, ".pdf")
pdf(outfig, width = opt$width, height = opt$height)
par(lwd = 1, cex.main = 1.8, cex.lab = 1.5)
if(opt$isgrid == F) {
	s3d <- scatterplot3d(data.p[, 1:3], pch = "", grid = T, box = T, 
	main = opt$title, xlab = opt$xlab, ylab = opt$ylab, zlab = opt$zlab, 
	xlim = xlim, ylim = ylim, zlim = zlim, angle = angle, 
	mar = c(5.01, 4.01, 4.01, mar.right))
	s3d$points3d(data.p[, 1:3], pch = shapes, type="h", col = colors, cex = opt$pcex)
} else {
	s3d <- scatterplot3d(data.p[, 1:3], pch = "", grid = F, box = T, 
	main = opt$title, xlab = opt$xlab, ylab = opt$ylab, zlab = opt$zlab,
	xlim = xlim, ylim = ylim, zlim = zlim, angle = angle, col.axis = "grey50", 
	mar = c(5.01, 4.01, 4.01, mar.right))
	addgrids3d(data.p[, 1:3], grid = c("xy", "xz", "yz"), col.grid = "grey50",
	xlim = xlim, ylim = ylim, zlim = zlim, angle = angle)
	s3d$points3d(data.p[, 1:3], pch = shapes, type="h", col = colors, cex = opt$pcex) 
}
if(opt$islabel == T) {
	text(s3d$xyz.convert(data.p[, 1:3]), labels = samples, cex = opt$tcex, 
	pos = 3, offset = opt$vjust)
}
if(legend == "right") {
	legend('right', pch = unique(shapes), yjust = 0, inset = inset, legend = groups,
	text.font = 1, xpd = TRUE, cex = 1.2, bty = "n", xjust = 0.5, horiz = F,
	col = unique(colors))
}

dev.off()

