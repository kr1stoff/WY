#!/Bio/bin/Rscript-3.6.0
library(RColorBrewer)

args  <- commandArgs(TRUE)
color <- args[1]
count <- args[2]

colors <- unlist(strsplit(color, split = ","))
colors_pal <- colorRampPalette(colors)(count)

cat(paste(unlist(colors_pal), collapse=","))

