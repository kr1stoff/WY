#!/usr/bin/env Rscript
options(warn = -1)
options(bitmapType = "cairo")
suppressPackageStartupMessages({
  library(funr)
  library(stringr)
  library(yaml)
  BIN = dirname(sys.script())
  f_software <- str_c(BIN, "../config/software.yml", sep = '/')
  software <- yaml.load_file(f_software)
  library(vegan)
  library(ape)
  library(GUniFrac)
  library(logging)
})
basicConfig()


#### Main
args <- commandArgs(T)
if (length(args) != 3) {
  logerror("Rscript dist.r <file> <method|tree> <outdir>")
  logerror("method : bray, jaccard ")
  q()
}

data <-
  read.table(
    args[1],
    sep = "\t",
    head = T,
    row.names = 1,
    check.names = F
  )
if (nrow(data) < 2) {
  logerror("The taxa or otu count less than 2, can draw bray,but ugly , so quit")
  quit("no")
}
outdir <- args[3]
if (!dir.exists(outdir))
  dir.create(outdir, recursive = T)


# function plot tree
plot_upgma <- function(mt, outpre) {
  hc <- hclust(as.dist(mt), method = "average")
  hc_phylo <- as.phylo(hc)
  # output tree
  write.tree(hc_phylo, file = paste(outpre, ".nwk", sep = ""))
  # output matrix
  mt_out <- data.frame(lable = row.names(mt), mt, check.names = F)
  write.table(
    mt_out,
    file = paste0(outpre, ".matrix.txt"),
    quote = F,
    sep = "\t",
    row = F
  )

  # output upgma plot
  w <- ifelse(nrow(mt) > 40, nrow(mt) * 0.2, 8)
  png(
    paste0(outpre, ".png"),
    units = 'in',
    res = 300,
    width = w,
    height = 6
  )
  plot(hc, sub = "", xlab = "")
  invisible(dev.off())
  pdf(paste0(outpre, ".pdf"), width = w, height = 6)
  plot(hc, sub = "", xlab = "")
  invisible(dev.off())
}

# run dist
if (args[2] %in% c('bray', 'jaccard')) {
  print(args[2])
  dist <- vegdist(t(data), method = args[2])
  mt <- as.matrix(dist)
  row.names(mt) <- labels(dist)
  colnames(mt) <- labels(dist)
  outpre <- paste0(outdir, "/", args[2])
  plot_upgma(mt, outpre)
} else if (file.exists(args[2])) {
  tree <- read.tree(file = args[2])
  if (!is.rooted(tree)) {
    loginfo("unroot tree, multi2di ...")
    tree <- multi2di(tree)
  }

  label <- intersect(sort(tree$tip.label), row.names(data))
  if (nrow(data) > length(label)) {
    loginfo("filter otu ...")
    data <- data[label,]
    logdebug(dim(data))
  }
  set.seed(999)
  otu.tab.rff <- Rarefy(t(data))$otu.tab.rff
  unifracs <-
    GUniFrac(otu.tab.rff, tree, alpha = c(0, 0.5, 1))$unifracs
  # Weighted UniFrac
  wu <- unifracs[, , "d_1"]
  wu_pre <- paste0(outdir, "/weighted_unifrac")
  plot_upgma(wu, wu_pre)
  # Unweighted UniFrac
  uw <- unifracs[, , "d_UW"]
  if (sum(uw) == 0) {
    loginfo("rm unweighted unifrac ...")
    unlink(paste0(outdir, "/unweighted_unifrac*"))
  } else {
    uw_pre <- paste0(outdir, "/unweighted_unifrac")
    plot_upgma(uw, uw_pre)
  }

} else {
  logerror("%s was neither method nor tree!", args[2])
}
