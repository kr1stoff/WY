#!/usr/bin/Rscript
args <- commandArgs(T)
if (length(args) != 3) {
  print("Rscript summary_alpha.R <table> <groups> <outdir>")
  q()
}

data <-
  read.table(
    args[1],
    sep = "\t",
    header = T,
    stringsAsFactors = F,
    check = F,
    quote = '',
    comment = ''
  )
group <-
  read.table(
    args[2],
    sep = "\t",
    header = T,
    stringsAsFactors = F,
    check = F,
    quote = '',
    comment = ''
  )

# filter summary
used <- !grepl("_lci|_hci|label|group|index", colnames(data))
dat <- data[, used]
row.names(dat) <- data[, 1]
# split group

func <- function(x) {
  x <- as.character(x)
  if (length(x) == 1) {
    dat[x,]
  } else {
    colMeans(dat[x,])
  }
}

gsd <- function(x) {
  x <- as.character(x)
  if (length(x) == 1) {
    dat[x,]
  } else {
    apply(dat[x,], 2, sd)
  }
}

gro <- by(group[, 1], group[, 2], func)
gro_tb <- do.call(rbind, gro)
gro_tb <- data.frame(group = names(gro), gro_tb, check.names = F)
write.table(
  gro_tb,
  file = paste0(args[3], '/summary.alpha_diversity.gro.xls'),
  sep = "\t",
  quote = F,
  row = F
)

gro_sd <- by(group[, 1], group[, 2], gsd)
gsd_tb <- do.call(rbind, gro_sd)
gsd_tb <- data.frame(group = names(gro), gsd_tb, check.names = F)
write.table(
  gsd_tb,
  file = paste0(args[3], '/summary.alpha_diversity.sd.xls'),
  sep = "\t",
  quote = F,
  row = F
)
