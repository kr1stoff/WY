#!/usr/bin/Rscript
args <- commandArgs()
bin <- dirname(sub('--file=', '', args[grep('--file=', args)]))

library(optparse)
library(data.table)

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "input otu file [default %default]"
  ),
  make_option(
    c("-t", "--tree"),
    type = "character",
    default = NULL,
    help = "input phylo file (for Phylogenetic Diversity) [default %default]"
  ),
  make_option(
    c("-d", "--depth"),
    type = "character",
    default = 'min',
    help = "depth for rarefaction min/max [default %default]"
  ),
  make_option(
    c("-r", "--rarefaction"),
    action = "store_true",
    default = F,
    help = "whether to calculate rarefaction curve [default %default]"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = '.',
    help = "output dir [default %default]"
  )
)
opts <-
  parse_args(
    OptionParser(usage = "%prog [options]", option_list = option_list),
    positional_arguments = F
  )

if (is.null(opts$input)) {
  cat("Rscript rarefaction -i all.otus.exp_for_unifrac [-r]\n")
  q()
}
if (!dir.exists(opts$outdir))
  dir.create(opts$output, recursive = T)
data <-
  read.table(
    opts$input,
    sep = "\t",
    header = T,
    row = 1,
    check.names = F
  )
otu <- t(as.matrix(data))

rarefyR <- function(x, depth) {
  if (sum(x) > depth) {
    y <- sample(rep(1:length(x), x), depth)
    y.tab <- table(y)
    z <- numeric(length(x))
    z[as.numeric(names(y.tab))] <- y.tab
    z
  } else {
    rep(NA, length(x))
  }
}

aceR <- function(y, abu = 10) {
  if (any(is.na(y))) {
    NA
  } else {
    y <- y[y != 0]
    n1 <- sum(y == 1)
    rare <- numeric(abu)
    y.tab <- table(y)
    rare <- y.tab[as.character(1:10)]
    rare[is.na(rare)] <- 0
    n_rare <- sum(y[y <= abu])
    s_rare <- sum(rare)
    c_ace <- 1 - n1 / n_rare
    s_abund <- sum(y > abu)
    gama2ace <-
      max(s_rare / c_ace * sum(1:abu * 0:(abu - 1) * rare) /
            n_rare /
            (n_rare -
              1) - 1, 0)
    s_abund + s_rare / c_ace + n1 / c_ace * gama2ace
  }
}

callace <- aceR

diversity2 <- function(x) {
  t <- rowSums(x)       # tags num
  p <- x / t            # otu percentage
  n1 <- rowSums(x == 1) # singletons
  n2 <- rowSums(x == 2) # doubleteons
  sobs <- rowSums(p != 0)
  chao <- sobs + n1 * (n1 - 1) / 2 / (n2 + 1) # change chao1 to chao
  ace <- apply(x, 1, callace)
  shannon <- rowSums(-p * log(p, 2), na.rm = T)
  shannon[shannon == 0] <- NA
  simpson <- 1 - rowSums(p * p)
  goods_coverage <- 1 - (n1 / t)
  data.frame(sobs,
             shannon,
             simpson,
             chao,
             ace,
             goods_coverage,
             check.names = F)
}

# summary alpha diversity
otu_diversity <- diversity2(otu)
otu_diversity <-
  data.frame(index = row.names(otu_diversity), otu_diversity)
otu_diversity$pielou <-
  otu_diversity$shannon / log(otu_diversity$sobs, 2)

runPd <- F
if (!is.null(opts$tree)) {
  tree <- ape::read.tree(opts$tree)
  source(paste0(bin, '/pd.R'))
  pd_whole_tree <- pd(otu, tree, include.root = F)
  otu_diversity$pd <- pd_whole_tree$PD
  runPd <- T
}

diversity_out <-
  paste0(opts$outdir, '/', 'summary.alpha_diversity.xls')
fwrite(otu_diversity,
       file = diversity_out,
       sep = "\t",
       quote = F)
if (!opts$rarefaction)
  q()

# load cpp
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
sourceCpp(paste0(bin, '/rarefyC.cpp'))

rftb <- function(tb, depth) {
  tbr <- t(apply(tb, 1, rarefyC, depth = depth))
  colnames(tbr) <- colnames(tb)
  return(tbr)
}

# diversity for rarefaction
index <-
  c('sobs', 'shannon', 'simpson', 'chao', 'ace', 'goods_coverage')

## sample interval and sample frequency
tags_num <- rowSums(otu)
max_depth <- max(tags_num)
min_depth <- min(tags_num)
n <- 400
step <- 400
sp <- c()
if (opts$depth != 'min')
  min_depth <- max_depth
while (n < min_depth) {
  sp <- c(sp, n)
  if (n >= 10000)
    step <- 500
  if (n >= 20000)
    step <- 1000
  if (n >= 30000)
    step <- 2500
  if (n >= 40000)
    step <- 5000
  n <- n + step
}
freq <- 10

## Prepare
mt <-
  matrix(0, ncol = (nrow(otu) + 1), nrow = (length(sp) + 1) * freq)
colnames(mt) <- c('numsampled', row.names(otu))

rarefaction <- list()
for (i in index)
  rarefaction[[i]] <- mt

## Run Alpha Rarefaction
system.time({
  ## Default Rarefaction
  num <- freq
  for (i in sp) {
    print(paste(Sys.time(), i))
    rff_list <- lapply(rep(i, freq), rftb, tb = otu) # rep 400 10
    alpha_list <- lapply(rff_list, diversity2)

    for (y in index) {
      for (j in 1:freq)
        rarefaction[[y]][num + j,] <- c(i, alpha_list[[j]][, y])
    }
    num <- num + freq
    #print(str(rarefaction));q()
  }

  ## PD Rarefaction
  pd2 <- function(tab, tree, include.root) {
    isna <- apply(tab, 1, function(x)
      any(is.na(x)))
    tab2 <- tab[!isna, , drop = F]
    pd_out <- pd(tab2, tree, include.root)
    rownames(pd_out) <- rownames(tab2)
    pd_out <- pd_out[rownames(tab),]
  }

  if (runPd) {
    sp <-
      c(
        500,
        1000,
        2000,
        4000,
        8000,
        15000,
        25000,
        35000,
        45000,
        55000,
        65000,
        75000,
        85000,
        95000,
        105000,
        115000
      )
    sp <- sp[sp <= min_depth]
    rarefaction[['pd']] <-
      matrix(0,
             ncol = (nrow(otu) + 1),
             nrow = (length(sp) + 1) * freq)
    colnames(rarefaction[['pd']]) <- c('numsampled', row.names(otu))
    num <- freq
    for (i in sp) {
      rff_list <- lapply(rep(i, freq), rftb, tb = otu) # rep 400 10
      print(paste(Sys.time(), 'pd:', i))
      pd_list <-
        lapply(rff_list, pd2, tree = tree, include.root = FALSE)
      for (j in 1:freq)
        rarefaction[['pd']][num + j,] <- c(i, pd_list[[j]][, 'PD'])
      num <- num + freq
    }
  }

})

for (i in names(rarefaction)) {
  f <- paste0(opts$outdir, '/', i, '.rarefaction.tsv')
  dt <- as.data.table(rarefaction[[i]])
  fwrite(dt,
         file = f,
         sep = "\t",
         quote = F)
}
