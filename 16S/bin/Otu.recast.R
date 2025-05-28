#!/usr/bin/Rscript
options(warn = -1)
library(funr)
library(stringr)
library(yaml)
BIN = dirname(sys.script())
f_software <- str_c(BIN, "../config/software.yml", sep = '/')
software <- yaml.load_file(f_software)
# .libPaths(software$Rlibpath)
library(optparse)
library(magrittr)
library(data.table)
library(hash)
library(logging)
basicConfig()

#Load packages needed
option_list <- list(
  make_option(c("-j", "--json"),
              type = "character",
              default = NULL,
              help = "JSON file [default %default]"),
  make_option(c("-i", "--otu"),
              type = "character",
              default = NULL,
              help = "OTU table [default %default]"),
  make_option(c("-a", "--absolute"),
              type = "integer",
              default = NULL,
              help = "Filter tags num [default %default]"),
  make_option(c("-t", "--top"),
              type = "integer",
              default = NULL,
              help = "Pick top abundance OTU [default %default]"),
  make_option(c("-r", "--relative"),
              type = "numeric",
              default = NULL,
              help = "Filter relative abundance [default %default]"),
  make_option(c("-s", "--samples"),
              type = "character",
              default = NULL,
              help = "Samples list seperate by comma [default %default]"),
  make_option(c("-d", "--depth"),
              type = "integer",
              default = NULL,
              help = "OTU rarefarefy depth [default %default]"),
  make_option(c("-z", '--threshold'),
              type = "integer",
              default = 0,
              help = "Tags num threshold [default %default]"),
  make_option(c("-l", '--limit'),
              type = "integer",
              default = NULL,
              help = "Number of samples' pass threshold will remain [default %default]"),
  make_option(c("-o", "--output"),
              type = "character",
              default = ".",
              help = "Output directory [default %default]")
)
opts <- parse_args(OptionParser(usage = "%prog [options]",
                                option_list = option_list),
                   positional_arguments = F)

#### My Functions
rarefy <- function(x, depth) {
  # Rarefy OTU
  if (depth < sum(x)) {
    y <- sample(rep(1:length(x), x), depth)
    y.tab <- table(y)
    z <- numeric(length(x))
    z[as.numeric(names(y.tab))] <- y.tab
    return(z)
  } else {
    return(x)
  }
}

splitTax <- function(x) {
  x <- gsub(' ', '_', x)
  x <- gsub(';.__', ';', x)
  arr <- str_split(x, pattern = ';')[[1]]
  if (arr[1] == 'Root') {
    arr <- arr[-1]
  } else {
    arr[1] <- gsub('k__', '', arr[1])
  }
  arr[ifelse(is.na(arr[1:7]) | arr[1:7] == "", T, F)] <- 'Unclassified'
  return(arr)
}


#### Main
if (!dir.exists(opts$output)) dir.create(opts$output)
opts$output <- normalizePath(opts$output)

## read json
if (!is.null(opts$json)) {
  json <- rjson::fromJSON(file = opts$json)
  for (i in names(json)) if (is.null(opts[[i]])) opts[[i]] = json[[i]]
} else {
  json <- list(otu = NULL)
}

## read otu table with header and rownames
loginfo("Read OTU table")
if (is.null(opts$otu)) {
  stop("OTU table was required!")
} else if (is.null(json$otu)) {
  opts$otu <- fread(opts$otu, sep = "\t", check = F, header = T)
} else {
  opts$otu <- fread(opts$otu, sep = "\t", check = F, header = T)
}

## split table to tax & otu
colnum <- ncol(opts$otu)
rownum <- nrow(opts$otu)
colnames(opts$otu)[c(1, colnum)] <- c('otuid', 'taxonomy')

opts$rn <- opts$otu[, otuid]                  # vector
opts$tax <- opts$otu[, taxonomy]               # vector
opts$otu[, c('otuid', 'taxonomy') := NULL]     # data.table
loginfo("numOtus : %s", rownum)
used <- rep(T, rownum)

## choose samples to analysis
loginfo("Filter samples ...")
if (!is.null(opts$samples)) {
  opts$samples <- str_split(opts$samples, ",")[[1]]
  if (!is.null(opts$old_samples)) {
    opts$old_samples <- str_split(opts$old_samples, ",")[[1]]
    if (length(opts$samples) != length(opts$old_samples)) stop("Samples num no matched!")
    opts$otu <- opts$otu[, opts$old_samples, with = F]
    colnames(opts$otu) <- opts$samples
  } else {
    opts$otu <- opts$otu[, opts$samples, with = F]
  }
} else {
  opts$samples <- colnames(opts$otu)
}
loginfo("Samples : %s", opts$samples)

## read species
loginfo("Read tax")
tax.level <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

opts$tax.new <- lapply(opts$tax, splitTax) %>%
  do.call(what = rbind) %>%
  as.data.table()
setnames(opts$tax.new, 1:7, tax.level)

keys <- lapply(2:7, function(x) paste(tax.level[x], opts$tax.new[, x:x][[1]], sep = ':')) %>% unlist()
values <- lapply(2:7, function(x) apply(opts$tax.new[, 1:x], 1, paste, collapse = ';')) %>% unlist()
opts$short2long <- hash(keys, values)

loginfo("Read tax done!")

# Step : Filer Tax --> "retain":["Domain,k__Bacteria","Phylum,p__Proteobacteria"]
if (!is.null(opts$species) &&
  opts$species != "" &&
  !grepl(".*,$", opts$species[[1]][1])) {
  loginfo("[Filter] Choose species %s", opts$species)
  taxon <- c(Domain = "k__", Phylum = "p__", Class = "c__", Order = "o__",
             Family = "f__", Genus = "g__", Species = "s__")
  type <- names(opts$species[1])
  loginfo(type, opts$species[[1]])

  hit <- !used
  for (i in opts$species[[1]]) {
    arr <- str_split(i, pattern = ",")[[1]]
    level <- arr[1]
    tax <- sub(".__", "", arr[2])
    hit <- hit | opts$tax.new[[level]] == tax
    loginfo("numHits : %s", sum(hit))
  }
  if (type == "retain" || type == "remain") {
    used <- used & hit
  } else {
    used <- used & !hit
  }

  loginfo("numOtus : %s", sum(used))
  opts$rn <- opts$rn[used]
  opts$tax <- opts$tax[used]
  opts$otu <- opts$otu[used, , drop = F]
  opts$tax.new <- opts$tax.new[used, , drop = F]
}

## Step : Remove OTU
if (!is.null(opts$oturm)) {
  loginfo("[Filter] Remove OTU : %s", opts$oturm)
  oturm <- str_split(opts$oturm, ',')[[1]]
  used <- !(opts$rn %in% oturm)
  opts$rn <- opts$rn[used]
  opts$otu <- opts$otu[used,]
  opts$tax <- opts$tax[used]
  opts$tax.new <- opts$tax.new[used,]
}

if (!is.null(opts$depth) && opts$depth != 0) {
  opts$depth <- ifelse((is.na(opts$depth) | opts$depth == "min"),
                       min(colSums(opts$otu)), as.integer(opts$depth))
  loginfo("[Rarefy] depth : %s", opts$depth)
  opts$otu <- apply(opts$otu, 2, rarefy, depth = opts$depth)
}

# Step : Filter Tag Num
numOtus <- nrow(opts$otu)
loginfo("numOtus : %s", numOtus)
opts$absolute <- ifelse(is.null(opts$absolute), 0, as.numeric(opts$absolute))
loginfo("[Filter] Filter tag num ( equal and less than %s )", opts$absolute)
rowsum <- rowSums(opts$otu)
used <- rowsum > opts$absolute
if (!is.null(opts$limit)) {
  loginfo("[Filter] %s samples' tags >= %s will remain ...", opts$limit, opts$threshold)
  passnum <- apply(opts$otu, 1, function(x) sum(x >= opts$threshold))
  used <- used & (passnum >= opts$limit)
}
opts$rn <- opts$rn[used]
opts$tax <- opts$tax[used]
opts$otu <- opts$otu[used,]
opts$tax.new <- opts$tax.new[used,]

numOtus <- nrow(opts$otu)
loginfo("numOtus : %s", numOtus)

# choose top
rowsum <- rowSums(opts$otu)
order <- order(rowsum, decreasing = T)
opts$rn <- opts$rn[order]
opts$tax <- opts$tax[order]
opts$otu <- opts$otu[order,]
opts$tax.new <- opts$tax.new[order,]

if (!is.null(opts$top) &&
  opts$top != 0 &&
  numOtus > as.numeric(opts$top)) {
  loginfo("[Filter] Get top %s", opts$top)
  opts$top <- as.numeric(opts$top) # convert to numeric

  # filter 
  opts$rn <- opts$rn[1:opts$top]
  opts$tax <- opts$tax[1:opts$top]
  opts$otu <- opts$otu[1:opts$top,]
  opts$tax.new <- opts$tax.new[1:opts$top,]

  # stat
  numOtus <- nrow(opts$otu)
  loginfo("numOtus : %s", nrow(opts$otu))
}

## output table
# calculate abundance
loginfo("Calculate abundance")
colsum <- colSums(opts$otu)
opts$abundance <- round(t(t(opts$otu) / colsum * 100), 6) # matrix
opts$abundance[is.nan(opts$abundance)] <- 0

if (is.null(opts$name)) opts$name <- "all.otus"

## output tab
loginfo('Output tab')
abu_name <- paste0(opts$output, "/all.otus.abu.xls")
tag_name <- paste0(opts$output, "/all.otus.exp_for_unifrac")
tab_name <- paste0(opts$output, "/all.otus.tab")
fwrite(data.table(OTU = opts$rn, opts$abundance),
       file = abu_name, sep = "\t", quote = F)
fwrite(data.table(OTU = opts$rn, opts$otu),
       file = tag_name, sep = "\t", quote = F)
fwrite(data.table(OTU = opts$rn, opts$otu, Taxonomy = opts$tax),
       file = tab_name, sep = "\t", quote = F)

## output shared
loginfo('Output shared')
shared <- data.table(label = 0.03, Group = opts$samples, numOtus = nrow(opts$otu),
                     apply(opts$otu, 2, paste0, collapse = "\t"))
setnames(shared, 4, paste0(opts$rn, collapse = "\t"))
fwrite(shared, file = paste0(opts$output, "/all.otus.shared"), sep = "\t", quote = F)

## output profile
loginfo('Output profile')
profile.name <- paste0(opts$output, "/", "all.otus.profiling.xls")
tag.head <- paste0(opts$samples, "_tags")
abu.head <- paste0(opts$samples, "_relative_abundance")
tax.tmp <- as.matrix(opts$tax.new)
tax.tmp[tax.tmp == 'Unclassified'] <- ''

profile.out <- data.table(opts$rn, rowSums(opts$otu), opts$otu, opts$abundance, tax.tmp)
setnames(profile.out, c('OTU', 'Total_tags', tag.head, abu.head, tax.level))
fwrite(profile.out, file = profile.name, sep = "\t", quote = F)

## output min0.1
loginfo('OTU min0.1')
name <- paste0(opts$output, '/profiling.OTU.min0.1.xls')
is_min <- apply(opts$abundance, 1, function(x) max(x) > 0.1)
otu_min <- data.table(OTU = opts$rn, opts$abundance)[is_min,]
fwrite(otu_min, file = name, sep = "\t", quote = F)

## output taxa
loginfo('Stat taxa abundance')
all_tax.df <- data.frame()
tax_stat <- matrix(nrow = ncol(opts$otu), ncol = 7)
colnames(tax_stat) <- tax.level
rownames(tax_stat) <- colnames(opts$otu)

for (i in tax.level) {
  # stat tag
  tax <- opts$tax.new[[i]]
  tag.mt <- apply(opts$otu, 2, function(x) tapply(x, tax, sum)) # mt is matrix
  abu.mt <- apply(opts$abundance, 2, function(x) tapply(x, tax, sum))
  abu.mt <- round(abu.mt, 4)

  if (is.vector(tag.mt)) {
    tag.mt <- matrix(tag.mt, nrow = 1, dimnames = list(unique(tax), names(tag.mt)))
    abu.mt <- matrix(abu.mt, nrow = 1, dimnames = list(unique(tax), names(abu.mt)))
  }

  # order
  tag.df <- data.frame(tax = rownames(tag.mt), tag.mt, check.names = F, stringsAsFactors = F)
  abu.df <- data.frame(tax = rownames(abu.mt), abu.mt, check.names = F, stringsAsFactors = F)
  ord <- order(rowSums(abu.mt), decreasing = T) # order by abundance
  tag.df <- tag.df[ord, , drop = F]
  abu.df <- abu.df[ord, , drop = F]

  # get unclassified
  isUnclass <- tag.df$tax == 'Unclassified'
  tag.unclass <- tag.df[isUnclass, , drop = F]
  abu.unclass <- abu.df[isUnclass, , drop = F]
  tag.df <- tag.df[!isUnclass, , drop = F]
  abu.df <- abu.df[!isUnclass, , drop = F]

  # stat unclassified
  tax_stat[, i] <- colSums(tag.df[-1])

  # tags
  tag.df <- rbind(tag.df, tag.unclass)
  colnames(tag.df)[1] <- i
  name <- paste0(opts$output, "/", 'profiling.', i, '.tag.xls')
  fwrite(tag.df, file = name, sep = "\t", quote = F)

  # tag + abu
  tmp.df <- rbind(abu.df, abu.unclass)[-1]
  all.df <- data.frame(tag.df, tmp.df)
  colnames(all.df)[-1] <- c(tag.head, abu.head)
  opts$short2long[[paste(i, 'Unclassified', sep = ':')]] <- 'Unclassified'
  if (i != 'Domain') {
    all.df[, 1] <- lapply(paste(i, row.names(all.df), sep = ':'), function(x) opts$short2long[[x]]) %>% unlist()
  }
  name <- paste0(opts$output, "/", 'profiling.', i, '.xls')
  fwrite(all.df, file = name, sep = "\t", quote = F)

  # all_tax data.frame
  colnames(all.df)[1] <- 'taxonomy'
  all_tax.df <- rbind(all_tax.df, all.df[row.names(all.df) != 'Unclassified', , drop = F])

  # min0.1
  out <- abu.df[apply(abu.df[, -1, drop = F], 1, max) > 0.1, , drop = F]
  colnames(out)[1] <- i
  name <- paste0(opts$output, "/", 'profiling.', i, '.min0.1.xls')
  fwrite(out, file = name, sep = "\t", quote = F)

  # top10
  Other <- c()
  if (nrow(abu.df) > 10) {
    Other <- colSums(abu.df[11:nrow(abu.df), -1, drop = F])
    Other <- c('Other', Other)
    abu.df <- abu.df[1:10,]
  }
  abu.df <- rbind(abu.df, Other, abu.unclass)
  colnames(abu.df)[1] <- i
  name <- paste0(opts$output, "/", 'profiling.', i, '.top10.xls')
  fwrite(abu.df, file = name, sep = "\t", quote = F)

}
name <- paste0(opts$output, '/profiling.all.taxonomy.xls')
fwrite(all_tax.df, file = name, sep = "\t", quote = F)

tax_stat_file <- paste0(opts$output, "/all.taxonomy.stat.xls")
write.table(tax_stat, file = tax_stat_file, sep = "\t", quote = F, col.names = NA)
tax_stat_file2 <- paste0(opts$output, "/all.taxonomy.stat.forFig.xls")
tax_stat_only <- tax_stat[, 1:6] - tax_stat[, 2:7]

tax_stat2 <- data.frame(
  Level = rep(colnames(tax_stat), rep(nrow(tax_stat), 7)),
  Sample = rownames(tax_stat),
  Sum = as.vector(tax_stat),
  Number = c(as.vector(tax_stat_only), tax_stat[, 7])
)
fwrite(tax_stat2, file = tax_stat_file2, sep = "\t", quote = F)

loginfo('All done!')
