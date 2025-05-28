options(warn = -1,
        stringsAsFactors = F)
library(stringr)


args <- commandArgs(trailingOnly = TRUE)
blast.res <- read.delim(args[1], header = F)
blast.threshold <- args[2]
output.tb <- read.delim(args[3], header = T)

# Define marker gene
EPEC_eae <- c('eae')
EPEC_bfp <- paste0('bfp', LETTERS)
ETEC_elt <- c('eltA', 'eltB','AAA24093','AAA24094')
ETEC_est <- c('estIa','AAA23990')
ETEC_sta <- paste0('sta', LETTERS)
EAEC_aggR <- c('aggR')
EIEC_inv <- c('invA', 'invE')
EIEC_ipaH <- c('ipaH')
EHEC_stx1 <- c('stx1vA', 'stx1vB')
EHEC_stx2 <- c('stx2A',
               'stx2B',
               'stx2d1A',
               'stx2eA',
               'stx2eB',
               'stx2fA',
               'stx2fB')

# Main
blast.res <- blast.res[(blast.res$V3 > blast.threshold), ]
output.tb$Pathotype <- 'uncertain'
gene.names <-
  unique(str_remove(str_remove(
    str_split(blast.res$V2, " ", simplify = T)[, 2], "\\("
  ), "\\)"))

subtype.names <- c()
if (length(which(gene.names %in% EPEC_eae)) > 0 |
    length(which(gene.names %in% EPEC_bfp)) > 0) {
  subtype.names <- c(subtype.names, 'EPEC')
}
if (length(which(gene.names %in% ETEC_elt)) > 0 |
    length(which(gene.names %in% ETEC_est)) > 0 |
    length(which(gene.names %in% ETEC_sta)) > 0) {
  subtype.names <- c(subtype.names, 'ETEC')
}
if (length(which(gene.names %in% EAEC_aggR)) > 0) {
  subtype.names <- c(subtype.names, 'EAEC')
}
if (length(which(gene.names %in% EIEC_inv)) > 0 |
    length(which(gene.names %in% EIEC_ipaH)) > 0) {
  subtype.names <- c(subtype.names, 'EIEC')
}
if (length(which(gene.names %in% EHEC_stx1)) > 0 |
    length(which(gene.names %in% EHEC_stx1)) > 0) {
  subtype.names <- c(subtype.names, 'EHEC')
}
subtype.names <- paste0(subtype.names, collapse = ',')
if (nchar(subtype.names) > 0) {
  output.tb$Pathotype <- subtype.names
}
write.table(
  output.tb,
  args[3],
  quote = F,
  col.names = T,
  row.names = F,
  sep = "\t"
)
