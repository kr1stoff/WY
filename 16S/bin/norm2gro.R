#!/usr/bin/Rscript
args <- commandArgs(T)
if (length(args) != 3) {
  print("Rscript norm2gro.R <table> <groups> <output>")
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
    colClass = "character",
    check = F,
    quote = '',
    comment = ''
  )
print(str(group))

#system.time({
func <- function(x) {
  if (length(x) == 1) {
    data[, as.character(x)]
  } else {
    rowMeans(data[, as.character(x)])
  }
}

gro <- by(group[, 1], group[, 2], func)
gro_tb <-
  data[, !colnames(data) %in% as.character(group[, 1]), drop = F]
for (i in unique(group[, 2])) {
  gro_tb <- data.frame(gro_tb, gro[i], check.names = F)
  colnames(gro_tb)[ncol(gro_tb)] <- i
}

write.table(
  gro_tb,
  file = args[3],
  sep = "\t",
  quote = F,
  row = F
)
