# ANI结果处理脚本
options(stringsAsFactors = F)
args <- commandArgs(trailingOnly = TRUE)
f_ani <- args[1]
query_strain <- args[2]
f_out <- args[3]

ani.res <- read.delim(f_ani, row.names = 1, header = T)
strain.ani.res <- ani.res[, query_strain]
names(strain.ani.res) <- row.names(ani.res)
strain.ani.res <-
  strain.ani.res[!(names(strain.ani.res) == query_strain)]
output.tb <- as.data.frame(array('uncertain', dim = c(1, 2)))
colnames(output.tb) <- c('Sample', 'Subtype')
output.tb[1, 1] <- query_strain
if (length(which(strain.ani.res > 0.95)) == 0) {
  output.tb[1, 2] <- 'NA'
} else{
  output.tb[1, 2] <-
    names(strain.ani.res)[order(strain.ani.res, decreasing = T)[1]]
}
write.table(
  output.tb,
  f_out,
  quote = F,
  col.names = T,
  row.names = F,
  sep = "\t"
)
