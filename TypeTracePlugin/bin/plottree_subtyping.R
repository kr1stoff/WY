#!/use/bin/env Rscript
# 炭疽的绘图计算距离
library('ape')
options(stringsAsFactors = F)

args = commandArgs(trailingOnly = TRUE)
f_in <- args[1]
d_out <- args[2]
query_strain <- args[3]


phy <- ape::read.tree(f_in)
pairwise_dist <- cophenetic.phylo(phy)
pdf(paste(d_out, "phylogenetic_tree.pdf", sep = '/'), width = 10)
plot(phy)
dev.off()

sample.cluster.tb <- as.data.frame(phy$edge)
sample.cluster <-
  sample.cluster.tb$V1[which(sample.cluster.tb$V2 == which(phy$tip.label ==
                                                             query_strain))]
node.no.tempt <- sample.cluster
sample.cluster.tb.tempt <-
  sample.cluster.tb[(sample.cluster.tb$V1 %in% node.no.tempt),]
cluster.samples <- c(node.no.tempt)
while (nrow(sample.cluster.tb.tempt) != 0) {
  sample.cluster.tb.tempt <-
    sample.cluster.tb[(sample.cluster.tb$V1 %in% node.no.tempt),]
  cluster.samples <- c(cluster.samples, sample.cluster.tb.tempt$V2)
  node.no.tempt <- sample.cluster.tb.tempt$V2
}
cluster.samples <-
  cluster.samples[(cluster.samples %in% (1:length(phy$tip.label)))]
cluster.samples.names <- phy$tip.label[cluster.samples]
cluster.samples.names <-
  cluster.samples.names[!(cluster.samples.names == query_strain)]

output.tb <- as.data.frame(array('NA', dim = c(1, 2)))
colnames(output.tb) <- c('Sample_genome', 'Subtype')
output.tb[1, 1] <- query_strain
if (length(cluster.samples.names) == 1) {
  output.tb[1, 2] <- cluster.samples.names
}
write.table(
  output.tb,
  paste(d_out, "result.tsv", sep = '/'),
  quote = F,
  col.names = T,
  row.names = F,
  sep = "\t"
)
