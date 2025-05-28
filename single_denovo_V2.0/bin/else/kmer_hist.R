#! /usr/local/bin/Rscript --vanilla
#jellyfish histo mer_counts.jf > kmer_hist.tsv

#args[1]  input the
#args[2]  output the png
args <- commandArgs(T)

x <- read.table(args[1], header = F, sep = " ", stringsAsFactors = F)
png(args[2])
plot(x, type = "l",xlab ="kmer_frequency",ylab = "kmer_total")
dev.off()
