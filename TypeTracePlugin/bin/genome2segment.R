#!/usr/bin/env Rscript
# 将基因组按照指定步长切片
options(warn =-1)
suppressMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)
fasta.input <- args[1]
output.file <- args[2]
window.size <- as.numeric(args[3])
increment <- as.numeric(args[4])
myFastaFile <- readDNAStringSet(fasta.input)
lcm.incre.window <- (max(increment, window.size) - 1)
repeat {
  lcm.incre.window <- lcm.incre.window + 1
  if (lcm.incre.window %% increment == 0 &&
    lcm.incre.window %% window.size == 0)
    break
}

increment.multiple <- ((lcm.incre.window / increment) - 1)
for (i in 0:increment.multiple) {
  myFastaFile.increment <-
    subseq(myFastaFile,
           start = (increment * i + 1),
           end = length(myFastaFile[[1]]))
  increment.id <- paste0('increment_', increment * i)
  segment.no <- length(myFastaFile.increment[[1]]) %/% window.size
  for (j in 1:segment.no) {
    myFastaFile.segment <-
      subseq(
        myFastaFile.increment,
        start = (window.size * (j - 1) + 1),
        end = (j * window.size)
      )
    names(myFastaFile.segment) <-
      paste0(increment.id,
             '_pos_start_',
             (increment * i + (window.size * (j - 1) + 1)),
             '_pos_end_',
             (increment * i + (j * window.size)))
    if (i == 0 & j == 1) {
      segment.merge <- myFastaFile.segment
    } else {
      segment.merge <- c(segment.merge, myFastaFile.segment)
    }
  }
}
writeXStringSet(segment.merge, filepath = output.file)
