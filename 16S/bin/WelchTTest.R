#!/Bio/bin/Rscript
# @Author: MingJia
# @Date:   2019-08-30 09:47:56
# @Last Modified by:   MingJia
# @Last Modified time: 2021-08-30 17:20:45
options(warn = -1)
library(optparse)

option_list <- list(
  make_option(c("--table"),
              type = "character",
              help = "The abundance talbe"),
  make_option(c("--group"),
              type = "character",
              help = "The group info file"),
  make_option(c("-c", "--compare"),
              type = "character",
              help = "The compare group info like \"A&B,B&C\""),
  make_option(c("--top"),
              type = "character",
              default = NULL,
              help = "The Top N to calculate"),
  make_option(c("--plot"),
              type = "character",
              help = "The script to plot error bar"),
  make_option(c("-o", "--out"),
              type = "character",
              help = "The out put dir")
)
opts <- parse_args(
  OptionParser(
    usage = "%prog [options]",
    option_list = option_list,
    description = "\nWelch T test for meta analysis"
  ),
  positional_arguments = F
)
df <- read.csv(opts$table,
               sep = '\t',
               check.names = F,
               row.names = 1,
               quote = "")
group <- read.table(opts$group,
                    check.names = F,
                    sep = "\t",
                    quote = "")
colnames(group) <- c("Samples", "Groups")
all_list <- by(group$Samples, group$Groups, function(x) df[, as.character(x)])
compares <- unlist(strsplit(opts$compare, split = ","))

for (compare in compares) {
  print(compare)
  union <- unlist(strsplit(compare, split = "&"))
  if (length(union) == 2) {
    ## filter data
    mt_1 <- all_list[[union[1]]]
    mt_2 <- all_list[[union[2]]]
    used <- (apply(mt_1, 1, sd) + apply(mt_2, 1, sd)) != 0
    mt_1 <- mt_1[used,]
    mt_2 <- mt_2[used,]
    if (!is.null(opts$top)) {
      select_length <- min(sum(used), as.numeric(opts$top))
      print(as.numeric(select_length))
      used <-
        sort(rowSums(mt_1) + rowSums(mt_2),
             decreasing = TRUE,
             index = T)$ix[1:select_length]
      mt_1 <- mt_1[used,]
      mt_2 <- mt_2[used,]
    }

    ## do t.test
    labels <- rownames(mt_1)
    labels_num <- length(labels)
    fold_name <- paste("fold(", union[2], "/", union[1], ")", sep = "")

    ttest.draw <- matrix(nrow = labels_num, ncol = 7)
    colnames(ttest.draw) <- c("labels", union[1], union[2], 'diff', "conf_min", "conf_max", "p_value")
    ttest.out <- matrix(nrow = labels_num, ncol = 5)
    colnames(ttest.out) <- c("labels", union[1], union[2], fold_name, "p-value")

    for (j in 1:labels_num) {
      a <- as.numeric(mt_1[j,])
      b <- as.numeric(mt_2[j,])
      ## welch's t-test and wilcox test
      ## is paired is T, samples must be paired one by one
      ttest_out <- t.test(a, b)
      ## statistics
      estimate <- ttest_out$estimate
      conf_range <- ttest_out$conf.int * -1
      diff <- estimate[2] - estimate[1]
      fold <- estimate[2] / estimate[1]
      ## output data
      ttest.draw[j,] <-
        c(labels[j], estimate, diff, conf_range, ttest_out$p.value)
      ttest.out[j,] <-
        c(labels[j], round(c(estimate, fold), 8), ttest_out$p.value)
    }
    ## output file name
    compare <- paste(union[1], "_vs_", union[2], sep = "")
    name.draw <-
      paste(opts$out, "/", compare, ".t-test.draw.xls", sep = "")
    name.pic <-
      paste(opts$out, "/", compare, '.t-test.extended_error_bar', sep = "")
    name.ttest <-
      paste(opts$out, "/", compare, ".t-test.xls", sep = "")

    ## t-test output file
    ttest.p_value <- as.numeric(ttest.out[, "p-value"])
    ttest.significant <- ifelse(ttest.p_value < 0.05, 'yes', 'no')
    ttest.out <-
      data.frame(ttest.out,
                 significant = ttest.significant,
                 check.names = F)
    write.table(
      ttest.out,
      file = name.ttest,
      sep = "\t",
      col.names = T,
      row.names = F,
      quote = F
    )
    ttest.draw <-
      data.frame(ttest.draw, check.names = F)[ttest.p_value < 0.05, , drop = F]
    if (nrow(ttest.draw) > 0) {
      write.table(
        ttest.draw,
        file = name.draw,
        sep = "\t",
        col.names = T,
        row.names = F,
        quote = F
      )
      cmd <-
        paste(
          "/Bio/User/renchaobo/software/miniconda3/envs/R3.6.1/bin/Rscript",
          opts$plot,
          name.draw,
          name.pic
        )
      print(cmd)
      system(cmd)
    } else {
      print("No Diff find")
    }
  } else if (length(union) > 2) {
    stop('T test not support group > 2')
  }
}
