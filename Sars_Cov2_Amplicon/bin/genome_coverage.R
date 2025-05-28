#!/usr/bin/env Rscript

######## Options ########
args <- commandArgs(trailingOnly = TRUE)
print_usage <- function() {
  cat("PROGRAM <depth file> <summary_numbers_file> <out_prefix>\n")
  quit(save = "no", status = 1)
}
if (length(args) != 3) {
  cat("Need 3 Arguments!\n")
  print_usage()
} else if (!file.exists(args[1]) | !file.exists(args[2])) {
  cat("<depth file> or <summary_numbers_file> not exists!\n")
  print_usage()
}

samtools_depth_file <- args[1]
samtools_summary_numbers_file <- args[2]
out_abs_prefix <- args[3]

######## Library ########
library(grid)
library(tidyverse)
library(infotheo) # 无监督分箱包

lancet_pal <- c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF","#ADB6B6FF","#1B1919FF")

# read depth file from `samtools depth -a`
depth_df <- read.table(samtools_depth_file, header = F)
colnames(depth_df) <- c("Ref","Position","Depth")

# 比对率mapping rate
summary_numbers <- read.table(samtools_summary_numbers_file, sep = ":")
raw_total_seqs <- summary_numbers[which(summary_numbers$V1 == 'raw total sequences'),]$V2
reads_mapped <- summary_numbers[which(summary_numbers$V1 == 'reads mapped'),]$V2
mapped_rate <- sprintf("%1.2f%%", reads_mapped/raw_total_seqs*100)

# bam stats
## basic 
whole_genome_length <- diff(range(depth_df$Position))
mean_depth <- mean(depth_df$Depth)
uniformity_threshold <- 0.2 * mean_depth
## coverage & uniformity
coverge_stats <- function(depth = 1) {
  tmp_filt_depth_df <- depth_df %>% 
    filter(Depth >= depth)
  tmp_coverage_percent <- sprintf("%1.2f%%", dim(tmp_filt_depth_df)[1]/whole_genome_length * 100)
  return(tmp_coverage_percent)
}
coverage_depth <- round(sum(depth_df$Depth)/whole_genome_length, digits = 1)
coverage_depth <- format(coverage_depth, big.mark=",", scientific=FALSE)
coverage_percent <- coverge_stats()
coverage10_percent <- coverge_stats(10)
coverage30_percent <- coverge_stats(30)
coverage100_percent <- coverge_stats(100)
uniformity <- coverge_stats(uniformity_threshold)
bam_stats_df <- data.frame(
                          "比对率" = mapped_rate,
                          "平均深度" = coverage_depth,
                          "覆盖度" = coverage_percent,
                          "深度≥10x" = coverage10_percent,
                          "深度≥30x" = coverage30_percent,
                          "深度≥100x" = coverage100_percent,
                          "均一性" = uniformity,
                          check.names = F
                          )
write.table(bam_stats_df, 
            file = paste0(out_abs_prefix, ".bam_stats.txt"),
            row.names = F,
            sep = "\t",
            quote = F,
            fileEncoding = "UTF-8")

# barplot bin width
bin_num <- 300  # 设计分成300分
bin_cols <- discretize(depth_df$Position, "equalwidth", bin_num)  # 等宽分箱，
max_x_dvd100 <- ceiling(whole_genome_length/bin_num)
depth_df$Bin <- bin_cols$X * max_x_dvd100

# 设置x轴
x_lables_num <- 30
x_breaks <- c(1:x_lables_num) * (max_x_dvd100 / x_lables_num * bin_num)
x_labels <- paste0(x_breaks/1000, "K")

# 设置y轴标签数量
# y_lables_num <- 5
# max_depth <- max(depth_df$Depth)


# mean depth of bins
mean_depth_df <- depth_df %>%
  group_by(Bin) %>%
  summarise(Mean=mean(Depth))

# main plot
p1 <- ggplot(mean_depth_df, aes(x=Bin, y=Mean)) +
  geom_bar(stat="identity", fill="#87CEEB" ,width = 90) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  geom_hline(aes(yintercept = uniformity_threshold), linetype=5, col="#FF0000", size = 1) +
  labs(x = "Position", y = "Depth") +
  theme_bw()
ggsave(paste0(out_abs_prefix, ".genome_coverage_depth.png"), plot = p1, dpi = 300 ,units = 'mm',width = 400)


#~ ylim 1000
p2 <- ggplot(mean_depth_df, aes(x=Bin, y=Mean)) +
  geom_bar(stat="identity", fill="#87CEEB") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  geom_hline(aes(yintercept = uniformity_threshold), linetype=5, col="#FF0000", size = 1) +
  labs(x = "Position", y = "Depth") +
  coord_cartesian(ylim = c(0,1000)) +
  theme_bw()
ggsave(paste0(out_abs_prefix, ".genome_coverage_depth_ylim1000.png"), plot = p2, dpi = 300,units = 'mm',width = 400)
