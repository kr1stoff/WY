#!/usr/bin/env Rscript
suppressPackageStartupMessages({
library(tidyverse)
library(infotheo)})
######## Options ########
args <- commandArgs(trailingOnly = TRUE)
print_usage <- function() {
  cat("PROGRAM <depth file> <out_dir>\n")
  quit(save = "no", status = 1)
}
if (length(args) != 2) {
  cat("Need 2 Arguments!\n")
  print_usage()
} else if (!file.exists(args[1])) {
  cat("<depth file> not exists!\n")
  print_usage()
} 

samtools_depth_file <- args[1]
out_abs_prefix <- args[2]
######## Function #######
## coverage & uniformity
coverge_stats <- function(depth = 1) {
  tmp_filt_depth_df <- depth_df %>%
    filter(Depth >= depth)
  tmp_coverage_percent <- sprintf("%1.2f%%", nrow(tmp_filt_depth_df)/whole_genome_length * 100)
  return(tmp_coverage_percent)
}
######## Library ########

lancet_pal <- c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF","#ADB6B6FF","#1B1919FF")

# read depth file from `samtools depth -a`
depth_df <- read.table(samtools_depth_file, header = F)
colnames(depth_df) <- c("Ref","Position","Depth")

# bam stats
## basic 
whole_genome_length = nrow(depth_df)
mean_depth <- mean(depth_df$Depth)
uniformity_threshold <- 0.2 * mean_depth
coverage_depth <- round(sum(depth_df$Depth)/whole_genome_length, digits = 1)
coverage_depth <- format(coverage_depth, big.mark=",", scientific=FALSE)
coverage_percent <- coverge_stats()
coverage10_percent <- coverge_stats(10)
coverage30_percent <- coverge_stats(30)
coverage100_percent <- coverge_stats(100)
uniformity <- coverge_stats(uniformity_threshold)
bam_stats_df <- data.frame(
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
bin_cols <- discretize(depth_df$Position, "equalwidth", 100)
max_x_dvd100 <- ceiling(whole_genome_length/100)
depth_df$Bin <- bin_cols$X * max_x_dvd100
# xtick label
x_breaks <- c(1:5) * (max_x_dvd100 / 5 * 100)
x_labels <- paste0(x_breaks/1000, "K")

# mean depth of bins
mean_depth_df <- depth_df %>%
  group_by(Bin) %>%
  summarise(Mean=mean(Depth))

# main plot
p1 <- ggplot(mean_depth_df, aes(x=Bin, y=Mean)) +
  geom_bar(stat="identity", fill="#87CEEB") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  geom_hline(aes(yintercept = uniformity_threshold), linetype=5, col="#FF0000", size = 1) +
  labs(x = "Position", y = "Depth") +
  theme_bw()
ggsave(paste0(out_abs_prefix, ".genome_coverage_depth.png"), plot = p1, dpi = 300)
#~ ylim 1000
p2 <- ggplot(mean_depth_df, aes(x=Bin, y=Mean)) +
  geom_bar(stat="identity", fill="#87CEEB") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  geom_hline(aes(yintercept = uniformity_threshold), linetype=5, col="#FF0000", size = 1) +
  labs(x = "Position", y = "Depth") +
  coord_cartesian(ylim = c(0,1000)) +
  theme_bw()
ggsave(paste0(out_abs_prefix, ".genome_coverage_depth_ylim1000.png"), plot = p2, dpi = 300)
