#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
print_usage <- function() {
  cat("PROGRAM <tree file>\n")
  quit(save = "no", status = 1)
}

#判断传入参数数量及文件是否存在
if (length(args) != 1) {
  cat("Just need 1 Arguments!\n")
  print_usage()
} else if (!file.exists(args[1])) {
  cat("<tree file> not exists!\n")
  print_usage()
} 

#newick格式
tree_file <- args[1]

# dirname()获取文件路径
output_dir <- dirname(tree_file)

library(ggplot2) # 加载ggplot2
library(ggtree) # 加载ggtree
tree=read.tree(tree_file) # 读取nwk文件
data=fortify(tree)

tregraph=ggtree(tree, layout="rectangular", size=0.8, col="deepskyblue3") +
  # 树体：树文件、树形、粗细、颜色
  geom_tiplab(size=3, color="purple4", hjust=-0.05) +
  # 枝名：大小、颜色、高度
  geom_tippoint(size=1.5, color="deepskyblue3") +
  # 端点：大小、颜色
  geom_nodepoint(color="orange", alpha=1/4, size=2) +
  # 末节点：颜色、透明度、大小
  theme_tree2() +
  # x轴标尺
  xlim(NA, max(data$x)*1.3)
  # x轴宽度
  
pdf("tree_graph.pdf")
tregraph
dev.off()
