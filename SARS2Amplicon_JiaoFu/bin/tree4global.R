#!/usr/bin/env Rscript

######## Options ########
args <- commandArgs(trailingOnly = TRUE)
print_usage <- function() {
  cat("PROGRAM <tree file>\n")
  quit(save = "no", status = 1)
}
if (length(args) != 1) {
  cat("Just need 1 Arguments!\n")
  print_usage()
} else if (!file.exists(args[1])) {
  cat("<tree file> not exists!\n")
  print_usage()
} 

tree_file <- args[1]
output_dir <- dirname(tree_file)
sample_name <- unlist(strsplit(basename(tree_file), split = "[.]"))[1]

# check tree file
if (file.info(tree_file)$size == 0) {
  print("Failed to create phylogenetic tree!")
  system(paste0("touch ",output_dir,"/Failed_to_create_phylogenetic_tree"), intern = TRUE)
}

######## Libraries & Settings ########
library(ggtree)
library(ggplot2)
library(cowplot)
library(hash)
nrc_pal <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")
lancet_pal <- c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF","#ADB6B6FF","#1B1919FF")

######## Main ########
tree <- read.tree(tree_file)
# [20220608 update] 新冠全球库，输入样本标记
groupInfo <- list(CUR=sample_name, REF=tree$tip.label[-match(sample_name, tree$tip.label)])
tree_withGroup <- groupOTU(tree, groupInfo)
data <- fortify(tree_withGroup)

######## Plots ########
## calculate longest node level
tree_hash <- hash(keys=data$node, values=data$parent)
node_length_vector <- c()
for (key in names(tree_hash)) {
  tmp_node <- key
  tmp_nodes <- c(tmp_node)
  while (tmp_node != as.character(values(tree_hash[tmp_node]))) {
    tmp_node <- as.character(values(tree_hash[tmp_node]))
    tmp_nodes <- c(tmp_nodes, tmp_node)
  }
  node_length_vector <- c(node_length_vector, length(tmp_nodes))
}
longest_node_level <- max(node_length_vector) + 1
# [20220606 update] 1.xlim n-fold, 2.ggsave width
#~ p1.rectangular none length, p2.rectangular, p3.circular, p4.slanted
max_nchar_label <- max(nchar(data$label), na.rm = T) # 名字长度
if (max_nchar_label <= 30) {
  ggtree_xlim_right <- longest_node_level * 1.2
  p2_xlim_fold <- 1.6
  p14_width <- 10
  p23_width <- 10
} else if (max_nchar_label > 30 & max_nchar_label <= 50) {
  ggtree_xlim_right <- longest_node_level * 1.4
  p2_xlim_fold <- 1.8
  p14_width <- 13
  p23_width <- 15
} else {
  ggtree_xlim_right <- longest_node_level * 1.4
  p2_xlim_fold <- 2
  p14_width <- 15
  p23_width <- 18
} 
# graph height
sample_number = length(which(data$isTip == TRUE))
if (sample_number <= 20) {
  graph_height <- 5.7
} else if (sample_number > 20 & sample_number <= 100) {
  graph_height <- 13
} else {
  graph_height <- 25
}

## branch.length scientific
branch_len_format <- format(data$branch, scientific=T, digits=2)
#~ 1 branch.length = "none"
p1 <- ggtree(tree_withGroup, layout="rectangular", 
             branch.length = "none",
             aes(color=group)) + #
  geom_tiplab(hjust=-0.05, size=4.5) +
  geom_tippoint(color="black", alpha=1/4) +
  geom_nodepoint(color="black", alpha=1/4) +
  scale_color_manual(values = c("red", "black")) +
  theme_tree2() +
  guides(color="none") +
  xlim(NA, ggtree_xlim_right) +
  geom_text(aes(x=branch, label=branch_len_format, vjust=-.3), color = "black")
ggsave(plot = p1, filename = "rectangular.svg",
       width = p14_width, height = graph_height, path = output_dir)
#~ 2 branch.length = "branch.length"
p2 <- ggtree(tree_withGroup, layout="rectangular",
			 aes(color=group)) +
  geom_tiplab(hjust=-0.05, size=4.5) +
  geom_tippoint(color="black", alpha=1/4) +
  geom_nodepoint(color="black", alpha=1/4) +
  scale_color_manual(values = c("red", "black")) +
  theme_tree2() +
  guides(color="none") +
  xlim(NA, max(data$branch)*p2_xlim_fold)
ggsave(plot = p2, filename = "rectangular_bl.svg",
       width = p23_width, height = graph_height, dpi = 300, path = output_dir) #width = 9
#~ 4 slanted
p4 <- ggtree(tree_withGroup, layout="slanted", branch.length = "none",
  			 aes(color=group)) + #
  geom_tiplab(hjust=-0.05, size=4.5) +
  geom_tippoint(color="black", alpha=1/4) +
  geom_nodepoint(color="black", alpha=1/4) +
  scale_color_manual(values = c("red", "black")) +
  theme_tree2() +
  guides(color="none") +
  xlim(NA, ggtree_xlim_right)
ggsave(plot = p4, filename = "slanted.svg",
       width = p14_width, height = graph_height, dpi = 300, path = output_dir) #width = 9
#~ 3 circular
p3 <- ggtree(tree_withGroup, layout="circular",
             branch.length = "none", root.position = 10,
             aes(color=group)) + #
  geom_tiplab(hjust=-0.05, size=4.5) + # angle=0
  geom_tippoint(color="black", alpha=1/4) +
  geom_nodepoint(color="black", alpha=1/4) + 
  xlim(NA, ggtree_xlim_right*1.4) +
  scale_color_manual(values = c("red", "black")) +
  theme_tree2() +
  guides(color="none")
ggsave(plot = p3, filename = "circular.svg",
       width = p23_width, height = p23_width, dpi = 300, path = output_dir) #width = 9
#~ multi graphs
# pm <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2)
# ggsave(plot = pm, filename = "jagsaw.svg",
#        width = pm_width, height = graph_height, dpi = 300, path = output_dir) #width = 15
