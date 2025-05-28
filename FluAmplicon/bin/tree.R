#!/use/bin/env Rscript
# Title     : tree.R
# Objective : Phylogenetic Tree Plot
# Created by: MingJia
# Created on: 2022/6/30
options(warn = -1)
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ggtree)
  library(logging)
  library(stringr)
  library(dplyr)
})
basicConfig()

#### Some Functions ####
plot_tree <- function(obj,
                      branch.length = "none",
                      layout = "circular",
                      tree.size = 2,
                      tiplap.size = 2,
                      ...) {
  tiplap.hjust <- ifelse(branch.length == "none", 0,-.08)
  p <- ggtree(
    obj,
    branch.length = branch.length,
    layout = layout,
    size = tree.size,
    aes(color = color)
  ) +
    geom_nodelab(hjust = 1.6,
                 vjust = -.5,
                 size = 1.5) +
    hexpand(0.618) +
    geom_tiplab(hjust = tiplap.hjust,
                size = tiplap.size,
                fontface = "italic") +
    scale_color_manual(values = c(red = "red", black = "black"))
  return(p)
}

#### Commond line options ####
option_list <- list(
  make_option(c("--tree"),
              type = "character",
              help = "The tree file to plot"),
  make_option(c("--name"),
              type = "character",
              help = "The branch name to add color"),
  make_option(c("--title"),
              type = "character",
              help = "The title for the plot"),
  make_option(
    c("--fontsize"),
    default = 2,
    type = "numeric",
    help = "The fontsize for the label[default = %default]"
  ),
  make_option(
    c("--branchlength"),
    action = "store_true",
    default = FALSE,
    help = "whether show the branch length[default = %default]"
  ),
  make_option(
    c("--layout"),
    default = "roundrect",
    type = "character",
    help = "The tree plot layout[default = %default]"
  ),
  make_option(
    c("-p", "--prefix"),
    type = "character",
    default = "result",
    help = "The out prefix[default = %default]"
  )
)
opts <- parse_args(
  OptionParser(
    usage = "%prog[options]",
    option_list = option_list,
    description = "\nScript of Flu Phylogenetic Tree plot"
  ),
  positional_arguments = F
)

#### Main ####
if (!opts$branchlength) {
  opts$branchlength <- "none"
} else {
  opts$branchlength <- "branch.length"
}

loginfo("Load the tree file %s", opts$tree)
tree <- read.tree(opts$tree)

loginfo("Start to draw")
df.color <- data.frame(label = tree$tip.label, color = "black")
if (!is.null(opts$name)) {
  df.color <-
    df.color %>% mutate(color = replace(color, label == opts$name, "red"))
}
df.color$color <- as.factor(df.color$color)
tree <- full_join(tree, df.color, by = "label")
tree@data$color[is.na(tree@data$color)] <- "black"
# set min bootstrap value
tree@phylo$node.label <-
  ifelse(tree@phylo$node.label < 0.9, "", round(as.numeric(tree@phylo$node.label) *
                                                  100, 0))
p <- plot_tree(
  tree,
  branch.length = opts$branchlength,
  layout = opts$layout,
  tree.size = 0.618,
  tiplap.size = opts$fontsize,
)

if (!is.null(opts$title)) {
  p <- p + ggtitle(opts$title)
}

p <- p + theme_tree2(legend.position = 'none')
loginfo("Save the plot")
f_png <- str_c(opts$prefix, "png", sep = '.')
ggsave(f_png, p, limitsize = FALSE)