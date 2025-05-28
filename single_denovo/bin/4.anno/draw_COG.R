#!/usr/bin/Rscript
# func: draw cog function classification barplot

args <- commandArgs(TRUE)

if (length(args) != 2) {
  usage = "Usage:   Rscript draw_COG.R <all.COG.class.xls> <dir_out>\n"
  cat(usage)
  quit("no")
}


file_name <- args[1]
dir_out <- args[2]

library(ggplot2)

# parse the data
#Code         FunctionalCategories             GeneNumber      ColorCode
#  A       RNA processing and modification        0             #FF9900

data = read.table(file_name,header = T,comment.char = "",sep = "\t",check.names = F)
porder = factor(as.integer(rownames(data)), labels = data$Code)
#output 'porder':1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25

#paste:link the 'Code' and 'FunctionalCategories' with ':'
colors = paste(data$ColorCode)

data$title = paste(data$Code, data$FunctionalCategories, sep = ": ")
#output: 'A:RNA processing and modification'
desc = paste(data$title)


# draw
ggplot(data, aes(porder, GeneNumber)) +
  geom_bar(stat = "identity", aes(fill = porder)) +
  scale_fill_manual(values = colors, labels = desc) +
  labs(x = "",
       y = "Number of Genes",
       title = "COG function classification",
       fill = "") +
  geom_text(aes(label = GeneNumber),
            size = 2,
            vjust = -1,
            color = "black") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "#000000", size = 0.8),
        axis.text = element_text(color = "#000000", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(hjust = 1, vjust = 0.5),
        axis.title = element_text(color = "#000000", size = 14, face = "plain"),
        axis.title.x = element_text(margin = margin(2.5, 0, 2.5, 0, "mm")),
        axis.title.y = element_text(margin = margin(0, 2.5, 0, 2.5, "mm")),
        axis.ticks = element_line(color = "#000000", size = 0.5),
        axis.ticks.length = unit(0.1, 'cm'),

        legend.title = element_blank(),
        plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
        plot.margin = unit(c(5, 5, 5, 5), "mm")
  )


ggsave(
  paste(dir_out, "all.COG.bar.png", sep = '/'),
  width = 10,
  height = 7
)
# ggsave(
#   paste(dir_out, "all.COG.bar.pdf", sep = '/'),
#   width = 10,
#   height = 7
# )

