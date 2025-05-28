#!/usr/bin/Rscript
# func: draw CAZy function classification barplot

args <- commandArgs(TRUE)

if (length(args) != 2) {
  usage = "Usage:   Rscript draw_CAZy.R <CAZY.kind_stat.txt> <dir_out>\n"
  cat(usage)
  quit("no")
}

in_file <- args[1]
dir_out <- args[2]

library(ggplot2)

df <- read.table(in_file,header = TRUE ,quote = '"', check.names = F,sep = "\t")
names(df) <- c("kind","frequency")
len <- length(df$kind)
group <- c(seq(0.1,0.1*len,0.1))
scale_label <- c(df$kind)
color_list <- c("#db6968","#4d97cd","#f8984e","#459943",
                "#e8c559","#a3d393","#99cbeb","#fdc58f","#767676FF")


ggplot(df, aes(x = kind , y = frequency , fill = kind)) +
    geom_bar(stat = "identity" , hjust = 0.5,position = position_dodge(0),width =0.3,colour="black",size=0.2 )+
    scale_fill_manual(values= color_list)+


    #去除绘图区和X轴之间的gap
    scale_y_continuous(expand = c(0, 0),limits = c(0,max(df$frequency)*1.1)) +
    labs(x = "CAZy Functional Classification", y = "Genes Count")+
    theme_bw()+
    geom_text(aes(label = frequency), size=2.5, vjust = -0.4, position = position_dodge(0.5))+   #数据标签
    theme(panel.grid = element_blank(),           #网格线
          legend.title=element_blank(), # 图例标题
          legend.background =element_blank(),# 图例背景
          axis.title = element_text(size = 8),    #坐标轴标题
          axis.text = element_text(size = 5),      #刻度大小
          legend.text = element_text( size = 5),  #图例字体
          legend.key.size = unit(10, "pt"),      #图例大小
        )       

# Out put
ggsave(
  paste(dir_out, "CAZy_anno_stats.png", sep = '/'),
  width = 4.5,
  height = 4.5
)
