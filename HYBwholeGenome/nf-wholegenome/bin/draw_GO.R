#!/usr/bin/Rscript
# func: draw go function classification spot plot

args <- commandArgs(TRUE)
if (length(args) != 2) {
  usage = "Usage:   Rscript draw_GO.R <go_stat.txt> <out_dir>\n"
  cat(usage)
  quit("no")
}


in_file <- args[1]
dir_out <- args[2]

library(ggplot2)
library(plyr)

GO_table <- read.table(in_file,header = T,quote = '"', check.names = F,sep = "\t")
names(GO_table) <- c("V1","V2","V3","V4")

# order tab
GO_table <- arrange(GO_table, desc(V2), desc(V3))
level <- GO_table$V3
GO_table$V3<- factor(GO_table$V3, levels=level)


color_list <- c("#db6968","#4d97cd","#f8984e","#459943",
                "#e8c559","#a3d393","#99cbeb","#fdc58f","#767676FF")

ggplot(GO_table, aes(x = V3 , y = V4 , fill = V2)) +
  # 条形图函数：stat表明取用样本点对应纵轴值
  # position_dodge(0.5) 表示同组内不同系列间错位比例 ：0.1表示90%重叠 
  geom_bar(stat = "identity" ,position = position_dodge(0.9),width =0.8,colour="black",size=0.2 )+
  scale_fill_manual(values= color_list)+
  geom_text(aes(label = V4 , y = V4 ), size=2.5, hjust = -0.4, position = position_dodge(0.5))+   #数据标签
  #去除绘图区和X轴之间的gap
  scale_y_continuous(expand = c(0, 0),limits = c(0,max(GO_table$V4)*1.1)) +
  labs(x = "GO Functional Classification", y = "Genes Count")+
  coord_flip()+
  theme_bw()+
  guides(color = guide_legend(ncol = 1),fill = guide_legend(reverse = F))+  #显示单列
  theme(legend.title=element_blank(), # 图例标题
        legend.background =element_blank(),# 图例背景
        legend.position = c("right"), #top,right,left, bottom
        legend.text = element_text( size = 6),  #图例字体
        legend.key.size = unit(14, "pt"),      #图例大小
        panel.grid = element_blank(),
        axis.text = element_text(size = 6)) 

# Out put
ggsave(
  paste(dir_out, "GO_anno_stats_level2.png", sep = '/'),
  width = 12,
  height = 7
)
