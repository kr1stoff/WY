#!/usr/bin/Rscript
library(ggplot2)

args <- commandArgs(TRUE)

if (length(args) != 2) {
  usage = "Usage: Rscript insertsize.R  <picard insert_size.txt>  <out dir>\n"
  cat(usage)
  quit("no")
}


file_name <- args[1]
dir_out <- args[2]


data <- read.csv(file_name,sep = '\t',comment.char ="#",skip = 10)  #跳过前10行

# 确定x轴步长
group <- 30   # 分成30各组
record_diff <- max(data$insert_size) - min(data$insert_size)
x_gap <- round(record_diff/group)
x_step_num <- nchar(as.character(x_gap))
cat(x_step_num,'\n')
if (x_step_num == 4){
    x_gap <- round(x_gap/1000)*1000
}else if (x_step_num == 3){
    x_gap <- round(x_gap/100)*100
}else if (x_step_num == 2){
    x_gap <- round(x_gap/10)*10
}

# 确定y轴步长
group <- 10   # 分成10各组
record_diff <- max(data$All_Reads.fr_count) - min(data$All_Reads.fr_count)
y_gap <- round(record_diff/group)
y_step_num <- nchar(as.character(y_gap))

if (y_step_num == 5){
    y_gap <- round(y_gap/10000)*10000
}else if (y_step_num == 4){
    y_gap <- round(y_gap/1000)*1000
}else if (y_step_num == 3){
    y_gap <- round(y_gap/100)*100
}else if (y_step_num == 2){
    y_gap <- round(y_gap/10)*10
}


g <- ggplot(data, aes(x = insert_size, y = All_Reads.fr_count)) +
       geom_bar(stat = "identity",fill ="#4d97cd") +   # 条形图
  
       theme(panel.grid =element_blank()) +
       theme(axis.text.x = element_text(size = 3.5),axis.text.y = element_text(size = 3.5)) +  # 坐标标题大小
       theme(axis.title = element_text(size = 4)) +       # 标题大小
       theme(axis.ticks = element_blank()) +              # 删去所有刻度线
       theme(panel.border = element_blank()) +            # 删去外层边框
       theme(panel.border = element_rect(fill=NA,color="black", size=0.2))+   #重新添加外框线
       theme(panel.background = element_blank()) +        # 去掉画布背景灰色
       theme(axis.ticks=element_line(color="black",size=0.2,lineend = 1)) +   # 添加刻度线
       scale_x_continuous(name = "Insert Size Length(bp)" ,breaks = seq(0, max(data$insert_size)*1.1, x_gap)) +  # x轴刻度间隔
       scale_y_continuous(name = "Count" ,breaks = seq(0, max(data$All_Reads.fr_count)*1.2, y_gap))       # y轴刻度间隔


#图片保存
ggsave(
  paste(dir_out, "insertsize.png", sep = '/'),
  limitsize = FALSE,
  width = 4,
  height = 2
)

# dev.off()