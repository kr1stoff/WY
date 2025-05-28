library(ggplot2)
library(grid)
library("optparse")

option_list = list(
  make_option(c("-i","--input"), type = "character",help = "Profile of each level"),
  make_option(c("-o","--output"), type = "character",help = "The output PDF file")
)
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

input <- opt$options$input
output <- opt$options$output

options(stringsAsFactors=F)
color=c("#ce0351", "#6600ff","red", "#00cfff", "#00ff00", "#ffcf00","red","yellow")
library(ggplot2)
library(grid)
library(RColorBrewer)
options(stringsAsFactors=F)
pdf(output, width=15, height=10)
d=read.table(input, sep='	')
df=data.frame(table(d$V2, d$V1)) # a new matrix
df$Freq[df$Freq==0]=NA # value = 0 -> value = NA
df=na.omit(df) # filter NA, for more beautiful pdf, but not necessary
df$Var1=factor(df$Var1, levels=rev(df$Var1)) # 
order1<-df[order(-df$Freq),]
order2<-order1[order(order1$Var2),]
order2$Var1<-factor(order2$Var1,levels =order2$Var1)
colours=colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFE528","#A65628","#F781BF","#999999"))(9)[1:length(levels(order2$Var2))]
ggplot(df,aes(x=Var1,y=Freq,fill=Var2)) + scale_fill_manual(values=colours)+
  geom_bar(position="stack",stat="identity") +theme(panel.background=element_rect(fill=NA,colour="grey"), panel.grid=element_line(color='grey'), panel.border=element_rect(fill='transparent',color='black'),legend.title=element_blank(), legend.text = element_text(size = 15), plot.title=element_text(face='bold', size=20),axis.title=element_text(size=20), axis.text.x=element_text(color='black', size=6),axis.text.y=element_text(color='black', size=12))+coord_flip()+labs(x='', y='Number of Genes')+geom_text(aes(label=Freq), hjust=-0.5, vjust=0.5, size = 3)+scale_x_discrete(limits = rev(levels(order2$Var1)))
dev.off()