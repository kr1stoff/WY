args = commandArgs(trailingOnly=TRUE)
#rm(list=ls())
#args[1]
options(stringsAsFactors=F)
query_strain<-args[1]
#query_strain<-'test'

ctx.coverage<-read.delim(paste0(args[2],"/alignment_res/ctx.coverage"),header=T)
output.tb<-read.delim(paste0(args[2],"/result.tsv"),header=T)
#output.tb<-as.data.frame(array('ctx negative',dim=c(1,2)))
#colnames(output.tb)<-c('Sample','ctx')
#output.tb[1,1]<-query_strain
output.tb$ctx<-'uncertain'
if(length(which(ctx.coverage$coverage>80))>0){
  output.tb$ctx<-'positive'
}

write.table(output.tb,paste0(args[2],"/result.tsv"),quote=F,col.names=T,row.names=F,sep = "\t")

