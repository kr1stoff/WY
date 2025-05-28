args = commandArgs(trailingOnly=TRUE)
#rm(list=ls())
#args[1]
options(stringsAsFactors=F)
query_strain<-args[1]
query_strain<-strsplit(query_strain,'/')[[1]]
query_strain<-query_strain[length(query_strain)]
#query_strain<-'test'

wbeO1.coverage<-read.delim(paste0(args[2],"/alignment_res/wbeO1.coverage"),header=T)
wbfO139.coverage<-read.delim(paste0(args[2],"/alignment_res/wbeO139.coverage"),header=T)
output.tb<-as.data.frame(array('uncertain',dim=c(1,2)))
colnames(output.tb)<-c('Sample','Serogroup')
output.tb[1,1]<-query_strain
serogroup<-c()
if(wbeO1.coverage$coverage>80){
  serogroup<-c(serogroup,'O1')
}
if(wbfO139.coverage$coverage>80){
  serogroup<-c(serogroup,'O139')
}
serogroup<-paste0(serogroup,collapse=',')
if(nchar(serogroup)>0){
  output.tb[1,2]<-serogroup
}

write.table(output.tb,paste0(args[2],"/result.tsv"),quote=F,col.names=T,row.names=F,sep = "\t")

