args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
f_chewbbaca <- args[1]
f_scheme <- args[2]
f_out <- args[3]

chewbbaca_res<-read.delim(f_chewbbaca,header=T)
colnames(chewbbaca_res)<-gsub('.fasta','',colnames(chewbbaca_res))

profile.res<-read.delim(f_scheme,header=T)

profile.res.tempt<-profile.res
for(i in 2:ncol(chewbbaca_res)){
  if(is.numeric(chewbbaca_res[1,i])){
    if(nrow(profile.res.tempt)>0){
      profile.res.tempt<-profile.res.tempt[(profile.res.tempt[,colnames(chewbbaca_res)[i]]==chewbbaca_res[1,i]),]
    }
  }
}


output.tb<-as.data.frame(array('-',dim=c(1,2)))
colnames(output.tb)<-c('Sample','Possible cgST')
output.tb[1,1]<-"Query  "
if(nrow(profile.res.tempt)>0){
  output.tb[1,2]<-paste0(profile.res.tempt$cgST,collapse=',')
}

write.table(output.tb,f_out,quote=F,col.names=T,row.names=F,sep = "\t")

