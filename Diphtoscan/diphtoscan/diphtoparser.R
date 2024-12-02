
samples<-list.files()
if(length(which(samples=="diphtoscan_results.csv"))==1){
  samples<-samples[-which(samples=="diphtoscan_results.csv")]
  file.remove("diptoscan_results.csv")
}

if(exists("out")) rm(out)
for (i in 1:length(samples)) {
  dummy<-read.csv(paste(samples[i],"/",samples[i],".txt",sep = ""),sep = "\t")
  
  if(!exists("out")){
    out<-dummy
  }else{
    missing.dummy<-colnames(out)[-which(colnames(out) %in% colnames(dummy))]
    missing.out<-colnames(dummy)[-which(colnames(dummy) %in% colnames(out))]
    
    if(length(missing.dummy)>0){
      pad.d<-as.data.frame(matrix(NA, nrow = nrow(dummy), ncol = length(missing.dummy)))
      colnames(pad.d)<-missing.dummy
      dummy<-cbind(dummy, pad.d)
    }
    
    if(length(missing.out)>0){
      pad.o<-as.data.frame(matrix(NA, nrow = nrow(out), ncol = length(missing.out)))
      colnames(pad.o)<-missing.out
      out<-cbind(out, pad.o)
    }
    
    out<-rbind(out,dummy)
  }
  
}

colnames(out)[1]<-"Sample"

write.csv(out,"diphtoscan_results.csv", row.names = FALSE)



