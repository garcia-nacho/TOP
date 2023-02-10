files.to.take<-list.files(pattern = ".tsv")

for (i in 1:length(files.to.take)) {
  dummy<-read.csv(files.to.take[i],sep = "\t")
  
  if(nrow(dummy)>0){
    dummy$DATABASE<-gsub("_","",dummy$DATABASE)
    output.temp<-as.data.frame(matrix(paste(dummy$GENE[order(dummy$GENE)],collapse = " | "), nrow = 1, ncol = 1))
    colnames(output.temp)<-paste("Abricate_", unique(dummy$DATABASE),sep = "")
    output.temp$Sample<-gsub("_.*","", unique(dummy$X.FILE))
    
    if(dummy$DATABASE[1]== "ncbi"){
      output.temp$Abricate_Resistance<-paste(dummy$RESISTANCE[order(dummy$RESISTANCE)],collapse = " | ")
      
      for (l in 1:nrow(dummy)) {
       if(length(which(colnames(output.temp)==paste("AMR.",dummy$RESISTANCE[l], sep="") ))==0){
         output.temp$dummy<-NA
         output.temp$dummy<-dummy$GENE[l]
         colnames(output.temp)[which(colnames(output.temp)=="dummy")]<-paste("AMR.",dummy$RESISTANCE[l], sep="")
       }else{
         col.to.use<-which(colnames(output.temp)==paste("AMR.",dummy$RESISTANCE[l], sep="") )
         output.temp[,col.to.use]<-gsub(","," |",paste(output.temp[,col.to.use],dummy$GENE[l],sep = ", "))
       } 
      }
    }
    
    if(dummy$DATABASE[1]== "vfdb"){
      for (l in 1:nrow(dummy)) {
        factor<-paste("VF",gsub(" ","_",  gsub(" \\(.*","",gsub(".*- ","",gsub("\\] \\[.*","",dummy$PRODUCT[l])))),sep = "_")
        if(length(which(colnames(output.temp)==factor))==0){
          output.temp$dummy<-dummy$GENE[l]
          colnames(output.temp)[which(colnames(output.temp)=="dummy")]<-factor
        }else{
          col.to.use<-which(colnames(output.temp)==factor)
          output.temp[,col.to.use]<-gsub(","," |",paste(output.temp[,col.to.use],dummy$GENE[l],sep = ", "))
        }
      }
    }
    
    if(dummy$DATABASE[1] %in% c("HinfGyrSubA","HinfTopoIVSubA", "HinfFtsI" )){
      output.temp$dummy<-paste(dummy$RESISTANCE,collapse = " | ")
      colnames(output.temp)[which(colnames(output.temp)=="dummy")]<-paste("Mutations_",dummy$DATABASE[1],sep = "")
    }
    
    
    if(!exists("output")){
      output<-output.temp
    }else{
      
      output<-merge(output, output.temp,by="Sample")
    }

  }else{
    output.temp<-as.data.frame(matrix(NA, nrow = 1, ncol = 1))
    database<-gsub(".tsv","",gsub(".*_","",files.to.take[i]))
    colnames(output.temp)<-paste("Abricate_", database,sep = "")
    
    if(database== "ncbi"){
      output.temp$Abricate_Resistance<-NA
    }
    output.temp$Sample<-gsub("_.*","", files.to.take[i])
    
    if(!exists("output")){
      output<-output.temp
    }else{
      
      
    output<-merge(output, output.temp,by="Sample")
    }
}
}
if(length(which(output[1,]=="NA")>0) ) output[1,which(output[1,]=="NA")]<-NA
  write.csv(output,"Abricate.csv",row.names = FALSE)



