files.to.take<-list.files(pattern = ".tsv")

for (i in 1:length(files.to.take)) {
  dummy<-read.csv(files.to.take[i],sep = "\t")
  
  if(nrow(dummy)>0){
    dummy$DATABASE<-gsub("_","",dummy$DATABASE)
    output.temp<-as.data.frame(matrix(paste(dummy$GENE,collapse = ", "), nrow = 1, ncol = 1))
    colnames(output.temp)<-paste("Abricate_", unique(dummy$DATABASE),sep = "")
    output.temp$Sample<-gsub("_.*","", unique(dummy$X.FILE))
    if(dummy$DATABASE[1]== "ncbi"){
      output.temp$Abricate_Resistance<-paste(dummy$RESISTANCE,collapse = ", ")
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
  write.csv(output,"Abricate.csv",row.names = FALSE)



