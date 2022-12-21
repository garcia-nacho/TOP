files.to.take<-list.files(pattern = ".tsv")

if(length(files.to.take)>0){
for (i in 1:length(files.to.take)) {
  dummy<-read.csv(files.to.take[i],sep = "\t")
  
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

}
  write.csv(output, "Abricate.csv", row.names = FALSE)
  
}else{
  write.csv(NA,"Abricate.csv", row.names = FALSE )
}

