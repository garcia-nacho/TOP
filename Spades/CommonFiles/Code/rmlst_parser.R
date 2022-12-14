library(jsonlite)


input<-list.files(pattern = "*_rmlst.json")
df<-fromJSON(input)

output<-df$taxon_prediction

#output$species<-df$fields$species
if(!is.null(df$fields$genus)){ 
output$genus<-df$fields$genus
}else{
output$genus<-NA   
}
colnames(output)<-paste("rMLST_",colnames(output),sep="")

if(!is.null(df$fields$rST)){ 
output$rST<-df$fields$rST
}else{
output$rST<-NA
}

output$Sample<-gsub(".*/","",gsub("_.*","",input))
write.csv(output, gsub("_rmlst.json","_rmlst.csv", input),row.names = FALSE )