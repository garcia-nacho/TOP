library(jsonlite)


input<-list.files(pattern = "*_rmlst.json")
df<-fromJSON(input)

output<-df$taxon_prediction

#output$species<-df$fields$species
output$genus<-df$fields$genus

colnames(output)<-paste("rMLST_",colnames(output),sep="")
output$rST<-df$fields$rST
output$Sample<-gsub(".*/","",gsub("_.*","",input))
write.csv(output, gsub("_rmlst.json","_rmlst.csv", input),row.names = FALSE )