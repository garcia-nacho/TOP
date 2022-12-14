library(jsonlite)


#System to get MLST


input<-list.files(pattern = "*_mlst.json")
df<-fromJSON(input)


output<-as.data.frame(df$fields$ST )
colnames(output)<-"ST"
output$ClonalComplex<-df$fields$clonal_complex

for (i in 1:length(df$exact_matches)) {
  output$dummy<-df$exact_matches[[i]]$allele_id
  colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
}

output$Sample<-gsub(".*/","",gsub("_.*","",input))
write.csv(output, gsub("_mlst.json","_mlst.csv", input),row.names = FALSE )