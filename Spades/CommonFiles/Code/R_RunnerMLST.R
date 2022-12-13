library(jsonlite)
inputfasta<-list.files(pattern = "clean_contigs.fasta")
system(paste("/home/docker/CommonFiles/Code/mlstRunner.sh", inputfasta))

input<-list.files(pattern = "testsh_rmlst.json")
df<-fromJSON(input)


output<-as.data.frame(df$fields$ST )
colnames(output)<-"ST"
output$ClonalComplex<-df$fields$clonal_complex

for (i in 1:length(df$exact_matches)) {
  output$dummy<-df$exact_matches[[i]]$allele_id
  colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
}

output$Sample<-gsub(".*/","",gsub("_clean_contigs.fasta","",inputfasta))
output$MLST_Date<-gsub("-","",Sys.Date())

write.csv(output, gsub("_clean_contigs.fasta","_seqmlst.csv", inputfasta),row.names = FALSE )
file.rename(input, gsub("_clean_contigs.fasta","_seqmlst.json", inputfasta))
