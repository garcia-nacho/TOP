library(jsonlite)

#DB update
rmlst.res<-list.files(pattern = "_rmlst.csv")
try(system("GET  https://rest.pubmlst.org/ > current_db.json"))

if(file.exists("current_db.json")){
  
  document <- fromJSON(txt="current_db.json")
  date<-gsub("-","",Sys.Date())
  system("rm current_db.json")
}else{
  document <- fromJSON(txt="/home/docker/CommonFiles/pubmlstDB.json")
  date<-"20221213"
}

document<-as.list(document)
results<-read.csv(rmlst.res)

if(nrow(results)>1){
results$MLST_Conflict<-paste(paste(results$rMLST_support, results$rMLST_taxon,sep = "% "), collapse = " / ")
results<-results[which(results$rMLST_support==max(results$rMLST_support)),]
}

dumm.db<-document$databases[[which(document$name==tolower(gsub(" .*","",results$rMLST_taxon)))]]

isolatesurl<-dumm.db$href[grep(paste(results$rMLST_taxon,"isolates"),dumm.db$description)]
sequrl<-dumm.db$href[grep(paste(results$rMLST_taxon,"sequence"),dumm.db$description)]

#write(paste(sequrl,"/schemes/1/sequence",sep = ""), stdout())
inputfasta<-list.files(pattern = "clean_contigs.fasta")
#system(paste("/media/nacho/Data/DockerImages/TOP/Spades/CommonFiles/Code/REST_Runner.sh",inputfasta,paste(sequrl,"/schemes/1/sequence",sep = ""),"dummy.json"))
if(length(sequrl)==1){
system(paste("/home/docker/CommonFiles/Code/REST_Runner.sh",inputfasta,paste(sequrl,"/schemes/1/sequence",sep = ""),"dummy.json"))

input<-list.files(pattern = "dummy.json")
df<-fromJSON(input)

  if(!is.null(df$fields$ST)){
  output<-as.data.frame(df$fields$ST )
  colnames(output)<-"ST"
  output$ClonalComplex<-df$fields$clonal_complex
  for (i in 1:length(df$exact_matches)) {
    output$dummy<-df$exact_matches[[i]]$allele_id
    colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
  }
  }else{
  output<-as.data.frame(NA)  
  colnames(output)<-"ST"
  output$ClonalComplex<-df$fields$clonal_complex
  for (i in 1:length(df$exact_matches)) {
    output$dummy<-df$exact_matches[[i]]$allele_id
    colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
  }
  }

}else{
  output<-as.data.frame(NA)  
  colnames(output)<-"ST"
}


output$Sample<-gsub(".*/","",gsub("_clean_contigs.fasta","",inputfasta))
output$MLST_Date<-date

shortname<-results$rMLST_taxon
shortname<-paste(unlist(base::strsplit(gsub(" .*", "",shortname),""))[1] ,
      paste(unlist(base::strsplit(gsub(".* ", "",shortname),""))[c(1:3)],collapse = ""),sep = "")

colnames(output)[-which(colnames(output) %in% c("Sample","MLST_Date"))]<-paste(shortname, colnames(output)[-which(colnames(output) %in% c("Sample","MLST_Date"))],sep = "_")

write.csv(output, gsub("_clean_contigs.fasta","_seqmlst.csv", inputfasta),row.names = FALSE )
if(exists("input")) file.rename(input, gsub("_clean_contigs.fasta","_seqmlst.json", inputfasta))
write.table(shortname, paste(shortname,".agent",sep=""), row.names =FALSE, col.names=FALSE, quote = FALSE )
write.table(results$rMLST_genus, paste(results$rMLST_genus,".genus",sep=""), row.names =FALSE, col.names=FALSE, quote = FALSE )
