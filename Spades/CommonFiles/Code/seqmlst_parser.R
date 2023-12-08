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
if(nrow(results)>1) results<-results[1,]
}

if(length(which(document$name==tolower(gsub(" .*","",results$rMLST_taxon))))>0){
dumm.db<-document$databases[[which(document$name==tolower(gsub(" .*","",results$rMLST_taxon)))]]
}

if(exists("dumm.db")){
  
  if(length(grep("Escherichia",results$rMLST_taxon))==1){
    isolatesurl<-dumm.db$href[grep("isolates",dumm.db$description)]
    sequrl<-dumm.db$href[grep("sequence",dumm.db$description)]  
  }else{
    isolatesurl<-dumm.db$href[grep(paste(results$rMLST_taxon,"isolates"),dumm.db$description)]
    sequrl<-dumm.db$href[grep(paste(results$rMLST_taxon,"sequence"),dumm.db$description)]  
  }
  
  if(length(sequrl)==0){
    isolatesurl<-dumm.db$href[grep("isolates",dumm.db$description)]
    sequrl<-dumm.db$href[grep("sequence",dumm.db$description)]  
  }
  
  schemes.url<-paste(sequrl,"/schemes",sep = "")
  system(paste("GET ", schemes.url, " > ", "schemes.json",sep = ""))
  
  if(file.exists("schemes.json")){
    schemes <- fromJSON(txt="schemes.json")
    sequrl<- schemes$schemes$scheme[which(schemes$schemes$description=="MLST")]
    if(length(grep("Escherichia",results$rMLST_taxon))==1) sequrl<- schemes$schemes$scheme[which(schemes$schemes$description=="MLST (Achtman)")]
    
    system("rm schemes.json")
  }else{
    rm(sequrl)
  }
  
}else{
  sequrl<-vector()
}
  
#write(paste(sequrl,"/schemes/1/sequence",sep = ""), stdout())
inputfasta<-list.files(pattern = "clean_contigs.fasta")
#system(paste("/media/nacho/Data/DockerImages/TOP/Spades/CommonFiles/Code/REST_Runner.sh",inputfasta,paste(sequrl,"/sequence",sep = ""),"dummy.json"))
if(length(sequrl)==1){

system(paste("/home/docker/CommonFiles/Code/REST_Runner.sh",inputfasta,paste(sequrl,"/sequence",sep = ""),"dummy.json"))
input<-list.files(pattern = "dummy.json")
df<-fromJSON(input)

  if(!is.null(df$fields$ST)){
  output<-as.data.frame(df$fields$ST )
  colnames(output)<-"MLST.Type"
  output$ClonalComplex<-df$fields$clonal_complex
  sch<-vector()

  for (i in 1:length(df$exact_matches)) {
    sch<-c(sch,paste(names(df$exact_matches)[i],":",df$exact_matches[[i]]$allele_id,sep = ""))
    #output$dummy<-paste(df$exact_matches[[i]]$allele_id,collapse = "/")
    #colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
  }
  
  }else{
  output<-as.data.frame(NA)  
  colnames(output)<-"MLST.Type"
  output$ClonalComplex<-df$fields$clonal_complex
  sch<-vector()
  if(length(df$exact_matches)>0){ 

    for (i in 1:length(df$exact_matches)) {
      sch<-c(sch,paste(names(df$exact_matches)[i],":",df$exact_matches[[i]]$allele_id,sep = ""))
      #output$dummy<-paste(df$exact_matches[[i]]$allele_id,collapse = "/")
      #colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
    }
  }else{
      sch<-NA
    }
  }

output$MLST.Scheme<-paste(sch[order(sch)],collapse = " | ")
if(is.na(sch)) output$MLST.Scheme<-NA

}else{
  output<-as.data.frame(NA)  
  colnames(output)<-"MLST.Type"
}

output$Sample<-gsub(".*/","",gsub("_clean_contigs.fasta","",inputfasta))
output$MLST_Date<-date

shortname<-results$rMLST_taxon
shortname<-paste(unlist(base::strsplit(gsub(" .*", "",shortname),""))[1] ,
      paste(unlist(base::strsplit(gsub(".* ", "",shortname),""))[c(1:3)],collapse = ""),sep = "")


if(length(grep("Salmonella",results$rMLST_taxon))==1) shortname<-"Salmo"
if(length(grep("Mycobact",results$rMLST_taxon))==1) shortname<-"Myco"

#colnames(output)[-which(colnames(output) %in% c("Sample","MLST_Date"))]<-paste(shortname, colnames(output)[-which(colnames(output) %in% c("Sample","MLST_Date"))],sep = "_")

write.csv(output, gsub("_clean_contigs.fasta","_seqmlst.csv", inputfasta),row.names = FALSE )
if(exists("input")) file.rename(input, gsub("_clean_contigs.fasta","_seqmlst.json", inputfasta))
write.table(shortname, paste(shortname,".agent",sep=""), row.names =FALSE, col.names=FALSE, quote = FALSE )

write.table(results$rMLST_genus, paste(gsub("/","_",results$rMLST_genus),".genus",sep=""), row.names =FALSE, col.names=FALSE, quote = FALSE )
