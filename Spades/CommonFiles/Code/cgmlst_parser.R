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

dumm.db<-document$databases[[which(document$name==tolower(gsub(" .*","",results$rMLST_taxon)))]]

if(length(grep("Escherichia",results$rMLST_taxon))==1){
  isolatesurl<-dumm.db$href[grep("isolates",dumm.db$description)]
  sequrl<-dumm.db$href[grep("sequence",dumm.db$description)]  
}else{
  isolatesurl<-dumm.db$href[grep(paste(results$rMLST_taxon,"isolates"),dumm.db$description)]
  sequrl<-dumm.db$href[grep(paste(results$rMLST_taxon,"sequence"),dumm.db$description)]  
} 

schemes.url<-paste(sequrl,"/schemes",sep = "")
system(paste("GET ", schemes.url, " > ", "schemes.json",sep = ""))

if(file.exists("schemes.json")){
  schemes <- fromJSON(txt="schemes.json")
  sequrl<- schemes$schemes$scheme[which(schemes$schemes$description=="cgMLST")]
  system("rm current_db.json")
}else{
  rm(sequrl)
}

system("GET https://rest.pubmlst.org/db/pubmlst_escherichia_seqdef/schemes/6/profiles_csv > profiles.csv" )
a<-read.csv("/media/nacho/Data/temp/test/profiles.list",sep = "\t")

#write(paste(sequrl,"/schemes/1/sequence",sep = ""), stdout())
inputfasta<-list.files(pattern = "clean_contigs.fasta")
#system(paste("/media/nacho/Data/DockerImages/TOP/Spades/CommonFiles/Code/REST_Runner.sh",inputfasta,paste(sequrl,"/sequence",sep = ""),"dummy.json"))
if(length(sequrl)==1){
  system(paste("/home/docker/CommonFiles/Code/REST_Runner.sh",inputfasta,paste(sequrl,"/sequence",sep = ""),"dummy.json"))
  
  input<-list.files(pattern = "dummy.json")
  df<-fromJSON(input)
  db<-as.data.frame(db)
  
  intersector<-list()
  pb<-txtProgressBar(min = 1, max=length(df$exact_matches), initial = 1)
  for (i in 1:length(df$exact_matches)) {
    setTxtProgressBar(pb,i)
    intersector[[i]]<- db$ST[which(db[,which(colnames(db)==names(df$exact_matches)[[i]])]== as.numeric(df$exact_matches[[i]]$allele_id))]
  }
  
  intersector.tab<-as.data.frame(table(unlist(intersector)))
  #Test enterobase
  #Test how it works
  
  if(!is.null(df$fields$ST)){
    output<-as.data.frame(df$fields$ST )
    colnames(output)<-"cgMLST.Type"
    output$ClonalComplex<-df$fields$clonal_complex
    sch<-vector()
    for (i in 1:length(df$exact_matches)) {
      sch<-c(sch,paste(names(df$exact_matches)[i],":",df$exact_matches[[i]]$allele_id,sep = ""))
      #output$dummy<-paste(df$exact_matches[[i]]$allele_id,collapse = "/")
      #colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
    }
  }else{
    output<-as.data.frame(NA)  
    colnames(output)<-"cgMLST.Type"
    output$ClonalComplex<-df$fields$clonal_complex
    sch<-vector()
    for (i in 1:length(df$exact_matches)) {
      sch<-c(sch,paste(names(df$exact_matches)[i],":",df$exact_matches[[i]]$allele_id,sep = ""))
      #output$dummy<-paste(df$exact_matches[[i]]$allele_id,collapse = "/")
      #colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
    }
  }
  #output$MLST.Scheme<-paste(sch[order(sch)],collapse = " | ")
}else{
  output<-as.data.frame(NA)  
  colnames(output)<-"cgMLST.Type"
}

output$Sample<-gsub(".*/","",gsub("_clean_contigs.fasta","",inputfasta))
output$cgMLST_Date<-date

write.csv(output, gsub("_clean_contigs.fasta","_cgmlst.csv", inputfasta),row.names = FALSE )
if(exists("input")) file.rename(input, gsub("_clean_contigs.fasta","_cgmlst.json", inputfasta))

