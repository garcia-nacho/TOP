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

results<-results[1,]
dumm.db<-document$databases[[which(document$name==tolower(gsub(" .*","",results$rMLST_taxon)))]]

isolatesurl<-dumm.db$href[grep(paste(results$rMLST_taxon,"isolates"),dumm.db$description)]
sequrl<-dumm.db$href[grep(paste(results$rMLST_taxon,"sequence"),dumm.db$description)]

write(paste(sequrl,"/schemes/1/sequence",sep = ""), stdout())
