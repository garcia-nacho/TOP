library(readxl)

tsvs<-list.files(pattern = ".tsv$", recursive = TRUE)
csvs<-list.files(pattern = ".csv$",recursive = TRUE)
xls<-list.files(pattern = "xlsx$",recursive = TRUE)

files<-list.files(pattern = "*.fastq.gz", recursive = TRUE)
folders<-unique(gsub("/.*","",files))

map<-as.data.frame(folders)
colnames(map)<-"folder"
map$newname<-folders
map$clean<-NA

continue<-TRUE
if(length(tsvs)==1){
  
  dum<-read.csv(xls, sep = "\t")
  if(length(which(colnames(dum) %in% c("Barcode","SequenceID") ))==2){
    if(length(which(duplicated(dum$Barcode)))>0) dum<-dum[-which(duplicated(dum$Barcode)),]
    if(length(which(duplicated(dum$ID)))>0) dum<-dum[-which(duplicated(dum$SI)),]
    map$clean<-dum$ID [match(map$folder,dum$Barcode) ] 
  }
  continue<-FALSE
}

if(length(csvs)==1 & continue){
  dum<-read.csv(csvs)
  if(length(which(colnames(dum) %in% c("Barcode","SequenceID") ))==2){
    if(length(which(duplicated(dum$Barcode)))>0) dum<-dum[-which(duplicated(dum$Barcode)),]
    if(length(which(duplicated(dum$SequenceID)))>0) dum<-dum[-which(duplicated(dum$SequenceID)),]
    map$clean<-dum$SequenceID [match(map$folder,dum$Barcode) ] 
  }
  continue<-FALSE
}


if(length(xls)==1 & continue){
  dum<-read_xlsx(xls)
  if(length(which(colnames(dum) %in% c("Barcode","SequenceID") ))==2){
    if(length(which(duplicated(dum$Barcode)))>0) dum<-dum[-which(duplicated(dum$Barcode)),]
    if(length(which(duplicated(dum$SequenceID)))>0) dum<-dum[-which(duplicated(dum$SequenceID)),]
    map$clean<-dum$ID [match(map$folder,dum$Barcode) ] 
  }
  continue<-FALSE
}

if(length(which(is.na(map$clean)))>0){
  map$clean[which(is.na(map$clean))]<-map$newname[which(is.na(map$clean))]
}

dir.create("merged")
for (i in 1:nrow(map)) {
  
  system(paste("cat ", map$folder[i],"/*.fastq.gz > merged/",map$clean[i],"_merged.fastq.gz", sep = ""))

  if(file.size(paste("merged/",map$clean[i],"_merged.fastq.gz", sep = ""))/10^6 < 100 ){
    
    file.remove(paste("merged/",map$clean[i],"_merged.fastq.gz", sep = ""))
  }
  if(map$clean[i] == "unclassified" ){
    file.remove(paste("merged/",map$clean[i],"_merged.fastq.gz", sep = ""))
  }
    
  
}



