library(readxl)

df<-read_xlsx("/media/nacho/Data/temp/toptest/TopValidation/Summaries_20240111.xlsx")


colnames(df)[which(colnames(df)=="MLST_Conflict")]<-"rMLST_Conflict"

df$MLST.Type[which(is.na(df$MLST.Type))]<-"Non Detected"
df$MLST.Scheme[which(is.na(df$MLST.Scheme))]<-"Non Detected"

df$Stx1<-NA
df$Stx2<-NA

if(length(which(df$StxType_Contigs=="No hit found" &  df$StxType_Mapping=="No hit found"))>0){
  df$Stx1[which(df$StxType_Contigs=="No hit found" &  df$StxType_Mapping=="No hit found")]<-"Assemblies:Non Detected | Reads:Non Detected"
  df$Stx2[which(df$StxType_Contigs=="No hit found" &  df$StxType_Mapping=="No hit found")]<-"Assemblies:Non Detected | Reads:Non Detected"
}

if(length(which(df$StxType_Contigs=="No hit found" &  df$StxType_Mapping!="No hit found" & !is.na(df$StxType_Mapping)))>0){
  df$Stx1[which(df$StxType_Contigs=="No hit found" &  df$StxType_Mapping!="No hit found" & !is.na(df$StxType_Mapping))]<-"Assemblies:Non Detected"
  df$Stx2[which(df$StxType_Contigs=="No hit found" &  df$StxType_Mapping!="No hit found" & !is.na(df$StxType_Mapping))]<-"Assemblies:Non Detected"
}

if(length(which(df$StxType_Contigs!="No hit found" &  df$StxType_Mapping=="No hit found" & !is.na(df$StxIdentity_Contigs)))>0){
  df$Stx1[which(df$StxType_Contigs!="No hit found" &  df$StxType_Mapping=="No hit found" & !is.na(df$StxIdentity_Contigs))]<-"Reads:Non Detected"
  df$Stx2[which(df$StxType_Contigs!="No hit found" &  df$StxType_Mapping=="No hit found" & !is.na(df$StxIdentity_Contigs))]<-"Reads:Non Detected"

}

for (i in 1:nrow(df)) {
  if(exists("stxvar"))rm(stxvar)
  if(exists("stxvar2"))rm(stxvar2)
  
   if(!is.na(df$StxVariant_Contigs[i])){
     stxvar<- unique(unlist(strsplit(df$StxVariant_Contigs[i], " \\| ")))
     stxvar<-gsub("\\:variant ","",stxvar)
     
     if(length(grep("stx1", stxvar))>0){
       df$Stx1[i]<-paste("Assemblies:", paste(stxvar[grep("stx1", stxvar)], collapse = ";"),sep = "")
     }
     if(length(grep("stx2", stxvar))>0){
       df$Stx2[i]<-paste("Assemblies:", paste(stxvar[grep("stx2", stxvar)], collapse = ";"),sep = "")
     }
    
   }
  
  if(!is.na(df$StxVariant_Mapping[i])){
    stxvar2<- unique(unlist(strsplit(df$StxVariant_Mapping[i], " \\| ")))
    stxvar2<-gsub("\\:variant ","",stxvar2)
    
    if(length(grep("stx1", stxvar2))>0){
      if(is.na(df$Stx1[i])){
        df$Stx1[i]<-paste("Reads:", paste(stxvar[grep("stx1", stxvar2)], collapse = ";"),sep = "")  
      }else{
        df$Stx1[i]<-paste(df$Stx1[i], paste("Reads:", paste(stxvar[grep("stx1", stxvar2)], collapse = ";"),sep = ""), sep = " | ")
      }
    }
    if(length(grep("stx2", stxvar2))>0){
      if(is.na(df$Stx2[i])){
        df$Stx2[i]<-paste("Reads:", paste(stxvar[grep("stx2", stxvar2)], collapse = ";"),sep = "")  
      }else{
        df$Stx2[i]<-paste(df$Stx2[i], paste("Reads:", paste(stxvar[grep("stx2", stxvar2)], collapse = ";"),sep = ""), sep = " | ")
      }
    }
  }
  

}


df$Stx1
df$Stx2


df$Instrument<-"NB551198"

write.table(df, "/media/nacho/Data/temp/toptest/TopValidation/TOP_results19012024.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
