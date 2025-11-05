library(jsonlite)

if(file.exists("profiles_mlst.tsv")){
  mlst_profile<-read.csv("profiles_mlst.tsv",sep = "\t")
  mlst_profile_date<-paste(gsub("-","",Sys.Date()), "_PASTEUR.FR",sep = "")
  
}else{
  mlst_profile<-read.csv("/home/docker/diphtoscan/data/profiles_mlst.tsv",sep = "\t")
  mlst_profile_date<-gsub("20241127_LOCAL")
}



inputs<-list.files(pattern = "_mlst.json")

for (f in 1:length(inputs)) {
  if(exists("df")) rm(df)
  
  try(df<-fromJSON(inputs[f]))
  
  if(exists("df")){
    
    if(!is.null(df$fields$ST)){
      output<-as.data.frame(df$fields$ST )
      colnames(output)<-"MLST.Type"
      sch<-vector()
      
      for (i in 1:length(df$exact_matches)) {
        sch<-c(sch,paste(names(df$exact_matches)[i],":",df$exact_matches[[i]]$allele_id,sep = ""))
        #output$dummy<-paste(df$exact_matches[[i]]$allele_id,collapse = "/")
        #colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
      }
      
    }else{
      output<-as.data.frame(NA)  
      colnames(output)<-"MLST.Type"
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
    
    if( is.null(df$status) ){
      output$MLST_Diphto_Warning<-NA
    }else{
      output$MLST_Diphto_Warning<- paste(unlist(df), collapse = ". ")
    }
    
    output$MLST.Scheme<-paste(sch[order(sch)],collapse = " | ")
    if(is.na(sch[1])) output$MLST.Scheme<-NA
    
    output$Sample<-gsub("_mlst.json","",inputs[f])
    if(!exists("mlst_out")){
      mlst_out<-output
    }else{
      mlst_out<-rbind(mlst_out,output)
    }
    if(length(df)==0)mlst_out$MLST_Diphto_Warning<-"Server error"
  }else{
    
    mlst_out<-as.data.frame(NA)  
    colnames(mlst_out)<-"MLST.Type"
    mlst_out$MLST_Diphto_Warning<-"Server error"
    mlst_out$Sample<-gsub("_mlst.json","",inputs[f])
    
    
  }
  
}

inputs<-list.files(pattern = "_pbp.json")
outamr<-list()

for (f in 1:length(inputs)) {
  if(exists("df")) rm(df)
  
  try(df<-fromJSON(inputs[f]))
  
  if(exists("df")){
    
    if(!is.null(df$exact_matches)){
      
      
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
      sch<-paste(sch[order(sch)], collapse = " | ")
      
      output<- as.data.frame(sch)
      
      colnames(output)<-c("PBP_Scheme")
      
    }else{
      output<-as.data.frame( NA)  
      colnames(output)<-c("PBP_Scheme")
    }  
    
  }else{
    output<-as.data.frame(c("DBError"))
    colnames(output)<-c("PBP_Scheme")
  }
  
  output$Sample<-gsub("_pbp.json","",inputs[f])
  outamr<-c(outamr, list(output))
}

outamr<-do.call(rbind, outamr)



# Toxin -------------------------------------------------------------------

inputs<-list.files(pattern = "_toxinst.json")
outoxin<-list()

for (f in 1:length(inputs)) {
  if(exists("df")) rm(df)
  
  try(df<-fromJSON(inputs[f]))
  
  if(exists("df")){
    
    if(!is.null(df$exact_matches)){
      
      
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
      sch<-paste(sch[order(sch)], collapse = " | ")
      
      output<- as.data.frame(sch)
      
      colnames(output)<-c("Toxin_Scheme")
      
    }else{
      output<-as.data.frame( NA)  
      colnames(output)<-c("Toxin_Scheme")
    }  
    
  }else{
    output<-as.data.frame(c("DBError"))
    colnames(output)<-c("Toxin_Scheme")
  }
  
  output$Sample<-gsub("_toxinst.json","",inputs[f])
  outoxin<-c(outoxin, list(output))
}

outoxin<-do.call(rbind, outoxin)



# spuA --------------------------------------------------------------------


inputs<-list.files(pattern = "_spuA_island.json")
outspu<-list()

for (f in 1:length(inputs)) {
  if(exists("df")) rm(df)
  
  try(df<-fromJSON(inputs[f]))
  
  if(exists("df")){
    
    if(!is.null(df$exact_matches)){
      
      
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
      sch<-paste(sch[order(sch)], collapse = " | ")
      
      output<- as.data.frame(sch)
      
      colnames(output)<-c("spuA_Island_Scheme")
      
    }else{
      output<-as.data.frame( NA)  
      colnames(output)<-c("spuA_Island_Scheme")
    }  
    
  }else{
    output<-as.data.frame(c("DBError"))
    colnames(output)<-c("spuA_Island_Scheme")
  }
  
  output$Sample<-gsub("_spuA_island.json","",inputs[f])
  outspu<-c(outspu, list(output))
}

outspu<-do.call(rbind, outspu)


mlst_out<-merge(mlst_out, outamr,  by="Sample", all.x = TRUE, all.y = TRUE)
mlst_out<-merge(mlst_out, outoxin,  by="Sample", all.x = TRUE, all.y = TRUE)
mlst_out<-merge(mlst_out, outspu,  by="Sample", all.x = TRUE, all.y = TRUE)


