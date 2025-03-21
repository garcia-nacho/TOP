library(jsonlite)

if(file.exists("profiles_cgmlst.tsv")){
  cg_profile<-read.csv("profiles_cgmlst.tsv",sep = "\t")
  cg_profile_date<-paste(gsub("-","",Sys.Date()), "_PASTEUR.FR",sep = "")
  
}else{
  cg_profile<-read.csv("/home/docker/profiles_scheme_cgmlst.tsv",sep = "\t")
  cg_profile_date<-gsub("20250320_LOCAL")
}


if(file.exists("profiles_mlst.tsv")){
  mlst_profile<-read.csv("profiles_cgmlst.tsv",sep = "\t")
  mlst_profile_date<-paste(gsub("-","",Sys.Date()), "_PASTEUR.FR",sep = "")
  
}else{
  mlst_profile<-read.csv("/home/docker/profiles_scheme_mlst.tsv",sep = "\t")
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
    
    if( is.null(df$status) ){
      output$MLST_BPE_Warning<-NA
    }else{
      output$MLST_BPE_Warning<- paste(unlist(df), collapse = ". ")
    }
    
    output$MLST.Scheme<-paste(sch[order(sch)],collapse = " | ")
    if(is.na(sch[1])) output$MLST.Scheme<-NA
    
    output$Sample<-gsub("_mlst.json","",inputs[f])
    if(length(which(colnames(output)=="ClonalComplex"))==0) output$ClonalComplex<-NA
    if(!exists("mlst_out")){
      mlst_out<-output
    }else{
      mlst_out<-rbind(mlst_out,output)
    }
    if(length(df)==0)mlst_out$MLST_BPE_Warning<-"Server error"
  }else{

    mlst_out<-as.data.frame(NA)  
    colnames(mlst_out)<-"MLST.Type"
    mlst_out$ClonalComplex<-NA
    mlst_out$MLST_BPE_Warning<-"Server error"
    mlst_out$Sample<-gsub("_mlst.json","",inputs[f])


  }


}



#AMR

inputs<-list.files(pattern = "_amrst.json")
outamr<-list()
for (f in 1:length(inputs)) {
  if(exists("df")) rm(df)
  
  try(df<-fromJSON(inputs[f]))
  
  if(exists("df")){
    
    if(!is.null(df$exact_matches)){
      
      if(!is.null(df$exact_matches$`23S_rRNA`)){
        BPE23S<-df$exact_matches$`23S_rRNA`$allele_id
        if(length(BPE23S)<3) BPE23S<-c(BPE23S, "Unknown")
        BPE23S<-paste(paste("ID:", BPE23S,sep=""),collapse = " | ")
        df$exact_matches$`23S_rRNA`<-NULL
        
      }else{
        BPE23S<-"Non detected"
      }

      sch<-vector()
      if(length(length(df$exact_matches))>0){
      for (i in 1:length(df$exact_matches)) {
        sch<-c(sch,paste(names(df$exact_matches)[i],":",df$exact_matches[[i]]$allele_id,sep = ""))
        #output$dummy<-paste(df$exact_matches[[i]]$allele_id,collapse = "/")
        #colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
      }
      }else{
        sch<-NA
      }
      sch<-paste(sch, collapse = " | ")
      
      output<- as.data.frame(t(c(BPE23S, sch)))
      colnames(output)<-c("BPE_23S", "BPE_AMRscheme")
      
      }else{
      BPE23S<-"Non detected"  
      output<-as.data.frame(t(c(BPE23S, NA)))  
      colnames(output)<-c("BPE_23S", "BPE_AMRscheme")
    }  

  }else{
      output<-as.data.frame(t(c("DBError", "DBError")))
      colnames(output)<-c("BPE_23S", "BPE_AMRscheme")
  }
  output$Sample<-gsub("_amrst.json","",inputs[f])
  outamr<-c(outamr, list(output))
}

outamr<-do.call(rbind, outamr)


#Phase

inputs<-list.files(pattern = "_phasest.json")
outphase<-list()
for (f in 1:length(inputs)) {
  if(exists("df")) rm(df)
  
  try(df<-fromJSON(inputs[f]))
  
  if(exists("df")){
    
    if(!is.null(df$exact_matches)){
      
      sch<-vector()
      if(length(length(df$exact_matches))>0){
        for (i in 1:length(df$exact_matches)) {
          sch<-c(sch,paste(names(df$exact_matches)[i],":",df$exact_matches[[i]]$allele_id,sep = ""))
          #output$dummy<-paste(df$exact_matches[[i]]$allele_id,collapse = "/")
          #colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
        }
      }else{
        sch<-NA
      }
      sch<-paste(sch, collapse = " | ")
      
      output<- as.data.frame(c(sch))
      colnames(output)<-"BPE_PhaseScheme"
      
    }else{

      output<-as.data.frame(t(NA))  
      colnames(output)<-"BPE_PhaseScheme"
    }  
    
  }else{
    output<-as.data.frame(c("DBError"))
    colnames(output)<-"BPE_PhaseScheme"
  }
  output$Sample<-gsub("_phasest.json","",inputs[f])
  outphase<-c(outphase, list(output))
}

outphase<-do.call(rbind, outphase)


#Toxin

inputs<-list.files(pattern = "_toxinst.json")
outtoxin<-list()
for (f in 1:length(inputs)) {
  if(exists("df")) rm(df)
  
  try(df<-fromJSON(inputs[f]))
  
  if(exists("df")){
    
    if(!is.null(df$exact_matches)){
      
      sch<-vector()
      if(length(length(df$exact_matches))>0){
        for (i in 1:length(df$exact_matches)) {
          sch<-c(sch,paste(names(df$exact_matches)[i],":",df$exact_matches[[i]]$allele_id,sep = ""))
          #output$dummy<-paste(df$exact_matches[[i]]$allele_id,collapse = "/")
          #colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
        }
      }else{
        sch<-NA
      }
      sch<-paste(sch, collapse = " | ")
      
      output<- as.data.frame(c(sch))
      colnames(output)<-"BPE_ToxinScheme"
      
    }else{
      
      output<-as.data.frame(t(NA))  
      colnames(output)<-"BPE_ToxinScheme"
    }  
    
  }else{
    output<-as.data.frame(c("DBError"))
    colnames(output)<-"BPE_ToxinScheme"
  }
  output$Sample<-gsub("_toxinst.json","",inputs[f])
  outtoxin<-c(outtoxin, list(output))
}

outtoxin<-do.call(rbind, outtoxin)


#Vaccine profile


inputs<-list.files(pattern = "_bpagst.json")
outvac<-list()
for (f in 1:length(inputs)) {
  if(exists("df")) rm(df)
  
  try(df<-fromJSON(inputs[f]))
  
  if(exists("df")){
    
    if(!is.null(df$fields$BPagST)){
      agst<-df$fields$BPagST
    }else{
      agst<-NA
    }
  
    
    if(!is.null(df$exact_matches)){
      
      sch<-vector()
      if(length(length(df$exact_matches))>0){
        for (i in 1:length(df$exact_matches)) {
          sch<-c(sch,paste(names(df$exact_matches)[i],":",df$exact_matches[[i]]$allele_id,sep = ""))
          #output$dummy<-paste(df$exact_matches[[i]]$allele_id,collapse = "/")
          #colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
        }
      }else{
        sch<-NA
      }
      sch<-paste(sch, collapse = " | ")
      
      output<- as.data.frame(c(sch))
      colnames(output)<-"BPE_AgProfileScheme"
      output$BPE_AgProfileType<-agst
      
    }else{
      
      output<-as.data.frame(t(NA))  
      colnames(output)<-"BPE_AgProfileScheme"
      output$BPE_AgProfileType<-agst
      
    }  
    
  }else{
    output<-as.data.frame(c("DBError"))
    colnames(output)<-"BPE_AgProfileScheme"
    output$BPE_AgProfileType<-"DBError"
  }
  output$Sample<-gsub("_bpagst.json","",inputs[f])
  outvac<-c(outvac, list(output))
}

outvac<-do.call(rbind, outvac)



#TSS3


inputs<-list.files(pattern = "_tss3st.json")
outtss3<-list()
for (f in 1:length(inputs)) {
  if(exists("df")) rm(df)
  
  try(df<-fromJSON(inputs[f]))
  
  if(exists("df")){
    
    if(!is.null(df$fields$T3SStype)){
      agst<-df$fields$T3SStype
    }else{
      agst<-NA
    }
    
    
    if(!is.null(df$exact_matches)){
      
      sch<-vector()
      if(length(length(df$exact_matches))>0){
        for (i in 1:length(df$exact_matches)) {
          sch<-c(sch,paste(names(df$exact_matches)[i],":",df$exact_matches[[i]]$allele_id,sep = ""))
          #output$dummy<-paste(df$exact_matches[[i]]$allele_id,collapse = "/")
          #colnames(output)[which(colnames(output)=="dummy")]<-names(df$exact_matches)[i]
        }
      }else{
        sch<-NA
      }
      sch<-paste(sch, collapse = " | ")
      
      output<- as.data.frame(c(sch))
      colnames(output)<-"BPE_T3SSScheme"
      output$BPE_T3SSType<-agst
      
    }else{
      
      output<-as.data.frame(t(NA))  
      colnames(output)<-"BPE_T3SSScheme"
      output$BPE_T3SSType<-agst
      
    }  
    
  }else{
    output<-as.data.frame(c("DBError"))
    colnames(output)<-"BPE_T3SSScheme"
    output$BPE_T3SSType<-"DBError"
  }
  output$Sample<-gsub("_tss3st.json","",inputs[f])
  outtss3<-c(outtss3, list(output))
}

outtss3<-do.call(rbind, outtss3)



#cgMLST

inputs<-list.files(pattern = "_cgmlst.json")

for (f in 1:length(inputs)) {

  if(exists("df")) rm(df)

  try(df<-fromJSON(inputs[f]))

  if(exists("df")){
  
    if(is.null(df$exact_matches)){
      cgdf<-as.data.frame(NA)
    }else{
      cgdf<-as.data.frame(names(df$exact_matches))
    }
    
    colnames(cgdf)<-"Allele"
    cgdf$ID<-NA
    if(!is.null(df$exact_matches)){
      for (i in 1:nrow(cgdf)) {
        cgdf$ID[i]<-df$exact_matches[[i]]$allele_id
      }
    }
    all.count<-ncol(cg_profile)-1
    cg_profile.temp<-cg_profile[,which(colnames(cg_profile) %in% cgdf$Allele)]
    
    distance<-vector()
    cg_profile.temp<-cg_profile.temp[,order(colnames(cg_profile.temp))]
    cgdf<-cgdf[order(cgdf$Allele),]
    pb<-txtProgressBar(max = nrow(cg_profile.temp))
    for (r in 1:nrow(cg_profile.temp)) {
      setTxtProgressBar(pb,r)
      dumcount<-ncol(cg_profile.temp) - length(which(cg_profile.temp[r,]=="N")) - length(which(is.na(cg_profile.temp[r,])))
      distance<- c(distance, length(which(cg_profile.temp[r,]==cgdf$ID))/dumcount)       
    }
    
    index<-which(distance==max(distance))
    cg<-as.data.frame( paste(paste("cgST", cg_profile$cgST[index], sep = ""),collapse = "|") )
    colnames(cg)<-"cgMLST"
    cg$Score<-max(distance)
    cg$Sample<-gsub( "_cgmlst.json","", inputs[f])
    cg$AllelesFound<-paste(ncol(cg_profile.temp), " of ", (ncol(cg_profile)-1))
    cg$cgMLSTDate<-cg_profile_date
    cg$cgMLST_BPE_Warning<-NA
    
    if(!exists("cg_out")){
      cg_out<-cg
    }else{
      cg_out<-rbind(cg_out,cg)
    }
    if(length(df)==0)cg_out$cgMLST_BPE_Warning<-"Server error"
    if(length(df)==0)cg_out$cgMLST<-NA
    if(length(df)==0)cg_out$Score<-NA
  }else{
    cg_out<-as.data.frame(NA)
    colnames(cg_out)<-"cgMLST"
    cg_out$Score<-NA
    cg_out$Sample<-gsub( "_cgmlst.json","", inputs[f])
    cg_out$AllelesFound<-NA
    cg_out$cgMLSTDate<-NA
    cg_out$cgMLST_BPE_Warning<-"Server Not responding"

  }
}

mlst_out<-merge(mlst_out, outamr,  by="Sample", all.x = TRUE, all.y = TRUE)
mlst_out<-merge(mlst_out, outvac,  by="Sample", all.x = TRUE, all.y = TRUE)
mlst_out<-merge(mlst_out, outtoxin,  by="Sample", all.x = TRUE, all.y = TRUE)
mlst_out<-merge(mlst_out, outtss3,  by="Sample", all.x = TRUE, all.y = TRUE)


finalout<-merge(mlst_out, cg_out, by="Sample", all.x = TRUE, all.y = TRUE)


write.csv(finalout, "BPE_MLST.csv", row.names = FALSE)
