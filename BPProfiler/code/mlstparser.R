library(jsonlite)

if(file.exists("profiles_cgmlst.tsv")){
  cg_profile<-read.csv("profiles_cgmlst.tsv",sep = "\t")
  cg_profile_date<-paste(gsub("-","",Sys.Date()), "_PASTEUR.FR",sep = "")
  
}else{
  cg_profile<-read.csv("/home/docker/profiles_cgmlst.tsv",sep = "\t")
  cg_profile_date<-gsub("20241127_LOCAL")
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

finalout<-merge(mlst_out, cg_out, by="Sample", all.x = TRUE, all.y = TRUE)

write.csv(finalout, "BPE_MLST.csv", row.names = FALSE)
