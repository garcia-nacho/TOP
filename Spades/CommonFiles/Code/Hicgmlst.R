library(jsonlite)


cg_profile<-read.csv("profiles_cgmlst.tsv",sep = "\t")
cg_profile_date<-paste("15052025")

#cgMLST

inputs<-list.files(pattern = "_Hicgmlst.json")

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
    colnames(cg)<-"HicgMLST"
    cg$HicgMLST_Score<-max(distance)
    cg$Sample<-gsub( "_Hicgmlst.json","", inputs[f])
    cg$HicgMLST_AllelesFound<-paste(ncol(cg_profile.temp), " of ", (ncol(cg_profile)-1))
    cg$HicgMLST_Warning<-NA
    
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
    colnames(cg_out)<-"HicgMLST"
    cg_out$HicgMLST_Score<-NA
    cg_out$Sample<-gsub( "_Hicgmlst.json","", inputs[f])
    cg_out$HicgMLST_AllelesFound<-NA
    cg_out$HicgMLST_Warning<-"Server Not responding"

  }
}


write.csv(cg_out, "Hicgmlst.csv", row.names = FALSE)
