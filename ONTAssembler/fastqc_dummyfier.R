#Contig stats

nanoplot<-list.files(pattern = "NanoStats.txt")


for (i in 1:length(nanoplot)) {
  
  dir.create(paste(gsub("_NanoStats.txt","",nanoplot[i]), "_R1_fastqc", sep = ""))
  dir.create(paste(gsub("_NanoStats.txt","",nanoplot[i]), "_R2_fastqc", sep = ""))
  dir.create(paste(gsub("_NanoStats.txt","",nanoplot[i]), "_1P_fastqc", sep = ""))
  dir.create(paste(gsub("_NanoStats.txt","",nanoplot[i]), "_2P_fastqc", sep = ""))
  
  dum<-readLines(nanoplot[i])
  q10<-dum[grep(">Q10",dum)]
  q15<-dum[grep(">Q15",dum)]
  q20<-dum[grep(">Q20",dum)]
  q25<-dum[grep(">Q25",dum)]
  q30<-dum[grep(">Q30",dum)]
  
  dum2<-as.data.frame(rbind(q10, q15,q20,q25,q30))
  colnames(dum2)<-"Quality"
  dum2$Count<-gsub(".*:", "",dum2$Quality)
  dum2$Count<-as.numeric(gsub(" .*", "",dum2$Count))
  
  dum2$Quality<-as.numeric(c("10","15","20","25","30"))

  
  pad<-dum[grep("Number of reads:",dum)]
  pad<-as.data.frame(as.numeric(gsub(",","",gsub(".*:","",pad))))
  colnames(pad)<-"Count"
  pad$Quality<-1
  pad$Count<-pad$Count-max(dum2$Count)

    
  dum2<-rbind(dum2, pad)
  dum2<-dum2[order(dum2$Quality,decreasing = FALSE),]
  colnames(dum2)[1]<-"#Quality"
  
  write.table(dum2, paste(gsub("_NanoStats.txt","",nanoplot[i]), "_R1_fastqc", "/fastqc_data.txt",sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(dum2, paste(gsub("_NanoStats.txt","",nanoplot[i]), "_R2_fastqc", "/fastqc_data.txt",sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(dum2, paste(gsub("_NanoStats.txt","",nanoplot[i]), "_1P_fastqc", "/fastqc_data.txt",sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(dum2, paste(gsub("_NanoStats.txt","",nanoplot[i]), "_2P_fastqc", "/fastqc_data.txt",sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE)

  
  zip(paste(gsub("_NanoStats.txt","",nanoplot[i]), "_R1_fastqc.zip",sep = "") ,
      paste(gsub("_NanoStats.txt","",nanoplot[i]), "_R1_fastqc", "/fastqc_data.txt",sep = ""))
  
  zip(paste(gsub("_NanoStats.txt","",nanoplot[i]), "_R2_fastqc.zip",sep = "") ,
      paste(gsub("_NanoStats.txt","",nanoplot[i]), "_R2_fastqc", "/fastqc_data.txt",sep = ""))
  
  zip(paste(gsub("_NanoStats.txt","",nanoplot[i]), "_1P_fastqc.zip",sep = "") ,
      paste(gsub("_NanoStats.txt","",nanoplot[i]), "_1P_fastqc", "/fastqc_data.txt",sep = ""))
  
  zip(paste(gsub("_NanoStats.txt","",nanoplot[i]), "_2P_fastqc.zip",sep = "") ,
      paste(gsub("_NanoStats.txt","",nanoplot[i]), "_2P_fastqc", "/fastqc_data.txt",sep = ""))
  
}

system(paste("rm -rf ", paste(gsub("_NanoStats.txt","",nanoplot[i]), "_R1_fastqc", sep = ""), sep = ""))
system(paste("rm -rf ", paste(gsub("_NanoStats.txt","",nanoplot[i]), "_R2_fastqc", sep = ""), sep = ""))
system(paste("rm -rf ", paste(gsub("_NanoStats.txt","",nanoplot[i]), "_1P_fastqc", sep = ""), sep = ""))
system(paste("rm -rf ", paste(gsub("_NanoStats.txt","",nanoplot[i]), "_2P_fastqc", sep = ""), sep = ""))

