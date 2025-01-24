snpfiles<-list.files("/media/nacho/Data/OnGoingProjects/Pertusis/fastq/SnippyResults/", pattern = "snps.csv", full.names = TRUE, recursive = TRUE)

if(exists("snpout")) rm(snpout)
for (i in 1:length(snpfiles)) {
  dum<-read.csv(snpfiles[i])
  if(nrow(dum)>0){
    dum<-dum[,c(1:6)]
    dum$Sample<-gsub(".*/","",gsub( "/snps.csv","", snpfiles[i]))
    if(!exists("snpout")){
      snpout<-dum
    }else{
      snpout<-rbind(snpout,dum)
    }
  }
}

write.csv(snpout, "/media/nacho/Data/OnGoingProjects/Pertusis/fastq/snpsBPE.csv", row.names = FALSE)
