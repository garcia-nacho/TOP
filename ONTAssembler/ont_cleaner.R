library(seqinr)

size.co<-500
cov.co<-2


genomes<-read.fasta(list.files(pattern = "raw_contigs.fasta"))
id<-gsub("_raw_contigs.*","",list.files(pattern = "raw_contigs.fasta"))
names(genomes)<-paste(id,"_NODE",c(1:length(names(genomes))),sep = "" )
genomes<-genomes[order(unlist(lapply(genomes,length)), decreasing = TRUE)]

summaries<-data.frame(matrix(NA, ncol=11, nrow=1))
colnames(summaries)<-c("Sample", "N50", "N90", "L50","RawContigCount","CleanContigCount","Coverage","AverageDepth","SD.Depth","Normalized.Depth","ReadDepth")
dummyfa<-genomes

summaries$Sample<-unique(gsub( "_.*","", names(genomes)))

summaries$RawContigCount<-length(dummyfa)
coverage<- NA
sizes<-as.numeric(unlist(lapply(dummyfa, length)))

if(length(which(sizes<size.co))){
  dummyfa<-dummyfa[-which(sizes<size.co)]
    
}


summaries$CleanContigCount<-length(dummyfa)

summaries$Coverage<-sum(unlist(lapply(dummyfa, length)))
summaries$AverageDepth<-NA
summaries$SD.Depth<-NA
summaries$Normalized.Depth<- NA



found<-FALSE
k<-1
while(!found){
  dummy.sum<-sum(unlist(lapply(dummyfa[1:k], length)))
  if(dummy.sum>=summaries$Coverage[1]/2) {
    summaries$N50[1]<-length(dummyfa[[k]])
    summaries$L50[1]<-k
    found<-TRUE
  }
  k<-k+1
}


found<-FALSE
k<-1
while(!found){
  dummy.sum<-sum(unlist(lapply(dummyfa[1:k], length)))
  if(dummy.sum>=summaries$Coverage[1]*0.9) {
    summaries$L90[1]<-k
    found<-TRUE
  }
  k<-k+1
}




write.fasta(dummyfa, names = names(dummyfa), file.out ="clean_contigs.fasta")
write.csv(summaries, "clean_contigs.stats.csv", row.names = FALSE)
