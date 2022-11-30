library(seqinr)


size.co<-500
cov.co<-2

genomes<-list.files( pattern = "_contigs.fasta", full.names = TRUE)

summaries<-data.frame(matrix(NA, ncol=9, nrow=length(genomes)))
colnames(summaries)<-c("Sample", "N50", "L50","RawContigCount","CleanContigCount","Coverage","AverageDepth","SD.Depth","Normalized.Depth")
dummyfa<-read.fasta(genomes[1])
summaries$Sample<-gsub("_.*","",gsub("_contigs.fasta","",gsub(".*/","",genomes[1])))

summaries$RawContigCount<-length(dummyfa)
coverage<- as.numeric(gsub(".*_cov_","",names(dummyfa)))
sizes<-as.numeric(unlist(lapply(dummyfa, length)))
dummyfa<-dummyfa[-unique(c(which(sizes<size.co), which(coverage<cov.co)))]

summaries$CleanContigCount<-length(dummyfa)

summaries$Coverage<-sum(unlist(lapply(dummyfa, length)))
summaries$AverageDepth<-mean(coverage)
summaries$SD.Depth<-sd(coverage)
summaries$Normalized.Depth<- sum(sizes*coverage)/sum(sizes)



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

names(dummyfa)<-paste(summaries$Sample[1],gsub("_length.*","",names(dummyfa)))
write.fasta(dummyfa, names = names(dummyfa), file.out ="clean_contigs.fasta")
write.csv(summaries, "clean_contigs.stats.csv", row.names = FALSE)
