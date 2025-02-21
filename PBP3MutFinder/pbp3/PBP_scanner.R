
library(seqinr)

#system("/home/nacho/anaconda3/bin/makeblastdb -in /media/nacho/Data/DockerImages/TOP/PBP3MutFinder/pbp3.fasta -out /media/nacho/Data/DockerImages/TOP/PBP3MutFinder/pbp3db -dbtype nucl")

fastas<-list.files( pattern = "*\\.fa*" , full.names = TRUE)
ref<-read.fasta("/home/docker/pbp3/pbp3.fasta")

mutlist<-read.csv("/home/docker/pbp3/Mutationlist.csv", sep = ";")
reslist<-list()

for (i in 1:length(fastas)) {
  system(paste("/home/docker/miniconda3/bin/blastn -query ", fastas[i],
               " -outfmt \'6 qseqid sseqid qstart qend sstart send slen qlen length mismatch evalue\' -db  /home/docker/pbp3/pbp3db -out /home/docker/pbp3/blast_results.tsv -num_threads 4 -perc_identity 90",
               sep = ""))
  if(file.exists("/home/docker/pbp3/blast_results.tsv")){
  df<-read.csv("/home/docker/pbp3/blast_results.tsv", sep = "\t", header = FALSE)
  colnames(df)<- c("qseqid", "sseqid", "qstart", "qend", "sstart", "send" , "slen" , "qlen" ,"length", "mismatch" ,"evalue")
  sq<-read.fasta(fastas[i])
  contig<-sq[[which(names(sq)==df$qseqid[1])]]
  
  if(df$sstart>df$send){
    genet<-contig[c(df$qend:df$qstart)]
    gene<-genet
    gene[which(tolower(genet) =="a")]<-"t"
    gene[which(tolower(genet) =="c")]<-"g"
    gene[which(tolower(genet) =="t")]<-"a"
    gene[which(tolower(genet) =="g")]<-"c"
    
  }else{
    gene<-contig[c(df$qstart:df$qend)]  
  }
  
  prot<-translate(gene)
  prot.ref<-translate(unlist(ref))
  if(length(which(prot.ref!=prot))>0){
    
    muts<-paste(paste("PBP3:",prot.ref[which(prot.ref!=prot)],which(prot.ref!=prot) ,prot[which(prot.ref!=prot)],sep = ""), collapse = "/")
    clas<-match(paste("PBP3:",prot.ref[which(prot.ref!=prot)],which(prot.ref!=prot) ,prot[which(prot.ref!=prot)],sep = ""), mutlist$Mutation)
    clas<-unique(clas)
    if(length(which(is.na(clas)))>0) clas<-clas[-which(is.na(clas))]
    clas<-paste(unique(mutlist$Class[clas]), collapse = "+")
    
    
  }else{
    muts<-"None"
    clas<-NA
  }
  reslist<-c(reslist,list(as.data.frame(t(c(gsub("_.*","",gsub(".*/","",fastas[i])), muts, clas))) ))
  file.remove("/home/docker/pbp3/blast_results.tsv")
  }
}

results<-do.call(rbind,reslist )
colnames(results)<-c("Sample", "PBP3Mutations", "Class")

writexl::write_xlsx(results, "PBP3Mutations.xlsx")
write.csv(results, "PBP3Mutations.csv", row.names = FALSE)