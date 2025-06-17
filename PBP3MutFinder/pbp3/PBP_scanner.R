
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
  
  if(file.info("/home/docker/pbp3/blast_results.tsv")$size > 0){
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
    muts_nc<-paste("PBP3:",prot.ref[which(prot.ref!=prot)],which(prot.ref!=prot) ,prot[which(prot.ref!=prot)],sep = "")
    clas<-match(paste("PBP3:",prot.ref[which(prot.ref!=prot)],which(prot.ref!=prot) ,prot[which(prot.ref!=prot)],sep = ""), mutlist$Mutation)
    claslong<-mutlist$Class[clas]
    if(length(which(claslong=="Stage1"))>0){
      S1<-paste(muts_nc[which(claslong=="Stage1")]  , collapse = "/")  
      
    }else{
      S1<-"NonDetected"
    }
    
    if(length(which(claslong=="Stage2"))>0){
      S2<-paste(muts_nc[which(claslong=="Stage2")]  , collapse = "/") 
    }else{
      S2<-"NonDetected"
    }
    
    if(length(which(claslong=="Stage3"))>0){
      S3<-paste(muts_nc[which(claslong=="Stage3")]  , collapse = "/") 
    }else{
      S3<-"NonDetected"
    }
    
    if(length(which(claslong %in% c("Stage1","Stage2","Stage3")))>0){
      others<-paste(muts_nc[-which(claslong %in% c("Stage1","Stage2","Stage3"))]  , collapse = "/") 
    }else{
      others<-muts
    }
    
    
      
    clas<-unique(clas)
    if(length(which(is.na(clas)))>0) clas<-clas[-which(is.na(clas))]
    clas<-paste(unique(mutlist$Class[clas]), collapse = "+")
    
    
  }else{
    muts<-"None"
    clas<-"None"
    S1<-"NonDetected"
    S2<-"NonDetected"
    S3<-"NonDetected"
    others<-"None"
  }
  file.remove("/home/docker/pbp3/blast_results.tsv")
  reslist<-c(reslist,list(as.data.frame(t(c(gsub("_.*","",gsub(".*/","",fastas[i])), muts, clas, S1, S2, S3, others))) ))
  }else{
    
    reslist<-c(reslist,list(as.data.frame(t(c(gsub("_.*","",gsub(".*/","",fastas[i])), "NonDetectedPBP3", "NonDetectedPBP3", "NonDetectedPBP3", "NonDetectedPBP3", "NonDetectedPBP3", "NonDetectedPBP3"))) ))
    
  }
  
 
  
}

results<-do.call(rbind,reslist )
colnames(results)<-c("Sample", "PBP3Mutations", "Class","HinfPBP3_Stage1","HinfPBP3_Stage2","HinfPBP3_Stage3","HinfPBP3_Other")

writexl::write_xlsx(results, "PBP3Mutations.xlsx")
write.csv(results, "PBP3Mutations.csv", row.names = FALSE)