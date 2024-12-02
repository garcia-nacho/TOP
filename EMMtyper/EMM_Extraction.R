library(seqinr)
#setwd("/media/nacho/Data/temp/toptest/TopValidation/emm/21nov2024/")
fastas.list<-list.files(pattern = ".fasta")

#blastn -out 2265046-GA-24GB002298_merged_clean_contigs_emm_results_blast.tab -outfmt '6 qseqid sseqid qstart qend sstart send slen qlen length mismatch evalue' -query /Data/2265046-GA-24GB002298_merged_clean_contigs.fasta -db /emmdb07022024/emmdb.tfa


primersq<-"TATTCGCTTAGAAAATTAA"
primer<-unlist(strsplit(primersq, ""))

for (i in 1:length(fastas.list)) {
  system(paste("blastn -out ",
               gsub(".fasta", "_emm_results_blast.tsv",fastas.list[i]),
               " -outfmt '6 qseqid sseqid qstart qend sstart send slen qlen length mismatch evalue' -query ",
               fastas.list[i] ," -db /emmdb07022024/emmdb.tfa",sep = ""))
}

results<-list.files(pattern = "_emm_results_blast.tsv")

sq<-list()
for (i in 1:length(results)) {
  res<-read.csv(results[i], sep = "\t", header = FALSE)
  res<-res[which(res$V11==min(res$V11))[1],]
  fastas<-read.fasta(fastas.list[grep(gsub("_.*","",res$V1[1]),fastas.list)])
  fastas<-fastas[which(names(fastas)==res$V1[1])]
  
  if(res$V5<res$V6){
    #Direct
    start<-max(min(res$V3,res$V4)-500, 1)
    end<- min(max(res$V3,res$V4)+500, res$V8)
    tmp.sq<-unlist(fastas)[c(start:end)]
    
    score<-vector()
    for (scan in 1:(length(tmp.sq)-length(primer))) {
      score <- c(score,  length(which(toupper(as.character(tmp.sq[scan:(scan+length(primer)-1)])) ==
                                        toupper(primer))  )) 
    }
    score<-score/length(primer) 
    
    if(length(which(score==1))){
      tmp.sq<-tmp.sq[which(score==1): min(length(tmp.sq), which(score==1)+500)]
    }
    
    
    
  }else{
    #Reverse
    start<- max(min(res$V3,res$V4)-500, 1)
    end<- min(max(res$V3,res$V4)+500, res$V8)
    tmp.sq<-unlist(fastas)[c(start:end)]
    tmp.sq<-comp(tmp.sq)[c(length(tmp.sq):1)]
    
    score<-vector()
    for (scan in 1:(length(tmp.sq)-length(primer))) {
      score <- c(score,  length(which(toupper(as.character(tmp.sq[scan:(scan+length(primer)-1)])) ==
                                        toupper(primer))  )) 
    }
    score<-score/length(primer) 
    
    if(length(which(score==1))){
      tmp.sq<-tmp.sq[which(score==1): min(length(tmp.sq), which(score==1)+500)]
    }
    
  }
  
  
  sq<-c(sq,list(tmp.sq))
  names(sq)[length(sq)]<-gsub("_.*","",res$V1[1])
}

write.fasta(sq,"EMM_seqs_extended.fa", names = paste("EMM_",names(sq),sep = ""))
