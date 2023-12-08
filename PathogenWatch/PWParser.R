library("RcppTOML")
library(seqinr)

df<-parseTOML("/media/nacho/Data/temp/pathogenwatch/485.toml")
inputfasta<-read.fasta("/media/nacho/Data/temp/toptest/pathogenwatch/23GB001751.fasta")

amr_genes<-list()

for (i in 1:length(df$genes)) {
amr_genes[[i]]<- df$genes[[i]]$sequence  
names(amr_genes)[i]<-df$genes[[i]]$name
}

write.fasta(amr_genes, "/media/nacho/Data/temp/pathogenwatch/485.fasta", names = names(amr_genes))
subject<-read.fasta("/media/nacho/Data/temp/pathogenwatch/485.fasta")
#Blast
print("Running blast")


blastresults<-read.csv("/media/nacho/Data/temp/toptest/pathogenwatch/test.txt",sep = "\t", header = FALSE)

colnames(blastresults)<-c("qseqid", "sseqid"  ,  "pident" ,  "length_aln" , "mismatch" , "gapopen",  "qstart",  "qend" ,"sstart","send","evalue","bitscore",
                          "queryseq", "subjectseq")

genes<-unique(blastresults$sseqid)

hit.index<-vector()

for (i in 1:length(genes)) {
  hit.index<-c(hit.index,  which(blastresults$sseqid==genes[i] & blastresults$evalue==min(blastresults$evalue[which(blastresults$sseqid==genes[i])]) )) 
  
}

#Filter by Gene & PID


blasthits<-blastresults[hit.index,]
blasthits<-blasthits[which(blasthits$length_aln>100),]
blasthits$variants<-NA
blasthits$MutationsDNA<-NA
blasthits$MutationsProt<-NA

for (i in 1:nrow(blasthits)) {
  query.dna <- unlist(strsplit(blasthits$queryseq[i],"") )
  sub.dna<-unlist(strsplit(blasthits$subjectseq[i],"") )
  
  if(blasthits$sstart[i]>blasthits$send[i]){
    sub.dna<-toupper(rev(comp(sub.dna)))
    query.dna<-toupper(rev(comp(query.dna)))
    
    
  }
  if(paste(sub.dna[c(1:3)], collapse = "")=="ATG"){
    sub.pro<-translate(sub.dna)
    query.pro<-translate(query.dna)
    if(length(which(sub.pro != query.pro))>0){
      blasthits$MutationsProt[i]<-paste(paste(sub.pro[which(sub.pro != query.pro)], which(sub.pro != query.pro),
                                              query.pro[which(sub.pro != query.pro)],sep = ""),collapse = "/")
    } 
  
  }
  
  if(length(which(sub.dna != query.dna))>0){
    blasthits$MutationsDNA[i]<-paste(paste(sub.dna[which(sub.dna != query.dna)], which(sub.dna != query.dna), query.dna[which(sub.dna != query.dna)],sep = ""),collapse = "/")
  } 
  
}

#Get IDs


#Mechanisms
name.vec<-vector()
for (i in 1:length(df$mechanisms)) {
  name.vec<-c(name.vec, names(unlist(df$mechanisms[[i]])))
}

mechdf<- as.data.frame(matrix(NA, nrow =length(df$mechanisms), ncol =  2))


colnames(mechdf)<-c("phenotype", "amr_factors")

for (i in 1:nrow(mechdf)) {
  dum<-unlist(df$mechanisms[[i]],"")
  phen<- vector()
  for (tem in 1:length(df$mechanisms[[i]]$phenotypes)) {
    phen<-c(phen, paste(df$mechanisms[[i]]$phenotypes[[tem]]$effect,  df$mechanisms[[i]]$phenotypes[[tem]]$profile))
  }
  
  
  mechdf$phenotype[i]<-paste(phen, collapse = " / ")
  
  memb<- vector()
  for (mb in 1:length(df$mechanisms[[i]]$members)) {
    if(length(df$mechanisms[[i]]$members[[mb]])>1){
      memb<-c(memb, paste(df$mechanisms[[i]]$members[[mb]]$gene,df$mechanisms[[i]]$members[[mb]]$variants,sep = ":"  ))  
    }else{
      memb<-c(memb, df$mechanisms[[i]]$members[[mb]])
    }
      
    
  }
  mechdf$amr_factors[i]<-paste(memb,collapse = " / ")
}


