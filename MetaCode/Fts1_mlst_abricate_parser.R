library(seqinr)

dbname<-"Hinf_FtsI"

df<-read.csv("/media/nacho/Data/DockerImages/TOP/Abricate/Hi_Fts1.txt",sep = "\t")

seqs<-list()

for (i in 1:nrow(df)) {

  seqs[[i]]<-df$sequence[i]
  colnames(df)[c(6:17)]
  mut.to.report<-df[i,c(6:17)]
  mut.to.report<-mut.to.report[-which(mut.to.report[1,]=="")]
  
  names(seqs)[i]<-paste(dbname,"~~~" ,paste(df$locus[i],"_ID",df$allele_id[i],sep = ""),"~~~" ,paste(df$locus[i],"_ID",df$allele_id[i],sep = "") ,"~~",
                       paste(mut.to.report[1,],collapse = "/"),
                       sep = "")
  names(seqs)<-gsub("~~$","~~ND", names(seqs))
}

write.fasta(seqs, "/media/nacho/Data/DockerImages/TOP/Abricate/Hi_Fts1.fasta", names = names(seqs))
