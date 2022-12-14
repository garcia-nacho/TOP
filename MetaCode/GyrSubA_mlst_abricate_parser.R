library(seqinr)

dbname<-"Hinf_GyrSubA"

df<-read.csv("/media/nacho/Data/DockerImages/TOP/Abricate/Hi_GyraseSubA.txt",sep = "\t")

seqs<-list()

for (i in 1:nrow(df)) {

  seqs[[i]]<-df$sequence[i]
  
  names(seqs)[i]<-paste(dbname,"~~~" ,paste(df$locus[i],"_ID",df$allele_id[i],sep = ""),"~~~" ,paste(df$locus[i],"_ID",df$allele_id[i],sep = "") ,"~~ND",
                       sep = "")
  names(seqs)<-gsub("~~$","~~ND", names(seqs))
}

write.fasta(seqs, "/media/nacho/Data/DockerImages/TOP/Abricate/Hi_GyrSubA.fasta", names = names(seqs))
