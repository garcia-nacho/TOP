library(writexl)

results<-list.files( full.names = TRUE, recursive = TRUE,
                    pattern = "SeqSero_result.tsv")

results<-results[grep("kmer", results)]
try(rm(out))
for (i in 1:length(results)) {
  dum<-read.csv(results[i],sep = "\t")
  if(!exists("out")){
    out<-dum
  }else{
    out<-rbind(out,dum)
  }
}

out$SeroSeqMode<- "kmer"

write_xlsx(out, "Results.SeqSero2_kmer.xlsx")


results<-list.files( full.names = TRUE, recursive = TRUE,
                     pattern = "SeqSero_result.tsv")

results<-results[grep("allele", results)]
try(rm(out))
for (i in 1:length(results)) {
  dum<-read.csv(results[i],sep = "\t")
  if(!exists("out")){
    out<-dum
  }else{
    out<-rbind(out,dum)
  }
}

out$SeroSeqMode<- "allele"

write_xlsx(out, "Results.SeqSero2_allele.xlsx")

