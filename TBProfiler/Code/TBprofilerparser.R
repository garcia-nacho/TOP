library(jsonlite)

jsonfile<-list.files(pattern = "tb_profiler.json")

df<-fromJSON(jsonfile)

results<-c(df$id, df$tbprofiler_version, df$qc$pct_reads_mapped, df$qc$median_coverage,df$main_lin,df$sublin,df$drtype ,df$db_version$Date)
results<-as.data.frame(t(results))
colnames(results)<-c("Sample", "TB_profiler_Version", "TB_profiler_FractionReadsMapped","TB_profiler_MedianCoverage", "TB_profiler_MainLineage",
                     "TB_profiler_Sublineage", "TB_profiler_AMRType", "TB_profiler_DBDate")

amr<-df$dr_variants

for (i in 1:nrow(amr)) {
  mutation<-paste(amr$gene[i],"_",amr$nucleotide_change[i], " Ratio:", amr$freq[i],sep = "")
  amrdummy<-amr$drugs[[i]]
  amrdummy$confers<-paste(toupper(amrdummy$confers), mutation, sep = " / ")
  amrdummy$type<-NULL
  amrdummy<-as.data.frame(t(amrdummy))
  colnames(amrdummy)<-paste("TB_profilerAMR_",toupper(amrdummy[1,]),sep = "")
  amrdummy<-amrdummy[-1,,drop=FALSE]
  if(!exists("amrout")){
    amrout<-amrdummy
  }else{
    amrout<- cbind(amrout,amrdummy)
  }
}


if(length(which(duplicated(colnames(amrout))))>0){
  dup.index<-which(duplicated(colnames(amrout)))
  dups<-colnames(amrout)[which(duplicated(colnames(amrout)))]
  
  for (i in 1:length(dups)) {
    merged_amr<-paste(amrout[1,which(colnames(amrout)==dups[i])], collapse =   " | ")
    amrout[1, which(colnames(amrout)==dups[i])[order(which(colnames(amrout)==dups[i]))[1]]]<-merged_amr
  }
  amrout<-amrout[,-dup.index]
}



results<-cbind(results,amrout)

write.table(results, paste(results$Sample,"_tb_profiler.tsv",sep = ""), row.names = FALSE, sep = "\t")

