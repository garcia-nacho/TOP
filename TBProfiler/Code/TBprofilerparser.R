library(jsonlite)

jsonfile<-list.files(pattern = "tb_profiler.json")

df<-fromJSON(jsonfile)

results<-c(df$id, df$pipeline$software_version, df$qc$percent_reads_mapped,
           df$qc$genome_median_depth,df$main_lineage,df$sub_lineage,df$drtype ,df$pipeline$db_version$Date)

results<-as.data.frame(t(results))
colnames(results)<-c("Sample", "TB_profiler_Version", "TB_profiler_FractionReadsMapped","TB_profiler_MedianCoverage", "TB_profiler_MainLineage",
                     "TB_profiler_Sublineage", "TB_profiler_AMRType", "TB_profiler_DBDate")


if(length(df$dr_variants)>0){
  amr<-df$dr_variants
for (i in 1:nrow(amr)) {
  mutation<-paste(amr$gene_name[i],"_",amr$nucleotide_change[i], " Ratio:", amr$freq[i],sep = "")
  amrdummy<-amr$drugs[[i]]
  
  amrdummy$confers<-paste(toupper(amrdummy$type), mutation, sep = " / ")
  amrdummy<-amrdummy[,c("drug","confers")]
  
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
}


drug_list<-list.files(pattern ="tb_profiler.csv" )
druglines<-readLines(drug_list)
star.read<-which(druglines=="Drug,Genotypic Resistance,Mechanisms" )+1
end.read<-which(druglines=="Resistance variants report")-2

drug.vector<-vector()
for (i in c(star.read:end.read)) {
  drug.vector<- c(drug.vector,gsub(",.*","",druglines[i])  )
}


drugs<-paste("TB_profilerAMR_",toupper(as.character(drug.vector)),sep = "")

if(length(which(drugs %in% colnames(results))))drugs<-drugs[-which(drugs %in% colnames(results))]
if(length(drugs)>0){
  pad.drugs<-as.data.frame(matrix("SENSITIVE", ncol = length(drugs), nrow = 1))
  colnames(pad.drugs)<-drugs
  results<-cbind(results, pad.drugs)
}


write.table(results, paste(results$Sample,"_tb_profiler.tsv",sep = ""), row.names = FALSE, sep = "\t")

