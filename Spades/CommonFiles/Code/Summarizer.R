library(data.table)
library(writexl)
summaries<-list.files(pattern = "*contigs.stats.csv")
krakens<-list.files(pattern = "*resultskraken.csv")
fastqc<-list.files(pattern = "*fastqc.zip")


#Merging summaries
if(exists("summ")) rm(summ)
for (i in 1:length(summaries)) {
  dummy<-read.csv(summaries[i])
  if(!exists("summ")){
    summ<-dummy
  }else{
    summ<-rbind(summ,dummy)
  }
}

#Kraken
kraken.fasta<-krakens[grep("cleancontigs",krakens)]

summ$AgentContigs<-NA

for (i in 1:length(kraken.fasta)) {
  dummy<-read.csv(kraken.fasta[i])
  summ$AgentContigs[which(summ$Sample==gsub("_.*","",kraken.fasta[i]))]<-dummy$Specie[which(dummy$Ratio==max(dummy$Ratio))]
}

kraken.raw<-krakens[grep("Raw.resultskraken.csv",krakens)]
summ$AgentRawReads<-NA
summ$ReadsSupportingAgent<-NA
summ$RatioHumanReadsPreTrimming<-0
summ$RatioUnclassifiedReadsPreTrimming<-0

for (i in 1:length(kraken.raw)) {
  dummy<-read.csv(kraken.raw[i])
  summ$AgentRawReads[which(summ$Sample==gsub("_.*","",kraken.raw[i]))]<-dummy$Specie[which(dummy$Ratio==max(dummy$Ratio))]
  summ$ReadsSupportingAgent[which(summ$Sample==gsub("_.*","",kraken.raw[i]))]<-dummy$Ratio[which(dummy$Ratio==max(dummy$Ratio))]
  if(length(which(dummy$Specie=="unclassified (taxid 0)"))==1){
    summ$RatioUnclassifiedReadsPreTrimming[which(summ$Sample==gsub("_.*","",kraken.raw[i]))]<-dummy$Ratio[which(dummy$Specie=="unclassified (taxid 0)")]  
  }
  if(length(which(dummy$Specie=="Homo sapiens (taxid 9606)"))==1){
  summ$RatioHumanReadsPreTrimming[which(summ$Sample==gsub("_.*","",kraken.raw[i]))]<-dummy$Ratio[which(dummy$Specie=="Homo sapiens (taxid 9606)")]
  }
}

kraken.after<-krakens[grep("Trimmed.resultskraken.csv",krakens)]

summ$RatioHumanReadsAfterTrimming<-0
summ$RatioUnclassifiedReadsAfterTrimming<-0

for (i in 1:length(kraken.after)) {
  dummy<-read.csv(kraken.after[i])
  if(length(which(dummy$Specie=="unclassified (taxid 0)"))==1){
    summ$RatioUnclassifiedReadsAfterTrimming[which(summ$Sample==gsub("_.*","",kraken.fasta[i]))]<-dummy$Ratio[which(dummy$Specie=="unclassified (taxid 0)")]  
  }
  if(length(which(dummy$Specie=="Homo sapiens (taxid 9606)"))==1){
    summ$RatioHumanReadsAfterTrimming[which(summ$Sample==gsub("_.*","",kraken.fasta[i]))]<-dummy$Ratio[which(dummy$Specie=="Homo sapiens (taxid 9606)")]
  }
}


#Fastqc parsing

if(exists("out")) rm(out)
 if(!dir.exists("temp")) dir.create("temp")
for (i in 1:length(fastqc)) {

  unzip(fastqc[i], exdir ="temp" )  
  to.load<-list.files("temp", full.names = TRUE, pattern = "fastqc_data.txt", recursive = TRUE)
  dum.df<-fread(to.load, skip = "#Quality	Count", header = TRUE, sep = "\t" )
  dum.df<-as.data.frame(dum.df)
  colnames(dum.df)<-c("Q", "count")
  #dum.df<-dum.df[-1,]
  
  dum.df$file<-gsub("_fastqc.*","",gsub(".*/","", fastqc[i]))
  dum.df$count<-as.numeric(dum.df$count)
  dum.df$ratio<-dum.df$count/sum(dum.df$count)
  dum.df$TotalCount<-sum(dum.df$count)
  dum.df$Q30<-sum(dum.df$count[which(dum.df$Q>=30)])/sum(dum.df$count)
  
  unlink(paste("temp/",gsub(".zip","",fastqc[i]),sep = ""), recursive = TRUE)
  
  if(!exists("out")){
    out<-dum.df
  }else{
    out<-rbind(out,dum.df)
  }
  
}
unlink("temp", recursive = TRUE)
#Sample finding and summarizing
samps<-unique(summ$Sample)

summ$TotalReadCount<-NA
summ$RawQ30_R1<-NA
summ$RawQ30_R2<-NA
summ$ReadCountAfterTrimming<-NA
summ$TrimmedQ30_R1<-NA
summ$TrimmedQ30_R2<-NA
summ$RatioReadsTrimmed<-NA



for (i in 1:length(samps)) {
  out.samp<-out[grep(samps[i],out$file ),]
  summ$TotalReadCount[which(summ$Sample==samps[i])]<- unique(out.samp$TotalCount[grep("_R1_", out.samp$file)])
  summ$RawQ30_R1[which(summ$Sample==samps[i])]<- unique(out.samp$Q30[grep("_R1_", out.samp$file)])
  summ$RawQ30_R2[which(summ$Sample==samps[i])]<- unique(out.samp$Q30[grep("_R2_", out.samp$file)])
  
  summ$ReadCountAfterTrimming[which(summ$Sample==samps[i])]<- unique(out.samp$TotalCount[grep("_1P$", out.samp$file)])
  summ$TrimmedQ30_R1[which(summ$Sample==samps[i])]<- unique(out.samp$Q30[grep("_1P$", out.samp$file)])
  summ$TrimmedQ30_R2[which(summ$Sample==samps[i])]<- unique(out.samp$Q30[grep("_2P$", out.samp$file)])
  
  summ$RatioReadsTrimmed[which(summ$Sample==samps[i])]<- (summ$TotalReadCount[which(summ$Sample==samps[i])]- summ$ReadCountAfterTrimming[which(summ$Sample==samps[i])])/
    (summ$TotalReadCount[which(summ$Sample==samps[i])])
    
  #Folder Organization
  
  
  }


rmlist<-list.files(pattern = "_rmlst.csv")
if(exists("out.mlst")) rm(out.mlst)

for (i in 1:length(rmlist)) {
  dummy<-read.csv(rmlist[i])
  
  if(nrow(dummy)>1){
    dummy$MLST_Conflict<-paste(paste(dummy$rMLST_support, dummy$rMLST_taxon,sep = "% "), collapse = " / ")
    dummy<-dummy[which(dummy$rMLST_support==max(dummy$rMLST_support)),]
    if(nrow(dummy)>1) dummy<-dummy[1,]
  }
  
  
  if(!exists("out.mlst")){
    out.mlst<-dummy
  }else{
    if(length(setdiff(colnames(out.mlst), colnames(dummy) ))>0){
      dummy[setdiff(names(out.mlst), names(dummy))] <- NA
      out.mlst[setdiff(names(dummy), names(out.mlst))] <- NA  
    }
    out.mlst<-rbind(out.mlst, dummy)
  }
}
out.mlst$Sample<-gsub("_.*","",out.mlst$Sample)

rmlist<-list.files(pattern = "_seqmlst.csv")
if(exists("out.seqmlst")) rm(out.seqmlst)

for (i in 1:length(rmlist)) {
  dummy<-read.csv(rmlist[i])
  
  if(!exists("out.seqmlst")){
    out.seqmlst<-dummy
  }else{
    if(length(setdiff(colnames(out.seqmlst), colnames(dummy) ))>0){
      dummy[setdiff(names(out.seqmlst), names(dummy))] <- NA
      out.seqmlst[setdiff(names(dummy), names(out.seqmlst))] <- NA  
    }
    out.seqmlst<-rbind(out.seqmlst, dummy)
  }
}
out.seqmlst$Sample<-gsub("_.*","",out.seqmlst$Sample)

summ<-merge(summ, out.mlst, by="Sample", all.x=TRUE)
summ<-merge(summ, out.seqmlst, by="Sample", all.x=TRUE)

#Abricate
abrilist<-list.files(pattern = "_Abricate.csv")
for (i in 1:length(abrilist)) {
  dummy<-read.csv(abrilist[i])
  
  if(!exists("out.abri")){
    out.abri<-dummy
  }else{
    if(length(setdiff(colnames(out.abri), colnames(dummy) ))>0){
      dummy[setdiff(names(out.abri), names(dummy))] <- NA
      out.abri[setdiff(names(dummy), names(out.abri))] <- NA  
    }
    out.abri<-rbind(out.abri, dummy)
  }
}
summ<-merge(summ, out.abri, by="Sample", all.x=TRUE)

#HiCap

hicaplist<-list.files(pattern = "_HiCap.tsv")
for (i in 1:length(hicaplist)) {
  dummy<-read.csv(hicaplist[i], sep = "\t")
  #Fix no hicap
  if(colnames(dummy)[1]=="NoHi"){ 
    dummy<-as.data.frame(matrix(NA))
    colnames(dummy)[1]<-"predicted_serotype"
  }
  if(!exists("out.hicap")){
    out.hicap<-dummy
  }else{
    if(length(setdiff(colnames(out.hicap), colnames(dummy) ))>0){
      dummy[setdiff(names(out.hicap), names(dummy))] <- NA
      out.hicap[setdiff(names(dummy), names(out.hicap))] <- NA  
    }
    out.hicap<-rbind(out.hicap, dummy)
  }
}
out.hicap$Sample<-gsub("_.*","",hicaplist)
out.hicap<-out.hicap[, which(colnames(out.hicap) %in% c("predicted_serotype", "Sample"))]
colnames(out.hicap)[which(colnames(out.hicap)=="predicted_serotype")]<-"HiCap_Predicted_Serotype"

summ<-merge(summ, out.hicap, by="Sample", all.x=TRUE)

#Seroba

serolist<-list.files(pattern = "_seroba.tsv")
if(exists("out.seroba")) rm(out.seroba)
for (i in 1:length(serolist)) {
  dummy<-read.csv(serolist[i], sep = "\t",header = FALSE)
  #Fix no hicap
  if(dummy[1,1]=="NoSpne"){ 
    dummy<-as.data.frame(matrix(NA))
    colnames(dummy)[1]<-"Seroba_Predicted_Serotype"
    dummy$Sample<-gsub("_.*","",serolist[i])
  }else{
    if(ncol(dummy)==2) colnames(dummy)<-c("Sample", "Seroba_Predicted_Serotype")
    if(ncol(dummy)==3) colnames(dummy)<-c("Sample", "Seroba_Predicted_Serotype","Seroba_warning")
  }
  
  if(!exists("out.seroba")){
    out.seroba<-dummy
  }else{
    if(length(setdiff(colnames(out.seroba), colnames(dummy) ))>0 | 
       length(setdiff(colnames(dummy), colnames(out.seroba) ))>0){
      dummy[setdiff(names(out.seroba), names(dummy))] <- NA
      out.seroba[setdiff(names(dummy), names(out.seroba))] <- NA  
    }
    out.seroba<-rbind(out.seroba, dummy)
  }
}

out.seroba$Sample<-gsub("_.*","",out.seroba$Sample)
summ<-merge(summ, out.seroba, by="Sample", all.x=TRUE)




write.csv(summ, paste("Summaries_",gsub("-","",Sys.Date()), ".csv",sep = ""), row.names = FALSE)
write_xlsx(summ, paste("Summaries_",gsub("-","",Sys.Date()), ".xlsx",sep = ""))
