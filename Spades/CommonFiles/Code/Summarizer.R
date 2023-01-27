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
  summ$AgentContigs[which(summ$Sample==gsub("_.*","",kraken.fasta[i]))]<-dummy$Specie[which(dummy$Ratio==max(dummy$Ratio))][1]
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
    if(length(setdiff(colnames(out.mlst), colnames(dummy) ))>0 | length(setdiff(colnames(dummy), colnames(out.mlst) ))>0){
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
    if(length(setdiff(colnames(out.abri), colnames(dummy) ))>0 |
       length(setdiff(colnames(dummy), colnames(out.abri) ))>0){
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



#EMM

emlist<-list.files(pattern = "_emmtyper.tsv")

for (i in 1:length(emlist)) {
  dummy<-read.csv(emlist[i],sep = "\t",header = FALSE)
  if(dummy[1,1]=="NoSpy"){ 
    dummy<-as.data.frame(matrix(NA,nrow = 1, ncol = 5))
    colnames(dummy)<-c("Sample","Cluster.N","emm-type","emm-like","emm-cluster")
    dummy$emm.warning<-NA
    dummy$Sample<-gsub("_.*","",emlist[i])
  }else{
  
  
  colnames(dummy)<-c("Sample","Cluster.N","emm-type","emm-like","emm-cluster")
  dummy$emm.warning<-NA
  if(as.numeric(dummy$Cluster.N)>3) dummy$emm.warning<-"Possible Contamination"

  }
  if(!exists("out.emm")){
    out.emm<-dummy
  }else{
    out.emm<-rbind(out.emm,dummy)
  }
}

out.emm$Sample<-gsub("_.*","",out.emm$Sample)
summ<-merge(summ, out.emm, by="Sample", all.x=TRUE)


# STXTyping ---------------------------------------------------------------
library(jsonlite)
stxlist<-list.files(pattern = "_virfinder.json")
for (i in 1:length(stxlist)) {
  if(exists("dum.json")) rm(dum.json)
  try(dum.json<-fromJSON(stxlist[i]),silent = TRUE)
  if(!exists("dum.json")) dum.json<-vector()
  if(exists("stx.table.out")) rm(stx.table.out)
  if(exists("output.stx"))rm(output.stx)
  if(length(dum.json)>0){
    output.stx<-dum.json$virulencefinder$results$`Escherichia coli`$stx  
  }else{
    output.stx<-vector()
  }
  
  if(length(output.stx)>0){
    for (j in 1:length(output.stx)) {
      stx.table<-as.data.frame(t(unlist(output.stx[[j]])))
      if(!exists("stx.table.out")){
        stx.table.out<-stx.table
      }else{
        stx.table.out<-rbind(stx.table.out,stx.table)
      }
    }
   stx.table.out$Sample<-gsub("_.*","",stxlist[i])
   stx.table.out<-stx.table.out[c(1,2,10,13)]
   colnames(stx.table.out)<-c("StxType","StxIdentity","StxInfo","Sample")
   write.csv(stx.table.out, paste(gsub("_.*","",stxlist[i]),"_STXType.csv",sep = ""),row.names = FALSE)
   if(nrow(stx.table.out)>1){
     stx.table.out$StxType[1]<-paste(stx.table.out$StxType, collapse ="/" )
     stx.table.out$StxIdentity[1]<-paste(stx.table.out$StxIdentity, collapse ="/" )
     stx.table.out$StxInfo[1]<-paste(stx.table.out$StxInfo, collapse ="/" )
     stx.table.out<-stx.table.out[1,,drop=FALSE]
   }
  }else{
    stx.table.out<-as.data.frame(matrix(NA, ncol =4, nrow = 1 ))
    colnames(stx.table.out)<-c("StxType","StxIdentity","StxInfo","Sample")
    stx.table.out$Sample<-gsub("_.*","",stxlist[i])
  }
 if(!exists("stx.out.final")){
   stx.out.final<-stx.table.out
 }else{
   stx.out.final<-rbind(stx.out.final, stx.table.out)
 }
  }


summ<-merge(summ, stx.out.final, by="Sample", all.x=TRUE)
rm(stx.out.final)
# VFyping ---------------------------------------------------------------
library(jsonlite)
stxlist<-list.files(pattern = "_virfinder.json")
for (i in 1:length(stxlist)) {
  if(exists("dum.json")) rm(dum.json)
  try(dum.json<-fromJSON(stxlist[i]),silent = TRUE)
  if(!exists("dum.json")) dum.json<-vector()
  if(exists("stx.table.out")) rm(stx.table.out)
  if(exists("output.stx"))rm(output.stx)
  if(length(dum.json)>0){
    output.stx<-dum.json$virulencefinder$results$`Escherichia coli`$virulence_ecoli 
  }else{
    output.stx<-vector()
  }
  
  if(length(output.stx)>0){
    for (j in 1:length(output.stx)) {
      stx.table<-as.data.frame(t(unlist(output.stx[[j]])))
      if(!exists("stx.table.out")){
        stx.table.out<-stx.table
      }else{
        stx.table.out<-rbind(stx.table.out,stx.table)
      }
    }
    stx.table.out$Sample<-gsub("_.*","",stxlist[i])
    stx.table.out<-stx.table.out[c(1,2,10,13)]
    colnames(stx.table.out)<-c("VirulenceFactorEcoli","VFEcoli_Identity","VFEcoli_Info","Sample")
    write.csv(stx.table.out, paste(gsub("_.*","",stxlist[i]),"_Virulencefactors.csv",sep = ""),row.names = FALSE)
    if(nrow(stx.table.out)>1){
      stx.table.out$VirulenceFactorEcoli[1]<-paste(stx.table.out$VirulenceFactorEcoli, collapse ="/" )
      stx.table.out$VFEcoli_Identity[1]<-paste(stx.table.out$VFEcoli_Identity, collapse ="/" )
      stx.table.out$VFEcoli_Info[1]<-paste(stx.table.out$VFEcoli_Info, collapse ="/" )
      stx.table.out<-stx.table.out[1,,drop=FALSE]
    }
    
  }else{
    stx.table.out<-as.data.frame(matrix(NA, ncol =4, nrow = 1 ))
    colnames(stx.table.out)<-c("VirulenceFactorEcoli","VFEcoli_Identity","VFEcoli_Info","Sample")
    stx.table.out$Sample<-gsub("_.*","",stxlist[i])
  }
  if(!exists("stx.out.final")){
    stx.out.final<-stx.table.out
  }else{
    stx.out.final<-rbind(stx.out.final, stx.table.out)
  }
}

summ<-merge(summ, stx.out.final, by="Sample", all.x=TRUE)
# Last Stage --------------------------------------------------------------


empty.col<-apply(summ, 2, function(x) length(which(is.na(x))))
empty.col<-as.numeric(empty.col)
empty.col<-as.numeric(which(empty.col==nrow(summ)))

if(length(empty.col)>0) summ<-summ[,-empty.col]

write.csv(summ, paste("Summaries_",gsub("-","",Sys.Date()), ".csv",sep = ""), row.names = FALSE)
write_xlsx(summ, paste("Summaries_",gsub("-","",Sys.Date()), ".xlsx",sep = ""))
