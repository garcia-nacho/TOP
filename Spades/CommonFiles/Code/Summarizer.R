library(data.table)
library(writexl)
library(jsonlite)
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

summ$GenomeSize<-paste(round(summ$Coverage/1000000,2),"Mbp",sep = "")

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
  out.samp<-out[grep(paste(samps[i],"_",sep = ""),out$file ),]
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
    dummy$rMLST_Conflict<-paste(paste(dummy$rMLST_support, dummy$rMLST_taxon,sep = "% "), collapse = " / ")
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
    if(length(setdiff(colnames(out.seqmlst), colnames(dummy) ))>0 | 
       length(setdiff(colnames(dummy), colnames(out.seqmlst) ))>0){
      dummy[setdiff(names(out.seqmlst), names(dummy))] <- NA
      out.seqmlst[setdiff(names(dummy), names(out.seqmlst))] <- NA  
    }
    out.seqmlst<-rbind(out.seqmlst, dummy)
  }
}
out.seqmlst$Sample<-gsub("_.*","",out.seqmlst$Sample)

summ<-merge(summ, out.mlst, by="Sample", all.x=TRUE)
summ<-merge(summ, out.seqmlst, by="Sample", all.x=TRUE)

#local mlst


localmlst<-list.files(pattern = "localmlst.tsv")

if(exists("out.localmlst")) rm(out.localmlst)

for (i in 1:length(localmlst)) {

  dummy<-read.csv(localmlst[i], header = FALSE,sep = "\t")
  dummy2<-as.data.frame(matrix(NA, nrow = 1, ncol = 3))
  colnames(dummy2)<-c("Sample","LocalMLST","LocalScheme")
  dummy2$LocalMLST<-dummy$V3
  
  if(ncol(dummy)>3){
  scheme<-as.character(dummy[1,c(4:ncol(dummy))])
  schemenames<-gsub("\\(.*","",scheme)
  schemetype<-gsub(".*\\(","",gsub("\\)","",scheme))
  schememerged<-paste(schemenames, schemetype, sep = ":")
  schememerged<-schememerged[order(schememerged)]
  schememerged<-paste(schememerged, collapse = " | ")
  dummy2$LocalScheme<-schememerged
  dummy2$Sample<-dummy$V1
  }else{
    dummy2$LocalScheme<-NA
  }
  if(dummy$V3=="-") dummy2$LocalMLST<-NA
  if(!exists("out.localmlst")){
    out.localmlst<-dummy2
  }else{
    out.localmlst<-rbind(out.localmlst, dummy2)
  }
}

out.localmlst$Sample<-gsub("_.*","",out.localmlst$Sample)

summ<-merge(summ, out.localmlst, by="Sample", all.x=TRUE)

to.replace<-which(is.na(summ$MLST.Type))
summ$MLST.Type[to.replace]<- summ$LocalMLST[to.replace]
summ$MLST.Scheme[to.replace]<- summ$LocalScheme[to.replace]
summ$MLST_Date[to.replace]<-"MLST_Database20240116"
summ$LocalMLST<-NULL
summ$LocalScheme<-NULL

if(length(which(is.na(summ$MLST.Type)))>0 ) summ$MLST.Type[which(is.na(summ$MLST.Type))]<-"Non Detected"
if(length(which(is.na(summ$MLST.Scheme)))>0 ) summ$MLST.Scheme[which(is.na(summ$MLST.Scheme))]<-"Non Detected"


#Abricate
abrilist<-list.files(pattern = "_Abricate.csv")
for (i in 1:length(abrilist)) {
  dummy<-read.csv(abrilist[i])
  if(length(grep("^VF_",colnames(dummy)))>0) dummy<-dummy[,-grep("^VF_",colnames(dummy))]
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
    if(length(setdiff(colnames(out.hicap), colnames(dummy) ))>0 | 
       length(setdiff(colnames(dummy), colnames(out.hicap) ))>0){
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
patterns<- c("_fastq_virfinder.json", "_contigs_virfinder.json")
lab<-c("Mapping", "Contigs")

for (hl in 1:length(patterns)) {
  
  if(exists("stx.table.out")) rm(stx.table.out)
  if(exists("stx.out.final"))rm(stx.out.final)
  stxlist<-list.files(pattern = patterns[hl])
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
          
          if(length(setdiff(colnames(stx.table.out), colnames(stx.table) ))>0 | 
             length(setdiff(colnames(stx.table), colnames(stx.table.out) ))>0){
            stx.table[setdiff(names(stx.table.out), names(stx.table))] <- NA
            stx.table.out[setdiff(names(stx.table), names(stx.table.out))] <- NA
          }
          
          
          stx.table.out<-rbind(stx.table.out,stx.table)
        }
      }
      stx.table.out$Sample<-gsub("_.*","",stxlist[i])
      
      if(stx.table.out[1,1]=="No hit found"){
        colnames(stx.table.out)[1]<-"StxType"
        stx.table.out$StxIdentity<-NA
        stx.table.out$StxInfo<-NA
        stx.table.out$StxVariant<-NA
        stx.table.out$StxCoverage<-NA
      }else{
        stx.table.out<-stx.table.out[c(1,2,10,11,13)]
        colnames(stx.table.out)<-c("StxType","StxIdentity","StxInfo", "StxCoverage", "Sample")
        stx.table.out$StxVariant<-paste(stx.table.out$StxType, gsub(".*, ","",stx.table.out$StxInfo ),sep = ":")
        stx.table.out$StxInfo<-gsub(", .*","",stx.table.out$StxInfo )
      }
      
      
      if(nrow(stx.table.out)>1){
        stx.table.out$StxType[1]<-paste(stx.table.out$StxType, collapse =" | " )
        stx.table.out$StxIdentity[1]<-paste(stx.table.out$StxIdentity, collapse =" | " )
        stx.table.out$StxInfo[1]<-paste(stx.table.out$StxInfo, collapse =" | " )
        stx.table.out$StxVariant[1]<-paste(stx.table.out$StxVariant, collapse =" | " )
        stx.table.out$StxCoverage[1]<-paste(stx.table.out$StxCoverage, collapse =" | " )
        stx.table.out<-stx.table.out[1,,drop=FALSE]
      }
    }else{
      stx.table.out<-as.data.frame(matrix(NA, ncol =5, nrow = 1 ))
      colnames(stx.table.out)<-c("StxType","StxIdentity","StxInfo", "StxCoverage", "Sample")
      stx.table.out$Sample<-gsub("_.*","",stxlist[i])
    }
    write.csv(stx.table.out, paste(gsub("_.*","",stxlist[i]),"_STXType.csv",sep = ""),row.names = FALSE)
    
    if(!exists("stx.out.final")){
      stx.out.final<-stx.table.out
    }else{
      if(length(setdiff(colnames(stx.table.out), colnames(stx.out.final) ))>0 | 
         length(setdiff(colnames(stx.out.final), colnames(stx.table.out) ))>0){
        stx.out.final[setdiff(names(stx.table.out), names(stx.out.final))] <- NA
        stx.table.out[setdiff(names(stx.out.final), names(stx.table.out))] <- NA
      }
      
      stx.out.final<-rbind(stx.out.final, stx.table.out)
    }
  }
  


#stx.out.final<-stx.out.final[,c(1,5,2,3,4)]

colnames(stx.out.final)[-which(colnames(stx.out.final)=="Sample")]<-
  paste(colnames(stx.out.final)[-which(colnames(stx.out.final)=="Sample")], lab[hl],sep = "_")

summ<-merge(summ, stx.out.final, by="Sample", all.x=TRUE)
rm(stx.out.final)

}

# VFyping ---------------------------------------------------------------

patterns<- c("_fastq_virfinder.json", "_contigs_virfinder.json")
lab<-c("Mapping", "Contigs")

for (hl in 1:length(patterns)) {

stxlist<-list.files(pattern = patterns[hl])
if(exists("stx.table.out")) rm(stx.table.out)
if(exists("stx.out.final"))rm(stx.out.final)
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
    
    if(stx.table.out[1,1]=="No hit found"){
      colnames(stx.table.out)[1]<-"VirulenceFactorEcoli"
      stx.table.out$VFEcoli_Identity<-NA
      stx.table.out$VFEcoli_Info<-NA
    }else{
      stx.table.out<-stx.table.out[c(1,2,10,13)]
      colnames(stx.table.out)<-c("VirulenceFactorEcoli","VFEcoli_Identity","VFEcoli_Info","Sample")
    }
    
    
    write.csv(stx.table.out, paste(gsub("_.*","",stxlist[i]),"_Virulencefactors.csv",sep = ""),row.names = FALSE)
    stx.table.out$VFEcoli_Info<-gsub("^ +","",stx.table.out$VFEcoli_Info) 
    
    if(nrow(stx.table.out)>1){
      stx.table.out$VirulenceFactorEcoli[1]<-paste(stx.table.out$VirulenceFactorEcoli, collapse =" | " )
      stx.table.out$VFEcoli_Identity[1]<-paste(stx.table.out$VFEcoli_Identity, collapse =" | " )
      stx.table.out$VFEcoli_Info[1]<-paste(stx.table.out$VFEcoli_Info, collapse =" | " )
      
      stx.table.out<-stx.table.out[1,,drop=FALSE]
    }
    
  }else{
    stx.table.out<-as.data.frame(matrix(NA, ncol =4, nrow = 1 ))
    colnames(stx.table.out)<-c("VirulenceFactorEcoli","VFEcoli_Identity","VFEcoli_Info","Sample")
    stx.table.out$Sample<-gsub("_.*","",stxlist[i])
    write.csv(stx.table.out, paste(gsub("_.*","",stxlist[i]),"_Virulencefactors.csv",sep = ""),row.names = FALSE)
    
  }
  if(!exists("stx.out.final")){
    stx.out.final<-stx.table.out
  }else{
    stx.out.final<-rbind(stx.out.final, stx.table.out)
  }
}

colnames(stx.out.final)[-which(colnames(stx.out.final)=="Sample")]<-
  paste(colnames(stx.out.final)[-which(colnames(stx.out.final)=="Sample")], lab[hl],sep = "_")

summ<-merge(summ, stx.out.final, by="Sample", all.x=TRUE)
}

#summary.txt in .fastq.gz. If Basic statistics = PASS -> OK else -> WARN
sum.txts<-list.files(recursive = FALSE,pattern = "Bowtie2summary.txt")
for (i in 1:length(sum.txts)) {
  sumdum<-read.csv(sum.txts[i],header = FALSE)
  qual<- sumdum[grep("overall alignment rate",sumdum$V1),]
  if(length(qual)>0){
    qual<-(100-as.numeric(gsub("%.*","",qual)))/100
    qual<-as.data.frame(qual)
    colnames(qual)<-"UnmappedFraction"
    qual$Sample<-gsub("_.*","",sum.txts[i])
  }else{
    qual<-as.data.frame(NA)
    colnames(qual)<-"UnmappedFraction"
    qual$Sample<-gsub("_.*","",sum.txts[i])
  }
  if(!exists("qualout")){
    qualout<-qual
  }else{
    qualout<-rbind(qualout,qual)
  }
}


summ<-merge(summ, qualout, by="Sample", all.x=TRUE)



# MeningoType ------------------------------------------------------------
meningo.list<-list.files(pattern = "_meningotype.txt")

for (i in 1:length(meningo.list)) {
  dummy<-read.csv(meningo.list[i],sep = "\t",header = FALSE)
  if(dummy[1,1]=="NoNmen"){ 
    dummy<-as.data.frame(matrix(NA,nrow = 1, ncol = 11))
    colnames(dummy)<-c("SAMPLE_ID", 
                       "SEROGROUP",
                       "ctrA",
                       "MLST",
                       "PorA",
                       "FetA",
                       "PorB",
                       "fHbp" ,
                       "NHBA",
                       "NadA"   ,
                       "BAST"   )

    dummy$SAMPLE_ID<-gsub("_.*","",meningo.list[i])
  }else{
    colnames(dummy)<-dummy[1,]
    dummy<-dummy[-1,]
    dummy$PorA<-paste("P1.", dummy$PorA,sep="")
  }
  if(!exists("out.meningotype")){
    out.meningotype<-dummy
  }else{
    out.meningotype<-rbind(out.meningotype,dummy)
  }
}

out.meningotype$Sample<-gsub("_.*","",out.meningotype$SAMPLE_ID)
out.meningotype$SAMPLE_ID<-NULL
colnames(out.meningotype)[-which(colnames(out.meningotype)=="Sample")]<-
  paste("Meningotype_",colnames(out.meningotype)[-which(colnames(out.meningotype)=="Sample")],sep = "")
summ<-merge(summ, out.meningotype, by="Sample", all.x=TRUE)


# NGMaster ----------------------------------------------------------------
gono.list<-list.files(pattern = "_ngmast_results.txt")
for (i in 1:length(gono.list)) {
  dummy<-read.csv(gono.list[i],sep = "\t",header = FALSE)
  if(dummy[1,1]=="NoNgon"){
    dummy<-as.data.frame(matrix(NA,nrow = 1, ncol = 4))
    colnames(dummy)<-c("ID",
                       "NG-MAST",
                       "POR",
                       "TBPB")
    
    dummy$ID<-gsub("_.*","",gono.list[i])
  }else{
    colnames(dummy)<-dummy[1,]
    dummy<-dummy[-1,]
  }
  if(!exists("out.gono")){
    out.gono<-dummy
  }else{
    out.gono<-rbind(out.gono,dummy)
  }
}

out.gono$Sample<-gsub("_.*","",out.gono$ID)
out.gono$ID<-NULL
colnames(out.gono)[-which(colnames(out.gono)%in% c("Sample","NG-MAST"))]<-
  paste("NG-MAST_",colnames(out.gono)[-which(colnames(out.gono)%in% c("Sample","NG-MAST"))],sep = "")
summ<-merge(summ, out.gono, by="Sample", all.x=TRUE)
rm(out.gono)



# NGStar ------------------------------------------------------------------

gono.list2<-list.files(pattern = "ngstar_results.tsv")
for (i in 1:length(gono.list2)) {
  dummy<-read.csv(gono.list2[i],sep = "\t",header = FALSE)
  if(dummy[1,1]=="NoNgon"){
    dummy<-as.data.frame(matrix(NA,nrow = 1, ncol = 1))
    colnames(dummy)<-c("strain")
    
    dummy$ID<-gsub("_.*","",gono.list[i])
  }else{
    colnames(dummy)<-dummy[1,]
    dummy<-dummy[-1,]
  }
  if(!exists("out.gono2")){
    out.gono2<-dummy
  }else{
    
    missing.dummy<- colnames(out.gono2)[-which(colnames(out.gono2) %in% colnames(dummy) )]
    if(length(missing.dummy)>0){
      padd<-as.data.frame(matrix(NA, nrow = nrow(dummy), ncol = length(missing.dummy)))
      colnames(padd)<-missing.dummy
      dummy<-cbind(dummy, padd)
    }
    
    missing.gono<- colnames(dummy)[-which(colnames(dummy) %in% colnames(out.gono2) )]
    if(length(missing.gono)>0){
      padd<-as.data.frame(matrix(NA, nrow = nrow(out.gono2), ncol = length(missing.gono)))
      colnames(padd)<-missing.gono
      out.gono2<-cbind(out.gono2, padd)
    }
    
    out.gono2<-rbind(out.gono2,dummy)
  }
}

out.gono2$Sample<-gsub("_.*","",out.gono2$strain)
out.gono2$strain<-NULL
colnames(out.gono2)[-which(colnames(out.gono2)%in% c("Sample"))]<-
  paste("NGstar_",colnames(out.gono2)[-which(colnames(out.gono2)%in% c("Sample"))],sep = "")
summ<-merge(summ, out.gono2, by="Sample", all.x=TRUE)
rm(out.gono2)



# Seqsero -----------------------------------------------------------------
salmo.list<-list.files(pattern = "_seqsero_results.tsv")
if(length(grep("kmer",salmo.list))>0) salmo.list<-salmo.list[-grep("kmer",salmo.list)]

for (i in 1:length(salmo.list)) {
  dummy<-read.csv(salmo.list[i],sep = "\t",header = FALSE)
  if(dummy[1,1]=="NoSalmo"){
    dummy<-as.data.frame(matrix(NA,nrow = 1, ncol = 8))
    colnames(dummy)<-c( "O antigen prediction",
                       "H1 antigen prediction(fliC)",
                       "H2 antigen prediction(fljB)",
                       "Predicted identification",
                       "Predicted antigenic profile",
                       "Predicted serotype",
                       "Potential inter-serotype contamination",
                       "Note" )
    
    dummy$ID<-gsub("_.*","",salmo.list[i])
  }else{
    colnames(dummy)<-dummy[1,]
    dummy<-dummy[-1,]
    dummy<-dummy[,-which(colnames(dummy) %in% c("Output directory","Input files"))]
    colnames(dummy)[which(colnames(dummy)=="Sample name")]<-"ID"
  }
  if(!exists("out.salmo")){
    out.salmo<-dummy
  }else{
    out.salmo<-rbind(out.salmo,dummy)
  }
}

out.salmo$Sample<-gsub("_.*","",out.salmo$ID)
out.salmo$ID<-NULL
colnames(out.salmo)<-gsub(" ","_",colnames(out.salmo))
colnames(out.salmo)[-which(colnames(out.salmo)=="Sample")]<-
  paste("SeqSeroAllele_",colnames(out.salmo)[-which(colnames(out.salmo)=="Sample")],sep = "")
summ<-merge(summ, out.salmo, by="Sample", all.x=TRUE)
rm(out.salmo)


# SeqseroKmer -------------------------------------------------------------
salmo.list<-list.files(pattern = "_kmer_seqsero_results.tsv")

for (i in 1:length(salmo.list)) {
  dummy<-read.csv(salmo.list[i],sep = "\t",header = FALSE)
  if(dummy[1,1]=="NoSalmo"){
    dummy<-as.data.frame(matrix(NA,nrow = 1, ncol = 7))
    colnames(dummy)<-c( "O antigen prediction",
                        "H1 antigen prediction(fliC)",
                        "H2 antigen prediction(fljB)",
                        "Predicted identification",
                        "Predicted antigenic profile",
                        "Predicted serotype",
                        "Note" )
    
    dummy$ID<-gsub("_.*","",salmo.list[i])
  }else{
    colnames(dummy)<-dummy[1,]
    dummy<-dummy[-1,]
    dummy<-dummy[,-which(colnames(dummy) %in% c("Output directory","Input files"))]
    colnames(dummy)[which(colnames(dummy)=="Sample name")]<-"ID"
  }
  if(!exists("out.salmo")){
    out.salmo<-dummy
  }else{
    out.salmo<-rbind(out.salmo,dummy)
  }
}

out.salmo$Sample<-gsub("_.*","",out.salmo$ID)
out.salmo$ID<-NULL
colnames(out.salmo)<-gsub(" ","_",colnames(out.salmo))
colnames(out.salmo)[-which(colnames(out.salmo)=="Sample")]<-
  paste("SeqSeroKmer_",colnames(out.salmo)[-which(colnames(out.salmo)=="Sample")],sep = "")
summ<-merge(summ, out.salmo, by="Sample", all.x=TRUE)
rm(out.salmo)


# Tartrate  -----------------------------------------------------------------
tart.list<-list.files(pattern = "_tartrate.txt")
if(exists("out.tart")) rm(out.tart)
for (i in 1:length(tart.list)) {
  dummy<-read.csv(tart.list[i],sep = "-", header = FALSE)
  
  
  if(dummy[1,1]=="NoSalmo"){
    dummy<-as.data.frame(matrix(NA,nrow = 1, ncol = 4))
    colnames(dummy)<-c( "ID",
                        "Tartratre_STM3356_Codon",
                        "Predicted_Tartrate_Fermentation",
                        "Tartrate_Notes")
    
    dummy$ID<-gsub("_.*","",tart.list[i])
  }else{
    if(ncol(dummy)>4){
      dummy2<-dummy[,c((ncol(dummy)-2):ncol(dummy))]
      dummy2<-cbind(gsub(" ","",paste(dummy[,c(1:(ncol(dummy)-3) )], collapse = "-")), dummy2) 
      dummy<-dummy2
    } 
    colnames(dummy)<-c( "ID",
                        "Tartratre_STM3356_Codon",
                        "Predicted_Tartrate_Fermentation",
                        "Tartrate_Notes")
  }
  if(!exists("out.tart")){
    out.tart<-dummy
  }else{
    out.tart<-rbind(out.tart,dummy)
  }
}

out.tart$Sample<-gsub("_.*","",out.tart$ID)
out.tart$ID<-NULL

summ<-merge(summ, out.tart, by="Sample", all.x=TRUE)
rm(out.tart)

summ$Tartratre_STM3356_Codon<-gsub("^functional", "Non functional", summ$Tartratre_STM3356_Codon)
summ$Tartratre_STM3356_Codon<-gsub("^ Intact", "Intact", summ$Tartratre_STM3356_Codon)
summ$Tartratre_STM3356_Codon<-gsub(" $", "", summ$Tartratre_STM3356_Codon)


# TB_Pipeline -------------------------------------------------------------
tbp.list<-list.files(pattern = "_tbp.csv")
mash.list<-list.files(pattern = "mashmykrobe.csv")
mashtb<-read.csv(mash.list,sep = "\t")

if(length(tbp.list)>0){
if(exists("out.tbp")) rm(out.tbp)
for (i in 1:length(tbp.list)) {
  dummy<-read.csv(tbp.list[i], header = TRUE)
  if(colnames(dummy)[1]=="NoMyco"){
    dummy<-as.data.frame(matrix(NA,nrow = 1, ncol = 1))
    colnames(dummy)<-"Sample"
    dummy$Sample<-gsub("_.*","",tbp.list[i])
  }
  if(!exists("out.tbp")){
    out.tbp<-dummy
  }else{
    
    missing.dummy<- colnames(out.tbp)[-which(colnames(out.tbp) %in% colnames(dummy) )]
    if(length(missing.dummy)>0){
      padd<-as.data.frame(matrix(NA, nrow = nrow(dummy), ncol = length(missing.dummy)))
      colnames(padd)<-missing.dummy
      dummy<-cbind(dummy, padd)
    }
    
    missing.tbp<- colnames(dummy)[-which(colnames(dummy) %in% colnames(out.tbp) )]
    if(length(missing.tbp)>0){
      padd<-as.data.frame(matrix(NA, nrow = nrow(out.tbp), ncol = length(missing.tbp)))
      colnames(padd)<-missing.tbp
      out.tbp<-cbind(out.tbp, padd)
    }
    
    out.tbp<-rbind(out.tbp,dummy)
  }
  
}

out.tbp$Sample<-gsub("_.*","",out.tbp$Sample)
out.tbp$TB_lineage<-NULL


if(mashtb[1,1]!="NonTB_in_the_run"){
  colnames(mashtb)[-1]<-paste("TB_",colnames(mashtb)[-1],sep = "")
  out.tbp<-merge(mashtb, out.tbp,by.x = "samples", by.y="Sample", all.x = TRUE)
}

summ<-merge(summ, out.tbp, by.x="Sample", by.y="samples", all.x=TRUE)
rm(out.tbp)
}

summ$TB_phylo_group<-NULL


# TBProfiler --------------------------------------------------------------

drugs<-toupper(c("Rifampicin",
"Isoniazid",
"Ethambutol",
"Pyrazinamide",
"Moxifloxacin",
"Levofloxacin",
"Bedaquiline",
"Delamanid",
"Pretomanid",
"Linezolid",
"Streptomycin",
"Amikacin",
"Kanamycin",
"Capreomycin",
"Clofazimine",
"Ethionamide",
"Para-aminosalicylic_acid",
"Cycloserine"))



tbp.list<-list.files(pattern = "_tb_profiler.tsv")
if(exists("out.tbp")) rm(out.tbp)

for (i in 1:length(tbp.list)) {
  dummy<-read.csv(tbp.list[i], sep = "\t")
  if(colnames(dummy)[1]!="dummy"){

    dummy$Sample<-gsub("_.*","",dummy$Sample)
    if(!exists("out.tbp")){
      out.tbp<-dummy
    }else{
      missing.dummy<- colnames(out.tbp)[-which(colnames(out.tbp) %in% colnames(dummy) )]
      if(length(missing.dummy)>0){
        padd<-as.data.frame(matrix(NA, nrow = nrow(dummy), ncol = length(missing.dummy)))
        colnames(padd)<-missing.dummy
        dummy<-cbind(dummy, padd)
      }
      
      missing.tbp<- colnames(dummy)[-which(colnames(dummy) %in% colnames(out.tbp) )]
      if(length(missing.tbp)>0){
        padd<-as.data.frame(matrix(NA, nrow = nrow(out.tbp), ncol = length(missing.tbp)))
        colnames(padd)<-missing.tbp
        out.tbp<-cbind(out.tbp, padd)
      }
      out.tbp<-rbind(out.tbp, dummy)
    }
  }
}

if(exists("out.tbp")){
  summ<-merge(summ, out.tbp, by="Sample", all.x=TRUE)
}


# EcoliPipeline -----------------------------------------------------------

#Raw 
eco.list<-list.files(pattern = "raw_ecopipeline.csv")

if(exists("out.eco")) rm(out.eco)
for (i in 1:length(eco.list)) {
  dummy<-read.csv(eco.list[i])
  if(colnames(dummy)[1]=="NoEcol"){
    dummy<-as.data.frame(matrix(NA,nrow = 1, ncol = 1))
    colnames(dummy)<-"Sample"
    dummy$Sample<-gsub("_.*","",eco.list[i])
  }else{
    dummy$Sample<-gsub("_.*","",eco.list[i])
  }
  if(!exists("out.eco")){
    out.eco<-dummy
  }else{
    
    missing.dummy<- colnames(out.eco)[-which(colnames(out.eco) %in% colnames(dummy) )]
    if(length(missing.dummy)>0){
      padd<-as.data.frame(matrix(NA, nrow = nrow(dummy), ncol = length(missing.dummy)))
      colnames(padd)<-missing.dummy
      dummy<-cbind(dummy, padd)
    }
    
    missing.eco<- colnames(dummy)[-which(colnames(dummy) %in% colnames(out.eco) )]
    if(length(missing.eco)>0){
      padd<-as.data.frame(matrix(NA, nrow = nrow(out.eco), ncol = length(missing.eco)))
      colnames(padd)<-missing.eco
      out.eco<-cbind(out.eco, padd)
    }
    
    out.eco<-rbind(out.eco,dummy)
    
  }
}
colnames(out.eco)[-which(colnames(out.eco)=="Sample")]<-paste("EcoPipeFastq:",colnames(out.eco)[-which(colnames(out.eco)=="Sample")],sep = "")

summ<-merge(summ, out.eco, by="Sample", all.x=TRUE)
rm(out.eco)

#Fasta 
eco.list<-list.files(pattern = "fasta_ecopipeline.csv")

if(exists("out.eco")) rm(out.eco)
for (i in 1:length(eco.list)) {
  dummy<-read.csv(eco.list[i])
  if(colnames(dummy)[1]=="NoEcol"){
    dummy<-as.data.frame(matrix(NA,nrow = 1, ncol = 1))
    colnames(dummy)<-"Sample"
    dummy$Sample<-gsub("_.*","",eco.list[i])
  }else{
    dummy$Sample<-gsub("_.*","",eco.list[i])
  }
  if(!exists("out.eco")){
    out.eco<-dummy
  }else{
    
    missing.dummy<- colnames(out.eco)[-which(colnames(out.eco) %in% colnames(dummy) )]
    if(length(missing.dummy)>0){
      padd<-as.data.frame(matrix(NA, nrow = nrow(dummy), ncol = length(missing.dummy)))
      colnames(padd)<-missing.dummy
      dummy<-cbind(dummy, padd)
    }
    
    missing.eco<- colnames(dummy)[-which(colnames(dummy) %in% colnames(out.eco) )]
    if(length(missing.eco)>0){
      padd<-as.data.frame(matrix(NA, nrow = nrow(out.eco), ncol = length(missing.eco)))
      colnames(padd)<-missing.eco
      out.eco<-cbind(out.eco, padd)
    }
    
    out.eco<-rbind(out.eco,dummy)
    
  }
}
colnames(out.eco)[-which(colnames(out.eco)=="Sample")]<-paste("EcoPipeAssemblies:",colnames(out.eco)[-which(colnames(out.eco)=="Sample")],sep = "")
summ<-merge(summ, out.eco, by="Sample", all.x=TRUE)
rm(out.eco)

# espk nhe1 ---------------------------------------------------------------

if(length(grep("EcoPipe", colnames(summ)))>0){
summ$Abricate_NleH1.1<-"Not detected"
summ$Abricate_NleH1.2<-"Not detected"

if(length(grep("nleH1", summ$Abricate_vfdb))>0) summ$Abricate_NleH1.1[grep("nleH1", summ$Abricate_vfdb)]<-"Detected"
if(length(grep("nleH2", summ$Abricate_vfdb))>0) summ$Abricate_NleH1.2[grep("nleH2", summ$Abricate_vfdb)]<-"Detected"

nleh<-which(is.na(summ[,c("EcoPipeAssemblies:Serotype")]))
summ$Abricate_espK<-"Not found assemblies"
if(length(grep("espK", summ$Abricate_vfdb))>0) summ$Abricate_espK[grep("espK", summ$Abricate_vfdb)]<-"Detected"
if(length(nleh)>0){
  summ$Abricate_NleH1.1[nleh]<-NA
  summ$Abricate_NleH1.2[nleh]<-NA
  summ$Abricate_espK[nleh]<-NA
}

summ$nleH1.1<-"Not Detected"
summ$nleH1.2<-"Not Detected"
#Code check
# nleh1index<- unique(c(which(summ$Abricate_NleH1.1=="Detected"), which(!is.na(summ$`EcoPipeFastq:nleH1.1`)),
#                       which(!is.na(summ$`EcoPipeAssemblies:nleH1.1`)))) 
# 
# nleh2index<- unique(c(which(summ$Abricate_NleH1.2=="Detected"), which(!is.na(summ$`EcoPipeFastq:nleH1.2`)),
#                       which(!is.na(summ$`EcoPipeAssemblies:nleH1.2`)))) 


nleh1index<- unique(c(which(summ$Abricate_NleH1.1=="Detected"), grep("^Detected", summ$`EcoPipeAssemblies:nleH1.1`), grep("^Detected", summ$`EcoPipeFastq:nleH1.1`)))
nleh2index<- unique(c(which(summ$Abricate_NleH1.2=="Detected"), grep("^Detected", summ$`EcoPipeAssemblies:nleH1.2`), grep("^Detected", summ$`EcoPipeFastq:nleH1.2`)))

if(length(nleh1index)>0)summ$nleH1.1[nleh1index] <-"Detected"
if(length(nleh2index)>0)summ$nleH1.2[nleh2index] <-"Detected"

non.eco<-which(is.na(summ$`EcoPipeFastq:nleH1.1`))

if(length(non.eco)>0){
  summ$nleH1.1[non.eco]<-NA
  summ$nleH1.2[non.eco]<-NA
} 

colnames(summ)[which(colnames(summ)=="nleH1.1")] <- "NleH1-1"
colnames(summ)[which(colnames(summ)=="nleH1.2")] <- "NleH1-2"
colnames(summ)[which(colnames(summ)=="Abricate_NleH1.1")] <- "Abricate:nleH1-1"
colnames(summ)[which(colnames(summ)=="Abricate_NleH1.2")] <- "Abricate:nleH1-2"
colnames(summ)[which(colnames(summ)=="Abricate_espK")] <- "Abricate:espK"
colnames(summ)<-gsub("nleH1\\.1", "nleH1-1", colnames(summ))  
}


# STXcolumns --------------------------------------------------------------

summ$Stx1<-NA
summ$Stx2<-NA

if(length(which(summ$StxType_Contigs=="No hit found" &  summ$StxType_Mapping=="No hit found"))>0){
  summ$Stx1[which(summ$StxType_Contigs=="No hit found" &  summ$StxType_Mapping=="No hit found")]<-"Assemblies:Non Detected | Reads:Non Detected"
  summ$Stx2[which(summ$StxType_Contigs=="No hit found" &  summ$StxType_Mapping=="No hit found")]<-"Assemblies:Non Detected | Reads:Non Detected"
}

if(length(which(summ$StxType_Contigs=="No hit found" &  summ$StxType_Mapping!="No hit found" & !is.na(summ$StxType_Mapping)))>0){
  summ$Stx1[which(summ$StxType_Contigs=="No hit found" &  summ$StxType_Mapping!="No hit found" & !is.na(summ$StxType_Mapping))]<-"Assemblies:Non Detected"
  summ$Stx2[which(summ$StxType_Contigs=="No hit found" &  summ$StxType_Mapping!="No hit found" & !is.na(summ$StxType_Mapping))]<-"Assemblies:Non Detected"
}

if(length(which(summ$StxType_Contigs!="No hit found" &  summ$StxType_Mapping=="No hit found" & !is.na(summ$StxIdentity_Contigs)))>0){
  summ$Stx1[which(summ$StxType_Contigs!="No hit found" &  summ$StxType_Mapping=="No hit found" & !is.na(summ$StxIdentity_Contigs))]<-"Reads:Non Detected"
  summ$Stx2[which(summ$StxType_Contigs!="No hit found" &  summ$StxType_Mapping=="No hit found" & !is.na(summ$StxIdentity_Contigs))]<-"Reads:Non Detected"
}

if(length(summ$StxVariant_Contigs)==0) summ$StxVariant_Contigs<-NA
if(length(summ$StxVariant_Mapping)==0) summ$StxVariant_Mapping<-NA

for (i in 1:nrow(summ)) {
  if(exists("stxvar"))rm(stxvar)
  if(exists("stxvar2"))rm(stxvar2)
  
  if(!is.na(summ$StxVariant_Contigs[i])){
    stxvar<- unique(unlist(strsplit(summ$StxVariant_Contigs[i], " \\| ")))
    stxvar<-gsub("\\:variant ","",stxvar)
    
    if(length(grep("stx1", stxvar))>0){
      summ$Stx1[i]<-paste("Assemblies:", paste(stxvar[grep("stx1", stxvar)], collapse = ";"),sep = "")
    }
    if(length(grep("stx2", stxvar))>0){
      summ$Stx2[i]<-paste("Assemblies:", paste(stxvar[grep("stx2", stxvar)], collapse = ";"),sep = "")
    }
    
  }
  
  if(!is.na(summ$StxVariant_Mapping[i])){
    stxvar2<- unique(unlist(strsplit(summ$StxVariant_Mapping[i], " \\| ")))
    stxvar2<-gsub("\\:variant ","",stxvar2)
    
    if(length(grep("stx1", stxvar2))>0){
      if(is.na(summ$Stx1[i])){
        summ$Stx1[i]<-paste("Reads:", paste(stxvar[grep("stx1", stxvar2)], collapse = ";"),sep = "")  
      }else{
        summ$Stx1[i]<-paste(summ$Stx1[i], paste("Reads:", paste(stxvar[grep("stx1", stxvar2)], collapse = ";"),sep = ""), sep = " | ")
      }
    }
    if(length(grep("stx2", stxvar2))>0){
      if(is.na(summ$Stx2[i])){
        summ$Stx2[i]<-paste("Reads:", paste(stxvar2[grep("stx2", stxvar2)], collapse = ";"),sep = "")  
      }else{
        summ$Stx2[i]<-paste(summ$Stx2[i], paste("Reads:", paste(stxvar2[grep("stx2", stxvar2)], collapse = ";"),sep = ""), sep = " | ")
      }
    }
  }
}


if(length(which(is.na(summ$Stx2) & !is.na(summ$Stx1)))>0) summ$Stx2[which(is.na(summ$Stx2) & !is.na(summ$Stx1))]<-"Assemblies:Non Detected | Reads:Non Detected"
if(length(which(is.na(summ$Stx1) & !is.na(summ$Stx2)))>0) summ$Stx1[which(is.na(summ$Stx1) & !is.na(summ$Stx2))]<-"Assemblies:Non Detected | Reads:Non Detected"

colnames(summ)[which(colnames(summ)=="Stx1")]<-"Stx-1"
colnames(summ)[which(colnames(summ)=="Stx2")]<-"Stx-2"



# SequencerID -------------------------------------------------------------

sqid.list<-list.files(pattern = "sequencerID.tsv")
if(exists("out.sqid")) rm(out.sqid)
for (i in 1:length(sqid.list)) {
  dum<-read.csv(sqid.list[i], sep = "\t", header = FALSE)
  
  dum2<-as.data.frame(t(c(gsub(":.*","",dum[1,1]), gsub("_sequencerID.tsv","",sqid.list[i]))))
  
  colnames(dum2)<-c("Instrument", "Sample")
  if(!exists("out.sqid")){
    out.sqid<-dum2
  }else{
    out.sqid<-rbind(out.sqid,dum2)
  }
}
out.sqid$Sample<-gsub("_.*","",out.sqid$Sample)
out.sqid$Instrument<-gsub("@","",out.sqid$Instrument)
summ<-merge(summ, out.sqid, by="Sample", all.x=TRUE)
rm(out.sqid)




# AMRFInderPlus -----------------------------------------------------------

amfplus<-list.files(pattern = "_amrfinderplus.tsv")
if(exists("outamr")) rm(outamr)

for (i in 1:length(amfplus)) {
  dummy<-read.csv(amfplus[i],sep = "\t")
  
  if(length(which(dummy$Element.type=="AMR"))>0){
  
  dummyamr<-dummy[which(dummy$Element.type=="AMR"),]
  if(nrow(dummyamr)>0){
    amrclasses<-unique(dummyamr$Class)
    amrpad<-as.data.frame(matrix(NA,nrow = 1, ncol = length(amrclasses)))
    colnames(amrpad)<-amrclasses
    for ( c in 1:ncol(amrpad)) {
      amrpad[,c]<-paste(dummyamr$Gene.symbol[which(dummyamr$Class==colnames(amrpad)[c])],collapse = " | ")
    }
    
  }
  amrpad$Sample<-unique(gsub("_.*","",dummy$Contig.id))
  
  if(!exists("outamr")){
    outamr<-amrpad
  }else{
    missing.out<-colnames(amrpad)[-which(colnames(amrpad) %in% colnames(outamr))]
    missing.dumm<-colnames(outamr)[-which(colnames(outamr) %in% colnames(amrpad))]
    
    if(length(missing.out)>0){
      padding<-as.data.frame(matrix(NA, nrow = nrow(outamr), ncol = length(missing.out)))
      colnames(padding)<-missing.out
      outamr<-cbind(outamr,padding)
    }
  
    if(length(missing.dumm)>0){
      padding<-as.data.frame(matrix(NA, nrow = nrow(amrpad), ncol = length(missing.dumm)))
      colnames(padding)<-missing.dumm
      amrpad<-cbind(amrpad,padding)
    }
    
    
    outamr<-rbind(outamr,amrpad)
  }
  }
}

if(exists("outamr")){
  colnames(outamr)[-which(colnames(outamr)=="Sample")]<-paste("AMRFinder_",colnames(outamr)[-which(colnames(outamr)=="Sample")],sep = "")
  summ<-merge(summ, outamr,by="Sample", all.x = TRUE)  
}



# BPEprofiler -------------------------------------------------------------
bpe<-list.files(pattern = "bpe_mlst.csv")
if(exists("outbpe")) rm(outbpe)

for (i in 1:length(bpe)) {
  dummy<-read.csv(bpe[i])
  if(colnames(dummy)[1]!="NoBper"){
    if(!exists("outbpe")){
      outbpe<-dummy
    }else{
      outbpe<-rbind(outbpe,dummy)
    }  
  }
}

if(exists("outbpe")){
outbpe$Sample<-gsub("_.*","",outbpe$Sample)
colnames(outbpe)[6]<-"cgMLST_Score"
colnames(outbpe)[7]<-"cgMLST_AlellesFound"
colnames(outbpe)[-1]<-paste("BPE_",colnames(outbpe)[-1],sep = "")

summ<-merge(summ, outbpe, by="Sample",all.x = TRUE, all.y = FALSE)

summ$MLST.Type[which(summ$Sample %in% outbpe$Sample)]<-
  summ$BPE_MLST.Type[which(summ$Sample %in% outbpe$Sample)]


summ$MLST.Scheme[which(summ$Sample %in% outbpe$Sample)]<-
  summ$BPE_MLST.Scheme[which(summ$Sample %in% outbpe$Sample)]


summ$MLST_Date[which(summ$Sample %in% outbpe$Sample)]<-
  summ$BPE_cgMLSTDate[which(summ$Sample %in% outbpe$Sample)]

if(length(which(colnames(summ)=="ClonalComplex"))>0){
  summ$ClonalComplex<-NA
}
summ$ClonalComplex[which(summ$Sample %in% outbpe$Sample)]<-
  summ$BPE_ClonalComplex[which(summ$Sample %in% outbpe$Sample)]

summ$BPE_ClonalComplex<-NULL
summ$BPE_MLST.Scheme<-NULL
summ$BPE_MLST.Type<-NULL

}

# Diphtoscan --------------------------------------------------------------

cds.files<-list.files(pattern = "diphtoscan.csv")
if(exists("outcds")) rm(outcds)

for (i in 1:length(cds.files)) {
  dummy<-read.csv(cds.files[i])
  if(colnames(dummy)[1]!="NoDiphto"){
    if(!exists("outcds")){
      outcds<-dummy
    }else{
      outcds<-rbind(outcds,dummy)
    }  
  }
}

if(exists("outcds")){

outcds$MLST.Scheme<-paste("atpA:", outcds$atpA, "|dnaE:", outcds$dnaE,"|dnaK:",outcds$dnaK, 
                                 "|fusA:", outcds$fusA,  "|leuA:",outcds$leuA ,  "|odhA:",outcds$odhA ,"|rpoB:", outcds$rpoB, sep = "")

colnames(outcds)[-which(colnames(outcds)=="Sample")]<-paste("DiphtoScan_",colnames(outcds)[-which(colnames(outcds)=="Sample")], sep = "")

summ<-merge(summ, outcds, by="Sample",all.x = TRUE, all.y = FALSE)

summ$MLST.Scheme[which(summ$Sample %in% outcds$Sample)]<-
  summ$DiphtoScan_MLST.Scheme[which(summ$Sample %in% outcds$Sample)]

summ$DiphtoScan_MLST.Scheme<-NULL

summ$MLST.Type[which(summ$Sample %in% outcds$Sample)]<-
  summ$DiphtoScan_ST[which(summ$Sample %in% outcds$Sample)]

summ$DiphtoScan_ST<-NULL

summ$MLST_Date[which(summ$Sample %in% outcds$Sample)]<-"DiphtoScan.11.11.2024"

summ$DiphtoScan_dnaK<-NULL
summ$DiphtoScan_atpA<-NULL
summ$DiphtoScan_dnaE<-NULL
summ$DiphtoScan_odhA<-NULL
summ$DiphtoScan_rpoB<-NULL
summ$DiphtoScan_leuA<-NULL
summ$DiphtoScan_fusA<-NULL

}



# TBProfiler --------------------------------------------------------------
tbp.files<-list.files(pattern = "tb_profiler.tsv")

if(exists("outtbpro")) rm(outtbpro)

for (i in 1:length(tbp.files)) {

  dummy<-read.csv(tbp.files[i],sep = "\t")
  if(colnames(dummy)[1]!="dummy"){
    if(!exists("outtbpro")){
      outtbpro<-dummy
    }else{
      outtbpro<-rbind(outtbpro,dummy)
    }
  }
  
}

if(exists("outtbpro")){
  outtbpro<-as.data.frame(apply(outtbpro,2,function(x) gsub("^DRUG_RESISTANCE / ", "", x)))
  summ<-merge(summ, outtbpro, by="Sample",all.x = TRUE, all.y = FALSE)
}



# versions ------------------------------------------------------------
versions<-read.csv("/home/docker/CommonFiles/Versions.csv", sep = ";", header = FALSE)
#versions<-read.csv("/media/nacho/Data/DockerImages/TOP_dev/Spades/CommonFiles/Versions.csv", sep = ";", header = FALSE)

summ$Software_versions<-paste(paste(versions$V1,versions$V2, sep = ":"), collapse = " | ")

# Last Stage --------------------------------------------------------------


empty.col<-apply(summ, 2, function(x) length(which(is.na(x))))
empty.col<-as.numeric(empty.col)
empty.col<-as.numeric(which(empty.col==nrow(summ)))

if(length(empty.col)>0) summ<-summ[,-empty.col]




write_xlsx(summ, paste("Summaries_",gsub("-","",Sys.Date()), ".xlsx",sep = ""))



#Implementation for LW
summ$MLST.Scheme<-paste("| ", summ$MLST.Scheme," |",sep = "")
to.separate<-apply(summ, 2, function(x)length(grep(",",x)))
to.separate<-which(to.separate>0)
if(length(to.separate)>0){
  for (col.to.sep in to.separate) {
    summ[,col.to.sep]<-gsub(","," |",summ[,col.to.sep])
  }
}

to.delete<-apply(summ, 2, function(x)length(which(x=="| NA |")))
to.delete<-which(to.delete>0)
if(length(to.delete)>0){
  for (col.to.del in to.delete) {
    summ[which(summ[,col.to.del]=="| NA |"),col.to.del]<-NA
  }
}
write.csv(summ, paste("Summaries_",gsub("-","",Sys.Date()), ".csv",sep = ""), row.names = FALSE)

