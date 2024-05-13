#library(readxl)
mykrobe.files<-list.files("COPY_TO_TB_PIPELINE_DATABASE", pattern = "mykrobe_output.csv", recursive = TRUE, full.names = TRUE)

if(length(mykrobe.files)>0){
for (i in 1:length(mykrobe.files)) {
  dummy.myk<-read.csv(mykrobe.files[i], stringsAsFactors = FALSE)
  colnames(dummy.myk)[4]<-"Mutations"
  dummy.df<-as.data.frame(matrix(NA, nrow = 1, ncol = 7))
  colnames(dummy.df)<-c("Sample","TB_Myk_Phylogroup","TB_Myk_Species", "TB_Myk_Lineage", "TB_Myk_PredictedResistances", "TB_Myk_AMRMutations", "TB_Myk_Version")
  dummy.df$Sample<-unique(dummy.myk$sample)
  dummy.df$TB_Myk_Phylogroup<-unique(dummy.myk$phylo_group)
  dummy.df$TB_Myk_Species<-unique(dummy.myk$species)
  dummy.df$TB_Myk_Lineage<-unique(dummy.myk$lineage)
  dummy.df$TB_Myk_Version<-unique(dummy.myk$mykrobe_version)
  resistances<-vector()
  mutations<-vector()

  for (j in 1:nrow(dummy.myk)) {
    if(dummy.myk$susceptibility[j]=="R" | dummy.myk$susceptibility[j]=="r"){  
      resistances<- c(resistances, dummy.myk$drug[j])
      mut.dummy<-unlist(strsplit(dummy.myk$Mutations[j], ";"))
      gene.dummy<-gsub("-.*","", mut.dummy)
      ratios<-vector()
      for (k in 1:length(mut.dummy)) {
        kmers<-gsub(paste(":",gsub(".*:","",gsub(".*[A-Z]:", "", mut.dummy[k])),sep = ""), "", gsub(".*[A-Z]:", "", mut.dummy[k]))
        ratios[k]<-round(as.numeric(gsub(".*:","", kmers))/(as.numeric(gsub(":.*","", kmers))+as.numeric(gsub(".*:","", kmers)) ),3)
      }
      mutations.temp<-paste(paste(gene.dummy, rep(" [Ratio:", length(ratios)), ratios, rep("]", length(ratios)),sep = ""), collapse = "/")
      mutations<-c(mutations,mutations.temp)
     
    }
    if(dummy.myk$susceptibility[j]=="S"){
      resistances<- c(resistances, dummy.myk$drug[j])
      mutations<-c(mutations,"None")
    }
    
  }
  
  #dummy.df$TB_Myk_PredictedResistances<-paste( resistances, collapse = ";")
  #dummy.df$TB_Myk_AMRMutations<-paste(paste( resistances, mutations, sep = ":"), collapse = ";")
  dummy.df$TB_Myk_PredictedResistances<-NULL
  dummy.df$TB_Myk_AMRMutations<-NULL
  
  mutations<-mutations[order(resistances)]
  resistances<-resistances[order(resistances)]
  
  if(length(mutations)>0){
  add.drugs<-as.data.frame(t(mutations))
  colnames(add.drugs)<-paste("TB_",as.character(resistances),sep = "")
  dummy.df<-cbind(dummy.df, add.drugs)}else{
    add.drugs<-as.data.frame("Mykrobe Error")
    colnames(add.drugs)<-"TB_Warning"
    dummy.df<-cbind(dummy.df, add.drugs)
  }
  
  if(!exists("df.out")){
    df.out<-dummy.df
  }else{
    df.out<-rbind(df.out, dummy.df)
  }

}
}

  
#Dist
dist<-read.csv("TB_all_dists.csv")
colnames(dist)[1]<-"Sample"
df.out$TB_SimilarSequences_5_SNP<-NA
df.out$TB_SimilarSequences_12_SNP<-NA
for (i in 1:nrow(df.out)) {
  singlesamp<-dist[which(gsub("_merged","",dist$Sample) == df.out$Sample[i]),]
  df.out$TB_SimilarSequences_5_SNP[i]<-length(which(as.numeric(singlesamp[1,])<=5))
  df.out$TB_SimilarSequences_12_SNP[i]<-length(which(as.numeric(singlesamp[1,])<=12))
  
}

df.out$TB_SimilarSequences_5_SNP<-df.out$TB_SimilarSequences_5_SNP-1
df.out$TB_SimilarSequences_12_SNP<-df.out$TB_SimilarSequences_12_SNP-1


col.files<-list.files("COPY_TO_TB_PIPELINE_DATABASE/",pattern = "colltype.txt", recursive = TRUE, full.names = TRUE)
for (i in 1:length(col.files)) {
  dummycoll<-read.csv(col.files[i], sep="\t")
  dummycoll$Lineage<-paste("lineage ", paste(dummycoll$Lineage, collapse = "/"),sep = "")
  dummycoll<-dummycoll[1,]
  dummycoll$Sample<-gsub(".*/","",gsub("/colltype.txt","",col.files[i]))
  if(!exists("collout")){
    collout<-dummycoll
  }else{
    collout<-rbind(collout,dummycoll)
  }
}

collout$Sample<-gsub("_.*","", collout$Sample)
df.out$Sample<-gsub("_.*","",df.out$Sample)
colnames(collout)[which(colnames(collout)=="Lineage")]<-"TB_CollType_Lineage"
df.out<-merge(df.out, collout[,c("Sample","TB_CollType_Lineage")], by="Sample")


#QC

depth.files<-list.files("COPY_TO_TB_PIPELINE_DATABASE/",pattern = "averagedepth.txt", recursive = TRUE, full.names = TRUE)

for (i in 1:length(depth.files)) {
  depth<-read.csv(depth.files[i], header = FALSE)
  colnames(depth)<-"TB_AverageDepth"
  depth$Sample<-gsub(".*/","",gsub("/averagedepth.txt","",depth.files[i]))
  if(!exists("depth.out")){
    depth.out<-depth
  }else{
    depth.out<-rbind(depth.out, depth)
  }
}
depth.out$Sample<-gsub("_.*","", depth.out$Sample)

df.out<-merge(df.out, depth.out, by="Sample")

df.out$TBDB_Size<-length(list.files("/mnt/global_collection/"))

info.latex<-list.files(pattern = "info.tex", recursive = TRUE, full.names = TRUE)

for (i in 1:length(info.latex)) {
  info<-read.csv(info.latex[i],sep = "&", header = FALSE)  
  info.c1<-info[,c(1,2)]
  colnames(info.c1)<-c("Item","Value")
  info.c2<-info[,c(3,4)]
  colnames(info.c2)<-c("Item","Value")
  info<-rbind(info.c1, info.c2)
  dumm2<-as.data.frame(matrix(NA,ncol = 1, nrow = 1))
  colnames(dumm2)<-"TB_Variants"
  dumm2$TB_Variants<-as.numeric(gsub(" ","",info$Value[grep("Variants", info$Item)]))
  dumm2$TB_QC<-gsub(" ","",info$Value[grep("Datakvalitet", info$Item)])
  dumm2$TB_PerctAligned<-as.numeric(gsub("\\\\.*","",gsub(" ","",info$Value[grep("aligned", info$Item)])))
  dumm2$TB_BasesLowCov<-as.numeric(gsub("\\\\.*","",gsub(" ","",info$Value[grep("low", info$Item)])))
  dumm2$Sample<-gsub("/.*","",gsub("\\./","",info.latex[i]))
  if(!exists("latx.out")){
    latx.out<-dumm2
  }else{
    latx.out<-rbind(latx.out, dumm2)
  }
}

latx.out$Sample<-gsub("_.*","",latx.out$Sample)
df.out<-merge(df.out, latx.out, by="Sample")


write.csv(df.out, "summary_tbp.csv", row.names = FALSE)




