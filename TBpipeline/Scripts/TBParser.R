#library(readxl)
mykrobe.files<-list.files(pattern = "mykrobe_output.csv", recursive = TRUE, full.names = TRUE)[1]

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


  
#Dist
dist<-read.csv("TB_all_dists.csv")
colnames(dist)[1]<-"Sample"
singlesamp<-dist[which(dist$Sample==df.out$Sample),]
df.out$TB_SimilarSequences_5_SNP<-length(which(as.numeric(singlesamp[1,])<=5))-1
df.out$TB_SimilarSequences_12_SNP<-length(which(as.numeric(singlesamp[1,])<=12))-1

col.files<-list.files(pattern = "collytype.txt", recursive = TRUE, full.names = TRUE)
col.files<-col.files[1]
coll<-read.csv(col.files, sep="\t")
coll$Lineage<-paste("lineage", coll$Lineage,sep = "")
df.out$TB_CollType_Lineage<-paste(coll$Lineage,collapse = "/")



#QC

depth.files<-list.files(pattern = "averagedepth.txt", recursive = TRUE, full.names = TRUE)
depth.files<-depth.files[1]
depth<-read.csv(depth.files, header = FALSE)
df.out$TB_AverageDepth<-depth[1,1]
df.out$TBDB_Size<-length(list.files("/mnt/global_collection/"))

info.latex<-list.files(pattern = "info.tex", recursive = TRUE, full.names = TRUE)
info<-read.csv(info.latex,sep = "&", header = FALSE)
info.c1<-info[,c(1,2)]
colnames(info.c1)<-c("Item","Value")
info.c2<-info[,c(3,4)]
colnames(info.c2)<-c("Item","Value")
info<-rbind(info.c1, info.c2)

df.out$TB_Variants<-as.numeric(gsub(" ","",info$Value[grep("Variants", info$Item)]))
df.out$TB_QC<-gsub(" ","",info$Value[grep("Datakvalitet", info$Item)])
df.out$TB_PerctAligned<-as.numeric(gsub("\\\\.*","",gsub(" ","",info$Value[grep("aligned", info$Item)])))
df.out$TB_BasesLowCov<-as.numeric(gsub("\\\\.*","",gsub(" ","",info$Value[grep("low", info$Item)])))

#Snippy integration

# #Catalog integration
# catalog<-read_xlsx("/media/nacho/Data/DockerImages/TOP/TBpipeline/WHO-TB_Catalog-2021.7-eng.xlsx", sheet = 2)
# catalog<-read_xlsx("/media/tbuser/WHO-TB_Catalog-2021.7-eng.xlsx", sheet = 2)
# catalog<-catalog[-1,]
# vcf.list<-list.files(pattern = "snps.tab", recursive = TRUE, full.names = TRUE)
# vcf<-read.csv(vcf.list[1],sep = "\t")

write.csv(df.out, paste(df.out$Sample, "_tbp.csv", sep = ""), row.names = FALSE)




